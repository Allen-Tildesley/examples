#!/usr/bin/env python3
# md_nvt_lj_le.py

#------------------------------------------------------------------------------------------------#
# This software was written in 2016/17                                                           #
# by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        #
# and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             #
# to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     #
# published by Oxford University Press ("the publishers").                                       #
#                                                                                                #
# LICENCE                                                                                        #
# Creative Commons CC0 Public Domain Dedication.                                                 #
# To the extent possible under law, the authors have dedicated all copyright and related         #
# and neighboring rights to this software to the PUBLIC domain worldwide.                        #
# This software is distributed without any warranty.                                             #
# You should have received a copy of the CC0 Public Domain Dedication along with this software.  #
# If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               #
#                                                                                                #
# DISCLAIMER                                                                                     #
# The authors and publishers make no warranties about the software, and disclaim liability       #
# for all uses of the software, to the fullest extent permitted by applicable law.               #
# The authors and publishers do not recommend use of this software for any purpose.              #
# It is made freely available, solely to clarify points made in the text. When using or citing   #
# the software, you should not imply endorsement by the authors or publishers.                   #
#------------------------------------------------------------------------------------------------#

"""Molecular dynamics, NVT ensemble, Lees-Edwards boundaries."""

def calc_variables ( ):
    """Calculates all variables of interest.
    
    They are collected and returned as a list, for use in the main program.
    """

    import numpy as np
    import math
    from averages_module import msd, VariableType

    # Preliminary calculations (n,r,v,f,total are taken from the calling program)
    vol = box**3                      # Volume
    rho = n / vol                     # Density
    kin = 0.5*np.sum(v**2)            # Kinetic energy
    fsq = np.sum ( f**2 )             # Total squared force
    tmp = 2.0 * kin / (3*n-3)         # Remove three degrees of freedom for momentum conservation
    kyx = np.sum(v[:,0]*v[:,1]) / vol # Kinetic part of off-diagonal pressure tensor
    eng = kin + total.pot             # Total energy

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Internal energy per atom
    # Total KE plus total PE divided by N
    e_s = VariableType ( nam = 'E/N', val = eng/n )

    # Pressure
    # Ideal gas contribution plus total virial divided by V
    p_s = VariableType ( nam = 'P', val = rho*tmp + total.vir/vol )

    # Kinetic temperature
    t_k = VariableType ( nam = 'T kinetic', val = tmp )

    # Configurational temperature
    # Total squared force divided by total Laplacian
    t_c = VariableType ( nam = 'T config', val = fsq/total.lap )

    # Shear viscosity
    if np.fabs(strain_rate)<tol: # Guard against simulation with zero strain rate
       eta = VariableType ( nam = 'Shear viscosity', val = 0.0 )
    else:
       eta = VariableType ( nam = 'Shear viscosity', val = -(kyx+total.pyx/vol) / strain_rate )

    # MSD of conserved kinetic energy
    conserved_msd = VariableType ( nam = 'Conserved MSD', val = kin/n,
                                   method = msd, e_format = True, instant = False )

    # Collect together into a list for averaging
    return [ e_s, p_s, t_k, t_c, eta, conserved_msd ]

def a_propagator ( t ):
    """A propagator. t is the time over which to propagate (typically dt/2)."""

    global r, strain
    import numpy as np

    x = t * strain_rate # Change in strain (dimensionless)

    r[:,0] = r[:,0] + x * r[:,1]         # Extra strain term
    r      = r      + t * v / box        # Drift half-step (positions in box=1 units)
    strain = strain + x                  # Advance strain and hence boundaries
    strain = strain - np.rint ( strain ) # Keep strain within (-0.5,0.5)

    r[:,0] = r[:,0] - np.rint ( r[:,1] ) * strain # Extra PBC correction (box=1 units)
    r = r - np.rint ( r )                         # Periodic boundaries (box=1 units)

def b1_propagator ( t ): 
    """B1 propagator. t is the time over which to propagate (typically dt/2)."""

    global v
    import numpy as np

    x = t * strain_rate # Change in strain (dimensionless)

    c1 = x * np.sum ( v[:,0]*v[:,1] ) / np.sum ( v**2 )
    c2 = ( x**2 ) * np.sum ( v[:,1]**2 ) / np.sum ( v**2 )
    g  = 1.0 / np.sqrt ( 1.0 - 2.0*c1 + c2 )

    v[:,0] = v[:,0] - x*v[:,1]
    v      = g * v

def b2_propagator ( t ):
    """B2 propagator. t is the time over which to propagate (typically dt)."""

    global v
    import numpy as np

    alpha = np.sum ( f*v ) / np.sum ( v**2 )
    beta  = np.sqrt ( np.sum ( f**2 ) / np.sum ( v**2 ) )
    h     = ( alpha + beta ) / ( alpha - beta )
    e     = np.exp ( -beta * t )

    dt_factor = ( 1 + h - e - h / e ) / ( ( 1 - h ) * beta )
    prefactor = ( 1 - h ) / ( e - h / e )

    v = prefactor * ( v + dt_factor * f )
  
# Takes in a configuration of atoms (positions, velocities)
# Cubic periodic boundary conditions, with Lees-Edwards shear
# Conducts molecular dynamics, SLLOD algorithm, with isokinetic thermostat
# Refs: Pan et al J Chem Phys 122 094114 (2005)

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in and we assume mass=1 throughout
# However, input configuration, output configuration, most calculations, and all results 
# are given in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1, epsilon = 1

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model is defined in md_lj_le_module

import json
import sys
import numpy as np
import math
from config_io_module import read_cnf_atoms, write_cnf_atoms
from averages_module  import run_begin, run_end, blk_begin, blk_end, blk_add
from md_lj_le_module  import introduction, conclusion, force, PotentialType

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'
tol = 1.0e-6

print('md_nvt_lj_le')
print('Molecular dynamics, constant-NVT ensemble, Lees-Edwards')
print('Particle mass=1 throughout')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":10000, "dt":0.005, "strain_rate":0.04}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
dt          = nml["dt"]          if "dt"          in nml else defaults["dt"]
strain_rate = nml["strain_rate"] if "strain_rate" in nml else defaults["strain_rate"]

introduction()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)       )
print( "{:40}{:15.6f}".format('Time step',                 dt)          )
print( "{:40}{:15.6f}".format('Strain rate',               strain_rate) )

# Insist that strain be zero (i.e. an integer) at end of each block
strain = strain_rate * dt * nstep
strain = strain - np.rint ( strain )
assert np.fabs(strain) < tol, 'Strain must be zero at end of block'

# Read in initial configuration
n, box, r, v = read_cnf_atoms ( cnf_prefix+inp_tag, with_v=True)
print( "{:40}{:15d}  ".format('Number of particles',          n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
print( "{:40}{:15.6f}".format('Density', n/box**3)  )
strain = 0.0                                  # Assume for simplicity that this is true
r      = r / box                              # Convert positions to box units
r[:,0] = r[:,0] - np.rint ( r[:,1] ) * strain # Extra correction (box=1 units)
r      = r - np.rint ( r )                    # Periodic boundaries
vcm    = np.sum ( v, axis=0 ) / n             # Centre-of mass velocity
v      = v - vcm                              # Set COM velocity to zero

# Initial forces, potential, etc plus overlap check
total, f = force ( box, strain, r )
assert not total.ovr, 'Overlap in initial configuration'

# Initialize arrays for averaging and write column headings
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        # Isokinetic SLLOD algorithm (Pan et al)

        a_propagator  ( dt/2 )
        b1_propagator ( dt/2 )

        total, f = force ( box, strain, r ) # Force evaluation
        assert not total.ovr, 'Overlap in configuration'

        b2_propagator ( dt )
        b1_propagator ( dt/2 )
        a_propagator  ( dt/2 )

        blk_add ( calc_variables() )

    blk_end(blk)                                             # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'       # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r*box, v ) # Save configuration

run_end ( calc_variables() )

total, f = force ( box, strain, r ) # Force evaluation
assert not total.ovr, 'Overlap in final configuration'

write_cnf_atoms ( cnf_prefix+out_tag, n, box, r*box, v ) # Save configuration
conclusion()
