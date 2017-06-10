#!/usr/bin/env python3
# bd_nvt_lj.py

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

"""Brownian dynamics, NVT ensemble."""

def calc_variables ( ):
    """Calculates all variables of interest.
    
    They are collected and returned as a list, for use in the main program.
    """

    from averages_module import msd, VariableType
    from lrc_module import potential_lrc, pressure_lrc
    import numpy as np
    import math

    # Preliminary calculations (n,v,f,total are taken from the calling program)
    vol = box**3                  # Volume
    rho = n / vol                 # Density
    kin = 0.5*np.sum(v**2)        # Kinetic energy
    fsq = np.sum ( f**2 )         # Total squared force

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Internal energy (cut-and-shifted) per atom
    # Total KE plus total cut-and-shifted PE divided by N
    e_s = VariableType ( nam = 'E/N cut&shifted', val = (kin+total.pot)/n )

    # Internal energy (full, including LRC) per atom
    # LRC plus total KE plus total cut (but not shifted) PE divided by N
    e_f = VariableType ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + (kin+total.cut)/n )

    # Pressure (cut-and-shifted)
    # Ideal gas contribution plus total virial divided by V
    p_s = VariableType ( nam = 'P cut&shifted', val = rho*temperature + total.vir/vol )

    # Pressure (full, including LRC)
    # LRC plus ideal gas contribution plus total virial divided by V
    p_f = VariableType ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*temperature + total.vir/vol )

    # Kinetic temperature
    # Momentum is not conserved, hence 3N degrees of freedom
    t_k = VariableType ( nam = 'T kinetic', val = 2.0*kin/(3*n) )

    # Configurational temperature
    # Total squared force divided by total Laplacian
    t_c = VariableType ( nam = 'T config', val = fsq/total.lap )

    # Heat capacity (cut-and-shifted)
    # Total energy divided by temperature and sqrt(N) to make result intensive
    c_s = VariableType ( nam = 'Cv/N cut&shifted', val = (kin+total.pot)/(temperature*math.sqrt(n)),
                         method = msd, instant = False )

    # Heat capacity (full)
    # Total energy divided by temperature and sqrt(N) to make result intensive; LRC does not contribute
    c_f = VariableType ( nam = 'Cv/N full', val = (kin+total.cut)/(temperature*math.sqrt(n)),
                         method = msd, instant = False )

    # Collect together into a list for averaging
    return [ e_s, p_s, e_f, p_f, t_k, t_c, c_s, c_f ]

def a_propagator ( t ):
    """A: drift step propagator.

    t is the time over which to propagate (typically dt/2).
    r, v, and box are accessed from the calling program.
    """

    global r
    import numpy as np
    
    r = r + t * v / box   # Positions in box=1 units
    r = r - np.rint ( r ) # Periodic boundaries

def b_propagator ( t ):
    """B: kick step propagator.

    t is the time over which to propagate (typically dt/2).
    v is accessed from the calling program.
    """

    global v
    v = v + t * f

def o_propagator ( t ):
    """O: friction and random contributions propagator.

    t is the time over which to propagate (typically dt).
    v, n, temperature, and gamma are accessed from the calling program.
    """

    global v
    import numpy as np

    x = gamma*t
    c = 1-np.exp(-2*x) if x > 0.0001 else np.polyval([-2/3,4/3,-2.0,2.0,0.0],x)
    c = np.sqrt(c)
    v = v*np.exp(-x) + c*np.sqrt(temperature)*np.random.randn(n,3)
 
# Takes in a configuration of atoms (positions, velocities)
# Cubic periodic boundary conditions
# Conducts molecular dynamics using BAOAB algorithm of BJ Leimkuhler and C Matthews
# Appl. Math. Res. eXpress 2013, 34–56 (2013); J. Chem. Phys. 138, 174102 (2013)
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in and we assume mass=1 throughout
# However, input configuration, output configuration, most calculations, and all results 
# are given in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1, epsilon = 1

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model is defined in md_lj_module

import json
import sys
import numpy as np
import math
from config_io_module import read_cnf_atoms, write_cnf_atoms
from averages_module import run_begin, run_end, blk_begin, blk_end, blk_add
from md_lj_module import introduction, conclusion, force, PotentialType

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('bd_nvt_lj')
print('Brownian dynamics, constant-NVT ensemble')
print('Particle mass=1 throughout')
introduction()

np.random.seed()

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "r_cut":2.5, "dt":0.005, "temperature":1.0, "gamma":1.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
r_cut       = nml["r_cut"]       if "r_cut"       in nml else defaults["r_cut"]
dt          = nml["dt"]          if "dt"          in nml else defaults["dt"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
gamma       = nml["gamma"]       if "gamma"       in nml else defaults["gamma"]

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)            )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)             )
print( "{:40}{:15.6f}".format('Potential cutoff distance', r_cut)             )
print( "{:40}{:15.6f}".format('Time step',                 dt)                )
print( "{:40}{:15.6f}".format('Friction coefficient',      gamma)             )
print( "{:40}{:15.6f}".format('Specified temperature',     temperature)       )
print( "{:40}{:15.6f}".format('Ideal diffusion coefft',    temperature/gamma) )

# Read in initial configuration
n, box, r, v = read_cnf_atoms ( cnf_prefix+inp_tag, with_v=True)
print( "{:40}{:15d}  ".format('Number of particles',          n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
print( "{:40}{:15.6f}".format('Density', n/box**3)  )
r = r / box                    # Convert positions to box units
r = r - np.rint ( r )          # Periodic boundaries

# Initial forces, potential, etc plus overlap check
total, f = force ( box, r_cut, r )
assert not total.ovr, 'Overlap in initial configuration'

# Initialize arrays for averaging and write column headings
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        b_propagator ( dt/2 ) # B kick half-step
        a_propagator ( dt/2 ) # A drift half-step
        o_propagator ( dt )   # O random velocities and friction step
        a_propagator ( dt/2 ) # A drift half-step

        total, f = force ( box, r_cut, r ) # Force evaluation
        assert not total.ovr, 'Overlap in configuration'

        b_propagator ( dt/2 ) # B kick half-step

        blk_add ( calc_variables() )

    blk_end(blk)                                             # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'       # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r*box, v ) # Save configuration

run_end ( calc_variables() )

total, f = force ( box, r_cut, r ) # Force evaluation
assert not total.ovr, 'Overlap in final configuration'

write_cnf_atoms ( cnf_prefix+out_tag, n, box, r*box, v ) # Save configuration
conclusion()
