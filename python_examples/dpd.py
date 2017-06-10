#!/usr/bin/env python3
# dpd.py

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

"""Dissipative particle dynamics."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    # The DPD potential is short ranged, zero at, and beyond, r_cut
    # so issues of shifted potentials and long-range corrections do not arise

    from averages_module import VariableType
    import numpy as np
    import math

    # Preliminary calculations (n,r,v,f,total are taken from the calling program)
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

    # Kinetic temperature
    # Momentum is conserved, hence 3N-3 degrees of freedom
    t_k = VariableType ( nam = 'T kinetic', val = 2.0*kin/(3*n-3) )

    # Internal energy per atom
    # Total KE plus total PE divided by N
    e_f = VariableType ( nam = 'E/N', val = (kin+total.pot)/n )

    # Pressure
    # Ideal gas contribution plus total virial divided by V
    p_f = VariableType ( nam = 'P', val = rho*temperature + total.vir/vol )

    # Configurational temperature
    # Total squared force divided by total Laplacian
    t_c = VariableType ( nam = 'T config', val = fsq/total.lap )

    # Collect together into a list for averaging
    return [ e_f, t_k, t_c, p_f ]

def drift_propagator ( t ):
    """velocity Verlet drift step propagator.

    t is the time over which to propagate (typically dt).
    r, v, and box are accessed from the calling program.
    """

    global r
    import numpy as np

    r = r + t * v / box   # Positions in box=1 units
    r = r - np.rint ( r ) # Periodic boundaries

def kick_propagator ( t ):
    """velocity Verlet kick step propagator.

    t is the time over which to propagate (typically dt/2).
    v is accessed from the calling program.
    """

    global v
    v = v + t * f

# Takes in a configuration of atoms (positions, velocities)
# Cubic periodic boundary conditions
# Conducts dissipative particle dynamics using Shardlow or Lowe-Andersen algorithm
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in and we assume mass=1 throughout
# However, input configuration, output configuration, most calculations, and all results
# are given in simulation units defined by the model
# The range parameter (cutoff distance) is taken as unity

# The model is defined in dpd_module
# The typical DPD model described by Groot and Warren, J Chem Phys 107, 4423 (1997)
# has temperature kT=1, density rho=3, noise level sigma=3, gamma=sigma**2/(2*kT)=4.5
# and force strength parameter a=25 (more generally 75*kT/rho).
# We recommend a somewhat smaller timestep than their 0.04.
# They also give an approximate expression for the pressure, written out at the end for comparison

import json
import sys
import numpy as np
import math
from config_io_module import read_cnf_atoms, write_cnf_atoms
from averages_module  import run_begin, run_end, blk_begin, blk_end, blk_add
from dpd_module       import introduction, conclusion, force, lowe, shardlow, p_approx, PotentialType

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('dpd')
print('Dissipative particle dynamics, constant-NVT ensemble')
print('Particle mass=1 and cutoff=1 throughout')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "dt":0.02, "temperature":1.0, "a":75.0,
            "gamma":4.5, "method":"Lowe"}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
dt          = nml["dt"]          if "dt"          in nml else defaults["dt"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
a           = nml["a"]           if "a"           in nml else defaults["a"]
gamma       = nml["gamma"]       if "gamma"       in nml else defaults["gamma"]
method      = nml["method"]      if "method"      in nml else defaults["method"]

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',              nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block',     nstep)       )
print( "{:40}{:15.6f}".format('Time step',                     dt)          )
print( "{:40}{:15.6f}".format('Specified temperature',         temperature) )
print( "{:40}{:15.6f}".format('Force strength a*rho/kT',       a)           )
print( "{:40}{:15.6f}".format('Friction / thermal rate gamma', gamma)       )

method = method.lower()
assert "lowe" in method or "shardlow" in method, 'Unrecognized thermalization method'
if "shardlow" in method:
    thermalize=shardlow
    print('Shardlow integration method')
    print( "{:40}{:15.6f}".format('DPD sigma parameter', np.sqrt(2*gamma*temperature)) )
else:
    thermalize=lowe
    print('Lowe thermalization method')
    assert gamma*dt<1.0, 'gamma*dt too large'

# Read in initial configuration
n, box, r, v = read_cnf_atoms ( cnf_prefix+inp_tag, with_v=True)
print( "{:40}{:15d}  ".format('Number of particles',          n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
rho = n/box**3
a = a * temperature / rho # Scale force strength accordingly
print( "{:40}{:15.6f}".format('Density', rho)  )
r = r / box                    # Convert positions to box units
r = r - np.rint ( r )          # Periodic boundaries

# Initial forces, potential, etc plus overlap check
total, f, pairs = force ( box, a, r )

# Initialize arrays for averaging and write column headings
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        v = thermalize ( box, temperature, gamma*dt, v, pairs )
        kick_propagator ( dt/2 )
        drift_propagator ( dt )

        total, f, pairs = force ( box, a, r ) # Force evaluation

        kick_propagator ( dt/2 )

        blk_add ( calc_variables() )

    blk_end(blk)                                             # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'       # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r*box, v ) # Save configuration

run_end ( calc_variables() )

total, f, pairs = force ( box, a, r ) # Force evaluation

print( "{:40}{:15.6f}".format('Approx pressure', p_approx ( a, rho, temperature ) )  )
write_cnf_atoms ( cnf_prefix+out_tag, n, box, r*box, v ) # Save configuration
conclusion()
