#!/usr/bin/env python3
# mc_npt_hs.py

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

"""Monte Carlo, NPT ensemble, hard spheres."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    from averages_module import VariableType
    import numpy as np
    import math

    # Preliminary calculations (m_ratio, v_ratio, box are taken from the calling program)
    vol = box**3  # Volume
    rho = n / vol # Density

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Move acceptance ratios

    m_r = VariableType ( nam = 'Move ratio',   val = m_ratio, instant = False )
    v_r = VariableType ( nam = 'Volume ratio', val = v_ratio, instant = False )

    # Density
    density = VariableType ( nam = 'Density', val = rho )

    # Collect together into a list for averaging
    return [ m_r, v_r, density ]

# Takes in a configuration of atoms (positions)
# Cubic periodic boundary conditions
# Conducts Monte Carlo at given NPT (the temperature is irrelevant)
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# We take kT=1 throughout defining the unit of energy
# Positions r are divided by box length after reading in
# However, input configuration, output configuration, most calculations, and all results
# are given in simulation units defined by the model
# in this case, for hard spheres, sigma = 1

# The logarithm of the box length is sampled uniformly

import json
import sys
import numpy as np
import math
from config_io_module  import read_cnf_atoms, write_cnf_atoms
from averages_module   import run_begin, run_end, blk_begin, blk_end, blk_add
from maths_module      import random_translate_vector, metropolis
from mc_hs_module      import introduction, conclusion, overlap, overlap_1

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('mc_npt_hs')
print('Monte Carlo, constant-NPT ensemble')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "dr_max":0.15, "db_max":0.005, "pressure":4.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock   = nml["nblock"]   if "nblock"   in nml else defaults["nblock"]
nstep    = nml["nstep"]    if "nstep"    in nml else defaults["nstep"]
dr_max   = nml["dr_max"]   if "dr_max"   in nml else defaults["dr_max"]
db_max   = nml["db_max"]   if "db_max"   in nml else defaults["db_max"]
pressure = nml["pressure"] if "pressure" in nml else defaults["pressure"]

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)   )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)    )
print( "{:40}{:15.6f}".format('Pressure',                  pressure) )
print( "{:40}{:15.6f}".format('Maximum displacement',      dr_max)   )
print( "{:40}{:15.6f}".format('Maximum box displacement',  db_max)   )

# Read in initial configuration
n, box, r = read_cnf_atoms ( cnf_prefix+inp_tag)
print( "{:40}{:15d}  ".format('Number of particles',          n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
r = r / box           # Convert positions to box units
r = r - np.rint ( r ) # Periodic boundaries

# Initial pressure and overlap check
assert not overlap ( box, r ), 'Overlap in initial configuration'

# Initialize arrays for averaging and write column headings
m_ratio = 0.0
v_ratio = 0.0
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        moves = 0

        for i in range(n): # Loop over atoms
            ri = random_translate_vector ( dr_max/box, r[i,:] ) # Trial move to new position (in box=1 units)
            ri = ri - np.rint ( ri )                            # Periodic boundary correction
            rj = np.delete(r,i,0)                               # Array of all the other atoms

            if not overlap_1 ( ri, box, rj ): # Test for non-overlapping configuration
                r[i,:] = ri                         # Update position
                moves = moves + 1                   # Increment move counter

        m_ratio = moves / n

        v_ratio   = 0.0                 # Zero volume move counter
        zeta      = np.random.rand()    # Uniform random number in range (0,1)
        zeta      = 2.0*zeta-1.0        # Now in range (-1,+1)
        box_scale = np.exp(zeta*db_max) # Sampling log(box) and log(vol) uniformly
        box_new   = box*box_scale       # New box (in sigma units)
        den_scale = 1.0 / box_scale**3  # Density scaling factor

        if not overlap ( box_new, r ):           # Test for non-overlapping configuration
            delta = pressure * ( box_new**3 - box**3 ) # PV term (temperature=1.0)
            delta = delta + (n+1) * np.log(den_scale)  # Factor (n+1) consistent with log(box) sampling

            if metropolis(delta): # Accept Metropolis test
                box     = box_new   # Update box
                v_ratio = 1.0       # Set volume move counter

        blk_add ( calc_variables() )

    blk_end(blk)                                          # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'    # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r*box ) # Save configuration

run_end ( calc_variables() )

assert not overlap ( box, r ), 'Overlap in final configuration'

write_cnf_atoms ( cnf_prefix+out_tag, n, box, r*box ) # Save configuration
conclusion()
