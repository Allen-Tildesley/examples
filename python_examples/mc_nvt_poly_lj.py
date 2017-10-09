#!/usr/bin/env python3
# mc_nvt_poly_lj.py

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

"""Monte Carlo, NVT ensemble, polyatomic molecule."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    # In this example we simulate using the shifted-force potential only
    # The values of < p_sf >, < e_sf > and density should be consistent (for this potential)
    # There are no long-range or delta corrections

    from averages_module import VariableType
    
    # Preliminary calculations
    vol = box**3  # Volume
    rho = n / vol # Density

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Move acceptance ratio
    m_r = VariableType ( nam = 'Move ratio', val = m_ratio, instant = False )

    # Internal energy per molecule (shifted-force potential)
    # Ideal gas contribution (assuming nonlinear molecules) plus total PE divided by N
    e_sf = VariableType ( nam = 'E/N shifted force', val = 3.0*temperature + total.pot/n )

    # Pressure (shifted-force potential)
    # Ideal gas contribution plus total virial divided by V
    p_sf = VariableType ( nam = 'P shifted force', val = rho*temperature + total.vir/vol )

    # Collect together into a list for averaging
    return [ m_r, e_sf, p_sf ]

# Takes in a configuration of polyatomic molecules (positions and quaternions)
# Cubic periodic boundary conditions
# Conducts Monte Carlo at the given temperature
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in
# However, input configuration, output configuration, most calculations, and all results
# are given in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1, epsilon = 1

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model (including the cutoff distance) is defined in mc_poly_lj_module

import json
import sys
import numpy as np
import math
from config_io_module  import read_cnf_mols, write_cnf_mols
from averages_module   import run_begin, run_end, blk_begin, blk_end, blk_add
from maths_module      import random_translate_vector, metropolis, random_rotate_quaternion, q_to_a
from mc_poly_lj_module import introduction, conclusion, potential, potential_1, PotentialType, na, db

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('mc_nvt_poly_lj')
print('Monte Carlo, constant-NVT ensemble, polyatomic molecule')
print('Simulation uses cut-and-shifted potential')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "temperature":1.0, "dr_max":0.05, "de_max":0.05}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
dr_max      = nml["dr_max"]      if "dr_max"      in nml else defaults["dr_max"]
de_max      = nml["de_max"]      if "de_max"      in nml else defaults["de_max"]

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)       )
print( "{:40}{:15.6f}".format('Specified temperature',     temperature) )
print( "{:40}{:15.6f}".format('Maximum r displacement',    dr_max)      )
print( "{:40}{:15.6f}".format('Maximum e displacement',    de_max)      )

# Read in initial configuration
n, box, r, e = read_cnf_mols ( cnf_prefix+inp_tag, quaternions=True )
print( "{:40}{:15d}  ".format('Number of particles', n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
print( "{:40}{:15.6f}".format('Density', n/box**3)  )
r = r / box           # Convert positions to box units
r = r - np.rint ( r ) # Periodic boundaries

# Calculate all bond vectors
d  = np.empty ( (n,na,3), dtype=np.float_ )
for i, ei in enumerate(e):
    ai = q_to_a ( ei ) # Rotation matrix for i
    d[i,:,:] = np.dot ( db, ai ) # NB: equivalent to ai_T*db, ai_T=transpose of ai

# Initial energy and overlap check
total = potential ( box, r, d )
assert not total.ovr, 'Overlap in initial configuration'

# Initialize arrays for averaging and write column headings
m_ratio = 0.0
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        moves = 0

        for i in range(n): # Loop over atoms
            rj = np.delete(r,i,0) # Array of all the other molecules
            dj = np.delete(d,i,0) # Array of all the other molecules
            partial_old = potential_1 ( r[i,:], d[i,:,:], box, rj, dj ) # Old potential, virial etc
            assert not partial_old.ovr, 'Overlap in current configuration'

            ri = random_translate_vector ( dr_max/box, r[i,:] ) # Trial move to new position (in box=1 units)
            ri = ri - np.rint ( ri )                            # Periodic boundary correction
            ei = random_rotate_quaternion ( de_max, e[i,:] )    # Trial rotation
            ai = q_to_a ( ei ) # Rotation matrix for i
            di = np.dot ( db, ai ) # NB: equivalent to ai_T*db, ai_T=transpose of ai
            partial_new = potential_1 ( ri, di, box, rj, dj ) # New atom potential, virial etc

            if not partial_new.ovr: # Test for non-overlapping configuration
                delta = partial_new.pot - partial_old.pot # Use cut (but not shifted) potential
                delta = delta / temperature
                if metropolis ( delta ): # Accept Metropolis test
                    total    = total + partial_new - partial_old # Update total values
                    r[i,:]   = ri                                # Update position
                    e[i,:]   = ei                                # Update quaternion
                    d[i,:,:] = di                                # Update bond vectors
                    moves    = moves + 1                         # Increment move counter

        m_ratio = moves / n

        blk_add ( calc_variables() )

    blk_end(blk)                                            # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'      # Number configuration by block
    write_cnf_mols ( cnf_prefix+sav_tag, n, box, r*box, e ) # Save configuration

run_end ( calc_variables() )

total = potential ( box, r, d ) # Double check book-keeping
assert not total.ovr, 'Overlap in final configuration'

write_cnf_mols ( cnf_prefix+out_tag, n, box, r*box, e ) # Save configuration
conclusion()
