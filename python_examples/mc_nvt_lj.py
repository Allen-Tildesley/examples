#!/usr/bin/env python3
# mc_nvt_lj.py

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

"""Monte Carlo, NVT ensemble."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    # In this example we simulate using the cut (but not shifted) potential
    # The values of < p_c >, < e_c > and density should be consistent (for this potential)
    # For comparison, long-range corrections are also applied to give
    # estimates of < e_f > and < p_f > for the full (uncut) potential
    # The value of the cut-and-shifted potential is not used, in this example

    import numpy as np
    import math
    from averages_module import msd, VariableType
    from lrc_module      import potential_lrc, pressure_lrc, pressure_delta
    from mc_lj_module    import force_sq
    
    # Preliminary calculations (n,r,total are taken from the calling program)
    vol = box**3                           # Volume
    rho = n / vol                          # Density
    fsq = force_sq ( box, r_cut, r ) # Total squared force

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Move acceptance ratio
    m_r = VariableType ( nam = 'Move ratio', val = m_ratio, instant = False )

    # Internal energy per atom for simulated, cut, potential
    # Ideal gas contribution plus cut (but not shifted) PE divided by N
    e_c = VariableType ( nam = 'E/N cut', val = 1.5*temperature + total.pot/n )

    # Internal energy per atom for full potential with LRC
    # LRC plus ideal gas contribution plus cut (but not shifted) PE divided by N
    e_f = VariableType ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + 1.5*temperature + total.pot/n )

    # Pressure for simulated, cut, potential
    # Delta correction plus ideal gas contribution plus total virial divided by V
    p_c = VariableType ( nam = 'P cut', val = pressure_delta(rho,r_cut) + rho*temperature + total.vir/vol )

    # Pressure for full potential with LRC
    # LRC plus ideal gas contribution plus total virial divided by V
    p_f = VariableType ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*temperature + total.vir/vol )

    # Configurational temperature
    # Total squared force divided by total Laplacian
    t_c = VariableType ( nam = 'T config', val = fsq/total.lap )

    # Heat capacity (full)
    # MSD potential energy divided by temperature and sqrt(N) to make result intensive; LRC does not contribute
    # We add ideal gas contribution, 1.5, afterwards
    c_f = VariableType ( nam = 'Cv/N full', val = total.pot/(temperature*math.sqrt(n)),
                         method = msd, add = 1.5, instant = False )

    # Collect together into a list for averaging
    return [ m_r, e_c, p_c, e_f, p_f, t_c, c_f ]

# Takes in a configuration of atoms (positions)
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
# The model is defined in mc_lj_module

import json
import sys
import numpy as np
import math
from config_io_module import read_cnf_atoms, write_cnf_atoms
from averages_module  import run_begin, run_end, blk_begin, blk_end, blk_add
from maths_module     import random_translate_vector, metropolis
from mc_lj_module     import introduction, conclusion, potential, potential_1, PotentialType

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('mc_nvt_lj')
print('Monte Carlo, constant-NVT ensemble')
print('Simulation uses cut (but not shifted) potential')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "temperature":1.0, "r_cut":2.5, "dr_max":0.15}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
r_cut       = nml["r_cut"]       if "r_cut"       in nml else defaults["r_cut"]
dr_max      = nml["dr_max"]      if "dr_max"      in nml else defaults["dr_max"]

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)       )
print( "{:40}{:15.6f}".format('Specified temperature',     temperature) )
print( "{:40}{:15.6f}".format('Potential cutoff distance', r_cut)       )
print( "{:40}{:15.6f}".format('Maximum displacement',      dr_max)      )

# Read in initial configuration
n, box, r = read_cnf_atoms ( cnf_prefix+inp_tag)
print( "{:40}{:15d}  ".format('Number of particles',          n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
print( "{:40}{:15.6f}".format('Density', n/box**3)  )
r = r / box           # Convert positions to box units
r = r - np.rint ( r ) # Periodic boundaries

# Initial energy and overlap check
total = potential ( box, r_cut, r )
assert not total.ovr, 'Overlap in initial configuration'

# Initialize arrays for averaging and write column headings
m_ratio = 0.0
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        moves = 0

        for i in range(n): # Loop over atoms
            rj = np.delete(r,i,0) # Array of all the other atoms
            partial_old = potential_1 ( r[i,:], box, r_cut, rj ) # Old atom potential, virial etc
            assert not partial_old.ovr, 'Overlap in current configuration'

            ri = random_translate_vector ( dr_max/box, r[i,:] ) # Trial move to new position (in box=1 units)
            ri = ri - np.rint ( ri )                            # Periodic boundary correction
            partial_new = potential_1 ( ri, box, r_cut, rj )    # New atom potential, virial etc

            if not partial_new.ovr: # Test for non-overlapping configuration
                delta = partial_new.pot - partial_old.pot # Use cut (but not shifted) potential
                delta = delta / temperature
                if metropolis ( delta ): # Accept Metropolis test
                    total = total + partial_new - partial_old # Update total values
                    r[i,:] = ri                               # Update position
                    moves = moves + 1                         # Increment move counter

        m_ratio = moves / n

        blk_add ( calc_variables() )

    blk_end(blk)                                          # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'    # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r*box ) # Save configuration

run_end ( calc_variables() )

total = potential ( box, r_cut, r ) # Double check book-keeping
assert not total.ovr, 'Overlap in final configuration'

write_cnf_atoms ( cnf_prefix+out_tag, n, box, r*box ) # Save configuration
conclusion()
