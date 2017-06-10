#!/usr/bin/env python3
# mc_chain_nvt_cbmc_lj.py

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

"""Monte Carlo, single chain, NVT, CBMC."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    import numpy as np
    import math
    from averages_module    import msd, VariableType
    from mc_chain_lj_module import potential, spring_pot, PotentialType

    # Preliminary calculations
    total = potential(r) # Nonbonded potential with overlap flag
    assert not total.ovr, 'Overlap in configuration'
    spr = spring_pot ( bond, k_spring, r ) # Total spring potential energy
    rcm = np.sum ( r, axis=0 ) / n         # Centre of mass
    rsq = np.sum( (r - rcm)**2 ) / n       # Mean-squared distance from CM

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Move acceptance ratio
    m_r = VariableType ( nam = 'Regrow ratio', val = m_ratio, instant = False )

    # Total potential energy per atom (excess, without ideal gas contribution)
    # Total PE of bond springs plus total LJ PE (not cut, nor shifted) divided by N
    e_x = VariableType ( nam = 'PE/N', val = (spr+total.pot)/n )

    # Radius of gyration
    r_g = VariableType ( nam = 'Rg', val = np.sqrt(rsq) )

    # Heat Capacity per atom (excess, without ideal gas contribution)
    # MSD of PE / (sqrt(N)*T)
    # Total PE of bond springs plus total LJ PE (not cut, nor shifted), divided by T
    c_x = VariableType ( nam = 'Cv(ex)/N', val = (spr+total.pot)/(np.sqrt(n)*temperature),
                         method = msd, instant = False )

    # Collect together into a list for averaging
    return [ m_r, e_x, r_g, c_x ]

# Takes in a configuration of atom positions in a linear chain
# NO periodic boundary conditions, no box
# Conducts Monte Carlo, NVT ensemble using CBMC regrowth moves
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Input configuration, output configuration, all calculations, and all results
# are given in simulation units defined by the model.
# E.g. for Lennard-Jones, atomic diameter sigma = 1, well-depth epsilon=1
# Configurational weights are calculated on the basis of the nonbonded interactions

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model is defined in mc_chain_lj_module

import json
import sys
import numpy as np
import math
from config_io_module   import read_cnf_atoms, write_cnf_atoms
from averages_module    import run_begin, run_end, blk_begin, blk_end, blk_add
from mc_chain_lj_module import introduction, conclusion, regrow

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('mc_chain_nvt_cbmc_lj')
print('Monte Carlo, constant-NVT ensemble, CBMC, chain molecule')
print('Simulation uses full nonbonded potential (no cutoff)')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":5000, "m_max":3, "k_max":32,
            "temperature":1.0, "k_spring":400.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
m_max       = nml["m_max"]       if "m_max"       in nml else defaults["m_max"]
k_max       = nml["k_max"]       if "k_max"       in nml else defaults["k_max"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
k_spring    = nml["k_spring"]    if "k_spring"    in nml else defaults["k_spring"]

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',                nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block',       nstep)       )
print( "{:40}{:15d}  ".format('Max atoms in regrow',             m_max)       )
print( "{:40}{:15d}  ".format('Random tries per atom in regrow', k_max)       )
print( "{:40}{:15.6f}".format('Specified temperature',           temperature) )
print( "{:40}{:15.6f}".format('Bond spring strength',            k_spring)    )

# Read in initial configuration
n, bond, r = read_cnf_atoms ( cnf_prefix+inp_tag)
print( "{:40}{:15d}  ".format('Number of particles',          n)     )
print( "{:40}{:15.6f}".format('Bond length (in sigma units)', bond)  )

# Initialize arrays for averaging and write column headings
m_ratio = 0.0
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        r, accepted = regrow ( temperature, m_max, k_max, bond, k_spring, r )
        m_ratio = 1.0 if accepted else 0.0

        blk_add ( calc_variables() )

    blk_end(blk)                                       # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav' # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, bond, r ) # Save configuration

run_end ( calc_variables() )

write_cnf_atoms ( cnf_prefix+out_tag, n, bond, r ) # Save configuration
conclusion()
