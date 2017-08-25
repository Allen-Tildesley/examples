#!/usr/bin/env python3
# mc_chain_nvt_sw.py

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

"""Monte Carlo, single chain, NVT, square wells."""
    
def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    import numpy as np
    from averages_module import msd, VariableType

    # Preliminary calculations
    rcm = np.sum ( r, axis=0 ) / n   # Centre of mass
    rsq = np.sum( (r - rcm)**2 ) / n # Mean-squared distance from CM

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Move acceptance ratios
    r_r = VariableType ( nam = 'Regrow ratio', val = r_ratio, instant = False )
    c_r = VariableType ( nam = 'Crank ratio',  val = c_ratio, instant = False )
    p_r = VariableType ( nam = 'Pivot ratio',  val = p_ratio, instant = False )

    # Total potential energy per atom
    # PE=negative of the number of square-well contacts divided by N
    e_x = VariableType ( nam = 'PE/N', val = -q/n )

    # Radius of gyration
    r_g = VariableType ( nam = 'Rg', val = np.sqrt(rsq) )

    # Heat Capacity per atom (excess, without ideal gas contribution, extensive)
    # MSD of PE / (sqrt(N)*T)
    # PE==negative of the number of square-well contacts
    c_x = VariableType ( nam = 'Cv(ex)/N', val = -q/(np.sqrt(n)*temperature),
                         method = msd, instant = False )

    # Collect together into a list for averaging
    return [ r_r, c_r, p_r, e_x, r_g, c_x ]

def update_histogram(q,q_min,q_max):
    """Updates the probability histogram of (negative) energies."""

    import numpy as np
    global h

    assert q >= 0 and q <= nq, 'q out of range'
    
    h[q]  = h[q] + 1.0
    q_min = min(q,q_min)
    q_max = max(q,q_max)
    return q_min, q_max

# Takes in a configuration of atom positions in a linear chain
# NO periodic boundary conditions, no box
# Conducts Monte Carlo, NVT ensemble using various moves
# such as CBMC regrowth, pivot, crankshaft
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Input configuration, output configuration, all calculations, and all results
# are given in simulation units defined by the model.
# atomic core diameter sigma = 1, well-depth epsilon=1
# Potential energy is -q, where q is the total number of attractive well interactions
# This is a negative integer, but for convenience we refer to q as energy
# Configurational weights are calculated on the basis of the athermal interactions
# only, i.e. weight=1 if no overlap, weight=0 if any overlap

import json
import sys
import numpy as np
from config_io_module   import read_cnf_atoms, write_cnf_atoms
from averages_module    import run_begin, run_end, blk_begin, blk_end, blk_add
from mc_chain_sw_module import introduction, conclusion, regrow, crank, pivot, qcount, weight

global h

cnf_prefix = 'cnf.'
his_prefix = 'his.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('mc_chain_nvt_sw')
print('Monte Carlo, constant-NVT ensemble, chain molecule, square wells')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":10000, "m_max":3, "k_max":32,
            "crank_max":0.5, "crank_fraction":0.5, "pivot_max":0.5, "pivot_fraction":0.5,
            "temperature":1.0, "q_range":1.5}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock         = nml["nblock"]         if "nblock"         in nml else defaults["nblock"]
nstep          = nml["nstep"]          if "nstep"          in nml else defaults["nstep"]
m_max          = nml["m_max"]          if "m_max"          in nml else defaults["m_max"]
k_max          = nml["k_max"]          if "k_max"          in nml else defaults["k_max"]
crank_max      = nml["crank_max"]      if "crank_max"      in nml else defaults["crank_max"]
crank_fraction = nml["crank_fraction"] if "crank_fraction" in nml else defaults["crank_fraction"]
pivot_max      = nml["pivot_max"]      if "pivot_max"      in nml else defaults["pivot_max"]
pivot_fraction = nml["pivot_fraction"] if "pivot_fraction" in nml else defaults["pivot_fraction"]
temperature    = nml["temperature"]    if "temperature"    in nml else defaults["temperature"]
q_range        = nml["q_range"]        if "q_range"        in nml else defaults["q_range"]

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',                nblock)         )
print( "{:40}{:15d}  ".format('Number of steps per block',       nstep)          )
print( "{:40}{:15d}  ".format('Max atoms in regrow',             m_max)          )
print( "{:40}{:15d}  ".format('Random tries per atom in regrow', k_max)          )
print( "{:40}{:15.6f}".format('Max move angle in crankshaft',    crank_max)      )
print( "{:40}{:15.6f}".format('Crank fraction',                  crank_fraction) )
print( "{:40}{:15.6f}".format('Max move angle in pivot',         pivot_max)      )
print( "{:40}{:15.6f}".format('Pivot fraction',                  pivot_fraction) )
print( "{:40}{:15.6f}".format('Temperature/well depth',          temperature)    )
print( "{:40}{:15.6f}".format('Attractive well range',           q_range)        )
if q_range<1.0:
    print('Warning, q_range < core diameter (1.0)')

# Read in initial configuration
n, bond, r = read_cnf_atoms ( cnf_prefix+inp_tag)
print( "{:40}{:15d}  ".format('Number of particles',          n)     )
print( "{:40}{:15.6f}".format('Bond length (in sigma units)', bond)  )
nq = 6*n

# Set number of crankshaft and pivot moves per step
n_crank = np.rint(crank_fraction*n).astype(np.int)
n_pivot = np.rint(pivot_fraction*n).astype(np.int)
print( "{:40}{:15d}".format('Number of crankshaft tries per step', n_crank) )
print( "{:40}{:15d}".format('Number of pivot tries per step',      n_pivot) )

# Calculate energy histogram and Boltzmann exponents, s(q), which are just values of E/kT
e = -np.arange(nq+1,dtype=np.float_)
e[0]=0.0
s = e /temperature
h = np.zeros_like(e)

# Initial energy calculation plus overlap check
assert weight ( r ) > 0, 'Overlap in initial configuration'
q = qcount ( r, q_range )
print( "{:40}{:15d}".format('Initial q', q)     )
q_min = q # Min q seen so far
q_max = q # Max q seen so far

# Initialize arrays for averaging and write column headings
r_ratio = 0.0
c_ratio = 0.0
p_ratio = 0.0
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        r, q, accepted = regrow ( s, m_max, k_max, bond, q_range, r, q )
        r_ratio = 1.0 if accepted else 0.0
        q_min, q_max = update_histogram(q,q_min,q_max)

        p_ratio = 0.0
        if n_pivot > 0:
            n_acc = 0
            for i_pivot in range(n_pivot):
                r, q, accepted = pivot ( s, pivot_max, q_range, r, q )
                if accepted:
                    n_acc = n_acc + 1
                q_min,q_max = update_histogram(q,q_min,q_max)
            p_ratio = n_acc / n_pivot

        c_ratio = 0.0
        if n_crank > 0:
            n_acc = 0
            for i_crank in range(n_crank):
                r, q, accepted = crank ( s, crank_max, q_range, r, q )
                if accepted:
                    n_acc = n_acc + 1
                q_min,q_max = update_histogram(q,q_min,q_max)
            c_ratio = n_acc / n_crank

        blk_add ( calc_variables() )

    blk_end(blk)                                       # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav' # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, bond, r ) # Save configuration

run_end ( calc_variables() )

write_cnf_atoms ( cnf_prefix+out_tag, n, bond, r ) # Save configuration
norm = np.sum(h)
h    = h / norm
np.savetxt(his_prefix+out_tag,np.column_stack((e[q_min:q_max+1],h[q_min:q_max+1])),fmt='%10.0f %20.8e')
conclusion()
