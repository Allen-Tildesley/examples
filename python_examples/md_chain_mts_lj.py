#!/usr/bin/env python3
# md_chain_mts_lj.py

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

"""Molecular dynamics, multiple timesteps, chain molecule."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    import numpy as np
    import math
    from averages_module import msd, cke, VariableType

    # Preliminary calculations
    kin = 0.5*np.sum(v**2)
    rcm = np.sum ( r, axis=0 ) / n  # Centre of mass
    rsq = np.sum ( (r-rcm)**2 ) / n # Mean-squared distance from CM
    eng = kin+total.pot+total_spr   # Total energy

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Internal energy per atom
    # Total KE plus total LJ nonbonded energy plus total spring energy divided by N
    e_f = VariableType ( nam = 'E/N', val = eng/n )

    # Kinetic temperature
    # Remove 6 degrees of freedom for conserved linear and angular momentum
    t_k = VariableType ( nam = 'T kinetic', val = 2.0*kin/(3*n-6) )

    # Radius of gyration
    r_g = VariableType ( nam = 'Rg', val = math.sqrt(rsq) )

    # MSD of kinetic energy, intensive
    # Use special method to convert to Cv/N
    c_f = VariableType ( nam = 'Cv/N', val = kin/math.sqrt(n), method = cke, instant = False )

    # Mean-squared deviation of conserved energy per atom
    conserved_msd = VariableType ( nam = 'Conserved MSD', val = eng/n,
                                   method = msd, e_format = True, instant = False )

    # Collect together for averaging
    return [ e_f, t_k, r_g, c_f, conserved_msd ]

# Takes in a configuration of atoms in a linear chain (positions, velocities)
# NO periodic boundary conditions, no box
# Conducts molecular dynamics with springs and multiple timesteps
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Input configuration, output configuration, all calculations, and all results 
# are given in mass = 1 units, and in simulation units defined by the model 
# For example, for Lennard-Jones, sigma = 1, epsilon = 1

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model is defined in md_chain_lj_module.py

import json
import sys
import numpy as np
import math
from config_io_module   import read_cnf_atoms, write_cnf_atoms
from averages_module    import run_begin, run_end, blk_begin, blk_end, blk_add
from md_chain_lj_module import introduction, conclusion, zero_cm, force, spring, worst_bond, PotentialType

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('md_chain_mts_lj')
print('Molecular dynamics, constant-NVE ensemble, chain molecule, multiple time steps')
print('Particle mass=1 throughout')
print('No periodic boundaries')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":10000, "dt":0.0002, "k_spring":10000.0, "n_mts":10}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock   = nml["nblock"]   if "nblock"   in nml else defaults["nblock"]
nstep    = nml["nstep"]    if "nstep"    in nml else defaults["nstep"]
dt       = nml["dt"]       if "dt"       in nml else defaults["dt"]
k_spring = nml["k_spring"] if "k_spring" in nml else defaults["k_spring"]
n_mts    = nml["n_mts"]    if "n_mts"    in nml else defaults["n_mts"]

introduction()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)   )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)    )
print( "{:40}{:15.6f}".format('Time step',                 dt)       )
print( "{:40}{:15.6f}".format('Bond spring constant',      k_spring) )
print( "{:40}{:15d}  ".format('Multiple time step factor', n_mts)    )
print( "{:40}{:15.6f}".format('Large time step',           dt*n_mts) )

# Read in initial configuration
n, bond, r, v = read_cnf_atoms ( cnf_prefix+inp_tag, with_v=True)
print( "{:40}{:15d}  ".format('Number of particles',          n) )
print( "{:40}{:15.6f}".format('Bond length (in sigma units)', bond)  )
r, v = zero_cm ( r, v )
print( "{:40}{:15.6f}".format('Worst bond length deviation', worst_bond(bond,r) )  )

# Initial forces, potential, etc plus overlap check
total, f = force ( r )
assert not total.ovr, 'Overlap in initial configuration'
total_spr, g = spring ( k_spring, bond, r )

# Initialize arrays for averaging and write column headings
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        # Single time step of length n_mts*dt

        v = v + 0.5 * n_mts * dt * f  # Kick half-step (long) with nonbonded forces f

        for stp_mts in range(n_mts): # Loop over n_mts steps of length dt
            v = v + 0.5 * dt * g                        # Kick half-step (short) with spring forces g
            r = r + dt * v                              # Drift step (short)
            total_spr, g = spring ( k_spring, bond, r ) # Evaluate spring forces g and potential
            v = v + 0.5 * dt * g                        # Kick half-step (short) with spring forces g

        total, f = force ( r ) # Evaluate nonbonded forces f and potential
        assert not total.ovr, 'Overlap in configuration'
        v = v + 0.5 * n_mts * dt * f # Kick half-step (long) with nonbonded forces f

        blk_add ( calc_variables() )

    blk_end(blk)                                             # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'       # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, bond, r, v ) # Save configuration

run_end ( calc_variables() )

print( "{:40}{:15.6f}".format('Worst bond length deviation', worst_bond(bond,r) )  )
write_cnf_atoms ( cnf_prefix+out_tag, n, bond, r, v ) # Save configuration
conclusion()

