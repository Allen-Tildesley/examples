#!/usr/bin/env python3
# smc_nvt_lj.py

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

"""Smart Monte Carlo, NVT ensemble."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """
    
    import numpy as np
    import math
    from averages_module import msd, VariableType
    from lrc_module      import potential_lrc, pressure_lrc
    
    # Preliminary calculations (n,r,total are taken from the calling program)
    vol = box**3       # Volume
    rho = n / vol      # Density
    fsq = np.sum(f**2) # Total squared force

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Move acceptance ratio

    m_r = VariableType ( nam = 'Move ratio', val = m_ratio, instant = False )

    # Internal energy per atom for simulated, cut-and-shifted, potential
    # Ideal gas contribution plus total cut-and-shifted PE divided by N
    e_s = VariableType ( nam = 'E/N cut&shifted', val = 1.5*temperature + total.pot/n )

    # Internal energy per atom for full potential with LRC
    # LRC plus ideal gas contribution plus total cut (but not shifted) PE divided by N
    e_f = VariableType ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + 1.5*temperature + total.cut/n )

    # Pressure for simulated, cut-and-shifted, potential
    # Ideal gas contribution plus total virial divided by V
    p_s = VariableType ( nam = 'P cut&shifted', val = rho*temperature + total.vir/vol )

    # Pressure for full potential with LRC
    # LRC plus ideal gas contribution plus total virial divided by V
    p_f = VariableType ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*temperature + total.vir/vol )

    # Configurational temperature
    # Total squared force divided by total Laplacian
    t_c = VariableType ( nam = 'T config', val = fsq/total.lap )

    # Heat capacity (excess, cut-and-shifted)
    # Total PE divided by temperature and sqrt(N) to make result intensive
    # We add ideal gas contribution, 1.5, afterwards
    c_s = VariableType ( nam = 'Cv/N cut&shifted', val = total.pot/(temperature*math.sqrt(n)),
                         method = msd, add = 1.5, instant = False )

    # Heat capacity (excess, full)
    # Total PE divided by temperature and sqrt(N) to make result intensive; LRC does not contribute
    # We add ideal gas contribution, 1.5, afterwards
    c_f = VariableType ( nam = 'Cv/N full', val = total.cut/(temperature*math.sqrt(n)),
                         method = msd, add = 1.5, instant = False )

    # Collect together into a list for averaging
    return [ m_r, e_s, p_s, e_f, p_f, t_c, c_s, c_f ]

# Takes in a configuration of atoms (positions)
# Cubic periodic boundary conditions
# Conducts Smart Monte Carlo using Hybrid Monte Carlo / Brownian Dynamics notation
# Uses no special neighbour lists
# Assume that a sweep consists of either
# (a) N successive single-particle moves
# (b) 1 multi-particle move involving a large fraction of atoms
# (large enough to justify calling the complete force routine)
# The ensemble corresponds to the shifted potential, not the simple cutoff potential

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in, and we assume mass=1 throughout
# However, input configuration, output configuration, most calculations, and all results 
# are given in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1, epsilon = 1

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model is defined in smc_module

import json
import sys
import numpy as np
import math
from config_io_module import read_cnf_atoms, write_cnf_atoms
from averages_module  import run_begin, run_end, blk_begin, blk_end, blk_add
from maths_module     import metropolis
from smc_lj_module    import introduction, conclusion, force, force_1, PotentialType

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('smc_nvt_lj')
print('Smart Monte Carlo, constant-NVT ensemble')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "temperature":1.0, "r_cut":2.5, "dt":0.1, "single_atom":True, "fraction":1.0 }
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
dt          = nml["dt"]          if "dt"          in nml else defaults["dt"]
single_atom = nml["single_atom"] if "single_atom" in nml else defaults["single_atom"]
fraction    = nml["fraction"]    if "fraction"    in nml else defaults["fraction"]

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)       )
print( "{:40}{:15.6f}".format('Specified temperature',     temperature) )
print( "{:40}{:15.6f}".format('Potential cutoff distance', r_cut)       )
print( "{:40}{:15.6f}".format('Time step',                 dt)          )
if single_atom:
    print( 'Single-atom moves' )
else:
    print( "{:40}{:15.6f}".format('Multi-atom moves with fraction',fraction) )
    assert 0.0 <= fraction <= 1.0, "Error: fraction out of range"

introduction()
np.random.seed()

v_rms = np.sqrt ( temperature ) # RMS value for velocity selection
print( "{:40}{:15.6f}".format('Typical dr', v_rms*dt) )

# Read in initial configuration
n, box, r = read_cnf_atoms ( cnf_prefix+inp_tag)
print( "{:40}{:15d}  ".format('Number of particles',n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
print( "{:40}{:15.6f}".format('Density', n/box**3)  )
r = r / box           # Convert positions to box units
r = r - np.rint ( r ) # Periodic boundaries

# Initial energy and overlap check
total, f = force ( box, r_cut, r )
assert not total.ovr, 'Overlap in initial configuration'

# Initialize arrays for averaging and write column headings
m_ratio = 0.0
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        if single_atom: # Single-atom moves
            
            n_move = 0

            for i in range(n): # Loop over atoms
                r_old = r[i,:].copy() # Store old position of this atom
                rj = np.delete(r,i,0) # Array of all the other atoms
                partial_old, f_old = force_1 ( r[i,:], box, r_cut, rj ) # Old forces, pot etc
                assert not partial_old.ovr, 'Overlap in current configuration'
                v                  = np.random.randn(3)*v_rms             # Choose 3 random momentum components
                kin_old            = 0.5*np.sum(v**2)                     # Old kinetic energy of this atom
                v                  = v + 0.5 * dt * np.sum(f_old,axis=0)  # Kick half-step for one atom with old force
                r[i,:]             = r[i,:] + dt * v / box                # Drift step (positions in box=1 units)
                r[i,:]             = r[i,:] - np.rint ( r[i,:] )          # Periodic boundaries (box=1 units)
                partial_new, f_new = force_1 ( r[i,:], box, r_cut, rj )   # New forces and pot etc for this atom

                if partial_new.ovr: # Test for overlap
                    r[i,:] = r_old # Restore position: this move is rejected
                else:
                    v  = v + 0.5 * dt * np.sum(f_new,axis=0) # Kick half-step for one atom with new force
                    kin_new = 0.5*np.sum(v**2)               # New kinetic energy of this atom

                    delta = partial_new.pot - partial_old.pot # Cut-and-shifted potential
                    delta = delta + kin_new - kin_old         # Include kinetic energy change
                    delta = delta / temperature               # Divide by temperature

                    if metropolis ( delta ): # Accept Metropolis test
                        total = total + partial_new - partial_old         # Update total values
                        f[:i,:]   = f[:i,:] - f_new[:i,:] + f_old[:i,:]   # change in forces due to i on other atoms j<i
                        f[i,:]    = f[i,:] + np.sum(f_new,axis=0) - np.sum(f_old,axis=0)
                        f[i+1:,:] = f[i+1:,:] - f_new[i:,:] + f_old[i:,:] # change in forces due to i on other atoms j>i
                        n_move    = n_move + 1                            # Update move counter
                    else:
                        r[i,:] = r_old # Restore position: this move is rejected

            m_ratio = n_move / n

        else: # Multi-atom moves

            move      = np.random.rand(n) < fraction         # Construct mask for moving atoms
            r_old     = r.copy()                             # Store old positions
            total_old = total                                # Store old totals
            f_old     = f                                    # Store old forces
            v         = np.random.randn(n,3)*v_rms           # Choose 3*n random momenta
            kin_old   = 0.5*np.sum(v**2)                     # Old kinetic energy
            v[move,:] = v[move,:] + 0.5 * dt * f_old[move,:] # Kick half-step with old forces
            r[move,:] = r[move,:] + dt * v[move,:] / box     # Drift step (positions in box=1 units)
            r[move,:] = r[move,:] - np.rint ( r[move,:] )    # Periodic boundaries (box=1 units)
            total, f  = force ( box, r_cut, r )              # New force and potential etc

            if total.ovr: # Test for overlap
              r       = r_old.copy() # Restore positions: this move is rejected
              total   = total_old    # Restore old totals
              f       = f_old        # Restore old forces
              m_ratio = 0.0          # Set move counter
            else:
                v[move,:] = v[move,:] + 0.5 * dt * f[move,:] # Kick half-step with new forces
                kin_new = 0.5*np.sum(v**2) # New kinetic energy

                delta = total.pot - total_old.pot # Cut-and-shifted potential
                delta = delta + kin_new - kin_old # Include kinetic energy change
                delta = delta / temperature       # Divide by temperature

                if metropolis ( delta ): # Accept Metropolis test
                    m_ratio  = 1.0 # Set move counter
                else:
                    r       = r_old.copy() # Restore positions: this move is rejected
                    total   = total_old    # Restore old values
                    f       = f_old        # Restore old forces
                    m_ratio = 0.0          # Set move counter

        # Calculate and accumulate variables for this step
        blk_add ( calc_variables() )

    blk_end(blk)                                          # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'    # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r*box ) # Save configuration

run_end ( calc_variables() ) # Output run averages

# Double check book-keeping for final overlap
total, f = force ( box, r_cut, r )
assert not total.ovr, 'Overlap in final configuration'

write_cnf_atoms ( cnf_prefix+out_tag, n, box, r*box ) # Save configuration
conclusion()

