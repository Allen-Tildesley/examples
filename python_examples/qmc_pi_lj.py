#!/usr/bin/env python3
# qmc_pi_lj.py

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

"""Quantum Monte Carlo, path-integral method."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    # In this example we simulate using the cut (but not shifted) potential
    # but we only report results which have had the long-range corrections applied
    # The value of the cut-and-shifted potential is not used, in this example

    import numpy as np
    import math
    from averages_module import VariableType
    from lrc_module      import potential_lrc, pressure_lrc
    
    # Preliminary calculations (n,r,total are taken from the calling program)
    vol   = box**3                    # Volume
    rho   = n / vol                   # Density
    kin   = 1.5 * n * p * temperature # Average kinetic energy for NP-atom system
    kin_q = kin - total_spr           # Quantum estimator for kinetic energy
    rad_g = rad_gyr ( r )

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Acceptance ratio of atomic moves
    r_r = VariableType ( nam = 'Atomic move ratio', val = r_ratio, instant = False )

    # Acceptance ratio of centre-of-mass moves
    c_r = VariableType ( nam = 'COM move ratio', val = c_ratio, instant = False )

    # Internal energy per atom for full potential with LRC
    # LRC plus cut (but not shifted) PE already divided by factor P
    # plus KE estimator: total classical KE for NP-atom system MINUS total spring potential
    # all divided by N
    e_f = VariableType ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + (kin_q+total.pot)/n )

    # Kinetic energy per atom, just for interest's sake
    k_q = VariableType ( nam = 'KE/N', val = kin_q/n )

    # Pressure for full potential with LRC
    # LRC plus ideal gas contribution plus total virial divided by V
    kin_q = kin_q / 1.5 # Convert KE estimator to kinetic energy part of pressure
    p_f = VariableType ( nam = 'P full', val = pressure_lrc(rho,r_cut) + (kin_q+total.vir)/vol )

    # Quantum spring energy per atom, just for interest's sake
    e_q = VariableType ( nam = 'Espring/N', val = total_spr/n )

    # Quantum polymer radius of gyration, just for interest's sake
    r_g = VariableType ( nam = 'Radius of gyration', val = rad_g )

    # Collect together into a list for averaging
    return [ r_r, c_r, e_f, p_f, e_q, k_q, r_g ]

def rad_gyr ( r ):
    """Calculate average radius of gyration of polymers."""

    import numpy as np
    
    # The formula we use involves a single sweep over atoms, and is origin-independent
    # To implement periodic boundaries, we take the origin on atom 1

    r_g = 0.0 # Zero function accumulator
    p, n, d = r.shape
    assert d==3, 'Dimension error for r'

    r_g = 0.0
    
    for i in range(n): # Loop over polymers
        ri = r[1:,i,:] - r[0,i,:]     # Get coordinates of beads 1.. in polymer relative to bead 0
        ri = ri - np.rint(ri)         # Position with PBC applied (box = 1 units)
        r_cm = np.sum(ri,axis=0) / p  # centre-of-mass vector
        r_sq = np.sum(ri**2) / p      # mean-squared distance
        r_sq = r_sq - np.sum(r_cm**2) # squared radius of gyration
        if r_sq<0.0:
            r_sq=0.0 # Guard against roundoff
        r_g = r_g + np.sqrt(r_sq)    # Accumulate root-mean-square radius of gyration

    return box * r_g / n # Average RMS Rg in sigma=1 units

# Takes in a configuration of atoms (positions)
# Cubic periodic boundary conditions
# Conducts path-integral Monte Carlo at the given temperature
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in
# However, input configuration, output configuration, most calculations, and all results
# are given in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1, epsilon = 1
# The importance of quantum effects is specified through the reduced de Boer length
# deboer = (hbar/sqrt(mass*epsilon))/sigma which takes typical values of
# 0.01 for Xe, 0.03 for Ar, and 0.095 for Ne.
# This means that the quantum spring constant may be expressed k_spring = P*(T/deboer)**2
# where T stands for the reduced temperature kB*T/epsilon

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model is defined in qmc_module

import json
import sys
import numpy as np
from config_io_module import read_cnf_atoms, write_cnf_atoms
from averages_module  import run_begin, run_end, blk_begin, blk_end, blk_add
from maths_module     import random_translate_vector, metropolis
from qmc_pi_lj_module import introduction, conclusion, potential, spring, potential_1, spring_1, PotentialType

cnf_prefix = 'cnf'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('qmc_pi_lj')
print('Path-integral Monte Carlo, constant-NVT ensemble')
print('Simulation uses cut (but not shifted) potential')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "temperature":0.701087, "r_cut":2.5,
            "dr_max":0.05, "dc_max":0.1, "p":4, "deboer":0.092}
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
dc_max      = nml["dc_max"]      if "dc_max"      in nml else defaults["dc_max"]
p           = nml["p"]           if "p"           in nml else defaults["p"]
deboer      = nml["deboer"]      if "deboer"      in nml else defaults["deboer"]

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Size of each ring polymer', p)           )
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)       )
print( "{:40}{:15.6f}".format('Temperature',               temperature) )
print( "{:40}{:15.6f}".format('Potential cutoff distance', r_cut)       )
print( "{:40}{:15.6f}".format('Maximum displacement',      dr_max)      )
print( "{:40}{:15.6f}".format('Maximum COM displacement',  dc_max)      )
print( "{:40}{:15.6f}".format('de Boer length',            deboer)      )
assert p>1 and p<100, 'p must lie between 2 and 99'
k_spring = p * ( temperature / deboer ) ** 2
print( "{:40}{:15.6f}".format('Quantum spring constant', k_spring) )

# Read in initial configuration
# Read each polymer index from a unique file; we assume that n and box are all the same!
n, box, r = read_cnf_atoms ( cnf_prefix+'00.'+inp_tag) #  Read first LJ configuration to get basic parameters n and box
print( "{:40}{:15d}  ".format('Number of particles',          n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
print( "{:40}{:15.6f}".format('Density', n/box**3)  )
r = np.empty( (p,n,3), dtype=np.float_ )
for k in range(p):
    n, box, r[k,:,:] = read_cnf_atoms ( cnf_prefix+str(k).zfill(2)+'.'+inp_tag)
r = r / box           # Convert positions to box units
r = r - np.rint ( r ) # Periodic boundaries

# Calculate classical LJ and quantum spring potential energies & check overlap
total = potential ( box, r_cut, r )
assert not total.ovr, 'Overlap in initial configuration'
total_spr = spring ( box, k_spring, r )

# Initialize arrays for averaging and write column headings
r_ratio = 0.0
c_ratio = 0.0
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        r_moves = 0
        c_moves = 0

        for i in range(n): # Loop over atoms

            # Centre of mass move
            dc = random_translate_vector ( dc_max/box, np.zeros(3,dtype=np.float_) )
            partial_old = PotentialType ( 0.0, 0.0, False )
            partial_new = PotentialType ( 0.0, 0.0, False )
            for k in range(p): # Loop over ring polymer indices
                rkj = np.delete(r[k,:,:],i,0) # Array of all other atoms for this k
                partial_old = partial_old + potential_1 ( r[k,i,:], box, r_cut, rkj, p ) # Old atom classical potential etc
                rki = r[k,i,:] + dc
                rki = rki - np.rint(rki)
                partial_new = partial_new + potential_1 ( rki, box, r_cut, rkj, p ) # New atom classical potential etc
            assert not partial_old.ovr, 'Overlap in current configuration'
            if not partial_new.ovr: # Test for non-overlapping configuration
                delta = partial_new.pot - partial_old.pot # Change in classical cut (but not shifted) potential
                delta = delta / temperature
                if metropolis ( delta ): # Accept Metropolis test
                    total = total + partial_new - partial_old # Update total values
                    r[:,i,:] = r[:,i,:] + dc                  # Update positions
                    r[:,i,:] = r[:,i,:] - np.rint(r[:,i,:])   # Periodic boundary conditions
                    c_moves  = c_moves + 1                    # Increment move counter

            # Individual atom moves
            for k in range(p): # Loop over ring polymer indices
                rkj = np.delete(r[k,:,:],i,0) # Array of all the other atoms
                partial_old = potential_1 ( r[k,i,:], box, r_cut, rkj, p ) # Old atom classical potential etc
                assert not partial_old.ovr, 'Overlap in current configuration'

                kp = (k+1)%p
                km = (k-1)%p
                partial_old_spr = ( spring_1 ( r[k,i,:], r[km,i,:], box, k_spring )
                                  + spring_1 ( r[k,i,:], r[kp,i,:], box, k_spring ) )
                
                rki = random_translate_vector ( dr_max/box, r[k,i,:] ) # Trial move to new position (in box=1 units)
                rki = rki - np.rint ( rki )                            # Periodic boundary correction
                partial_new = potential_1 ( rki, box, r_cut, rkj, p )  # New atom potential, virial etc

                if not partial_new.ovr: # Test for non-overlapping configuration
                    partial_new_spr =  ( spring_1 ( rki, r[km,i,:], box, k_spring )
                                       + spring_1 ( rki, r[kp,i,:], box, k_spring ) )
                    delta = partial_new.pot - partial_old.pot         # Change in classical cut (but not shifted) potential
                    delta = delta + partial_new_spr - partial_old_spr # Add change in quantum potential
                    delta = delta / temperature
                    if metropolis ( delta ): # Accept Metropolis test
                        total     = total + partial_new - partial_old             # Update total values
                        total_spr = total_spr + partial_new_spr - partial_old_spr # Update quantum system potential
                        r[k,i,:]  = rki                                           # Update position
                        r_moves = r_moves + 1                                     # Increment move counter

        r_ratio = r_moves / (n*p)
        c_ratio = c_moves / n

        blk_add ( calc_variables() )

    blk_end(blk)                                          # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'    # Number configuration by block
    for k in range(p):
        write_cnf_atoms ( cnf_prefix+str(k).zfill(2)+'.'+sav_tag, n, box, r[k,:,:]*box ) # Save configuration

run_end ( calc_variables() )

total = potential ( box, r_cut, r ) # Double check book-keeping
assert not total.ovr, 'Overlap in final configuration'
total_spr = spring ( box, k_spring, r )

for k in range(p):
    write_cnf_atoms ( cnf_prefix+str(k).zfill(2)+'.'+out_tag, n, box, r[k,:,:]*box ) # Save configuration

conclusion()
