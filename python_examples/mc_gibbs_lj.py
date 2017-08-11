#!/usr/bin/env python3
# mc_gibbs_lj.py

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

"""Monte Carlo, Gibbs ensemble."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    # In this example we simulate using the cut (but not shifted) potential
    # The values of < p_c >, < e_c > and < density > should be consistent (for this potential)
    # For simplicity, long-range corrections are not applied here to give estimates of
    # < e_f > and < p_f > for the full (uncut) potential, but this is straightforward to do.
    # The value of the cut-and-shifted potential is not used, in this example

    import numpy as np
    import math
    from averages_module import VariableType
    from lrc_module      import potential_lrc, pressure_lrc, pressure_delta

    # Preliminary calculations (n1,n2,r1,r2,etc are taken from the calling program)
    vol1 = box1**3    # Volume
    vol2 = box2**3    # Volume
    rho1 = n1 / vol1  # Density
    rho2 = n2 / vol2  # Density

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Move, swap, volume exchange acceptance ratios
    m1_r  = VariableType ( nam = 'Move ratio (1)',    val = m1_ratio,  instant = False )
    m2_r  = VariableType ( nam = 'Move ratio (2)',    val = m2_ratio,  instant = False )
    x12_r = VariableType ( nam = 'Swap ratio (1->2)', val = x12_ratio, instant = False )
    x21_r = VariableType ( nam = 'Swap ratio (2->1)', val = x21_ratio, instant = False )
    v_r   = VariableType ( nam = 'Volume ratio',      val = v_ratio,   instant = False )

    # Number of particles
    n_1 = VariableType ( nam = 'Number (1)', val = float(n1) )
    n_2 = VariableType ( nam = 'Number (2)', val = float(n2) )

    # Density
    density_1 = VariableType ( nam = 'Density (1)', val = rho1 )
    density_2 = VariableType ( nam = 'Density (2)', val = rho2 )

    # Internal energy per atom for simulated, cut, potential
    # Ideal gas contribution plus cut (but not shifted) PE divided by N
    e1_c = VariableType ( nam = 'E/N cut (1)', val = 1.5*temperature + total1.pot/n1 )
    e2_c = VariableType ( nam = 'E/N cut (2)', val = 1.5*temperature + total2.pot/n2 )

    # Pressure for simulated, cut, potential
    # Delta correction plus ideal gas contribution plus total virial divided by V
    p1_c = VariableType ( nam = 'P cut (1)', val = pressure_delta(rho1,r_cut) + rho1*temperature + total1.vir/vol1 )
    p2_c = VariableType ( nam = 'P cut (2)', val = pressure_delta(rho2,r_cut) + rho2*temperature + total2.vir/vol2 )

    # Collect together into a list for averaging
    return [ m1_r, m2_r, x12_r, x21_r, v_r, n_1, n_2, density_1, density_2, e1_c, e2_c, p1_c, p2_c ]

# Takes in a pair of configurations of atoms (positions)
# Cubic periodic boundary conditions
# Conducts Gibbs ensemble Monte Carlo at the given temperature, total volume and total N
# To avoid some inconvenient tests, we disallow configurations in which either box is empty
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in
# However, input configuration, output configuration, most calculations, and all results
# are given in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1, epsilon = 1

# Note that long-range corrections are not included in the acceptance/rejection
# of creation and destruction moves

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

cnf1_prefix = 'cnf1.'
cnf2_prefix = 'cnf2.'
inp_tag     = 'inp'
out_tag     = 'out'
sav_tag     = 'sav'

print('mc_gibbs_lj')
print('Monte Carlo, Gibbs ensemble')
print('Simulation uses cut (but not shifted) potential')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "nswap":20, "temperature":1.0,
            "r_cut":2.5, "dr_max":0.15, "dv_max":10.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
nswap       = nml["nswap"]       if "nswap"       in nml else defaults["nswap"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
r_cut       = nml["r_cut"]       if "r_cut"       in nml else defaults["r_cut"]
dr_max      = nml["dr_max"]      if "dr_max"      in nml else defaults["dr_max"]
dv_max      = nml["dv_max"]      if "dv_max"      in nml else defaults["dv_max"]

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',              nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block',     nstep)       )
print( "{:40}{:15d}  ".format('Swap attempts per step',        nswap)       )
print( "{:40}{:15.6f}".format('Specified temperature',         temperature) )
print( "{:40}{:15.6f}".format('Potential cutoff distance',     r_cut)       )
print( "{:40}{:15.6f}".format('Maximum displacement',          dr_max)      )
print( "{:40}{:15.6f}".format('Maximum volume change',         dv_max)      )

# Read in initial configurations
n1, box1, r1 = read_cnf_atoms ( cnf1_prefix+inp_tag )
n2, box2, r2 = read_cnf_atoms ( cnf2_prefix+inp_tag )
print( "{:40}{:15d}{:15d}    ".format('Number of particles',   n1,         n2         ) )
print( "{:40}{:15.6f}{:15.6f}".format('Simulation box length', box1,       box2       ) )
print( "{:40}{:15.6f}{:15.6f}".format('Density',               n1/box1**3, n2/box2**3 ) )
r1 = r1 / box1           # Convert positions to box units
r2 = r2 / box2           # Convert positions to box units
r1 = r1 - np.rint ( r1 ) # Periodic boundaries
r2 = r2 - np.rint ( r2 ) # Periodic boundaries

# Initial energy and overlap check
total1 = potential ( box1, r_cut, r1 )
assert not total1.ovr, 'Overlap in initial configuration 1'
total2 = potential ( box2, r_cut, r2 )
assert not total2.ovr, 'Overlap in initial configuration 2'

# Initialize arrays for averaging and write column headings
m1_ratio  = 0.0
m2_ratio  = 0.0
x12_ratio = 0.0
x21_ratio = 0.0
v_ratio   = 0.0
run_begin ( calc_variables() )

# Initialize histograms
nh = 300
rho_min, rho_max = 0.0, 0.9
eng_min, eng_max = -3.3, 1.2
rho_vals = np.zeros ( (nstep,2), dtype = np.float_ ) # Stores density values for both boxes
eng_vals = np.zeros ( (nstep,2), dtype = np.float_ ) # Stores energy values for both boxes
rho_hist = np.zeros ( nh, dtype = np.float_ ) # Density histogram
eng_hist = np.zeros ( nh, dtype = np.float_ ) # Energy histogram

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        m_acc = 0

        for i in range(n1): # Loop over atoms in system 1

            rj = np.delete(r1,i,0)    # Array of all the other atoms
            partial_old = potential_1 ( r1[i,:], box1, r_cut, rj ) # Old atom potential, virial etc
            assert not partial_old.ovr, 'Overlap in current configuration'

            ri = random_translate_vector ( dr_max/box1, r1[i,:] ) # Trial move to new position (in box=1 units)
            ri = ri - np.rint ( ri )                              # Periodic boundary correction
            partial_new = potential_1 ( ri, box1, r_cut, rj )     # New atom potential, virial etc

            if not partial_new.ovr: # Test for non-overlapping configuration
                delta = partial_new.pot - partial_old.pot # Use cut (but not shifted) potential
                delta = delta / temperature

                if metropolis ( delta ): # Accept Metropolis test
                    total1  = total1 + partial_new - partial_old # Update total values
                    r1[i,:] = ri                                 # Update position
                    m_acc   = m_acc + 1                          # Increment move counter

        m1_ratio = m_acc / n1

        m_acc = 0

        for i in range(n2): # Loop over atoms in system 2

            rj = np.delete(r2,i,0)    # Array of all the other atoms
            partial_old = potential_1 ( r2[i,:], box2, r_cut, rj ) # Old atom potential, virial etc
            assert not partial_old.ovr, 'Overlap in current configuration'

            ri = random_translate_vector ( dr_max/box2, r2[i,:] ) # Trial move to new position (in box=1 units)
            ri = ri - np.rint ( ri )                              # Periodic boundary correction
            partial_new = potential_1 ( ri, box2, r_cut, rj )     # New atom potential, virial etc

            if not partial_new.ovr: # Test for non-overlapping configuration
                delta = partial_new.pot - partial_old.pot # Use cut (but not shifted) potential
                delta = delta / temperature

                if metropolis ( delta ): # Accept Metropolis test
                    total2  = total2 + partial_new - partial_old # Update total values
                    r2[i,:] = ri                                 # Update position
                    m_acc   = m_acc + 1                          # Increment move counter

        m2_ratio = m_acc / n2

        x12_try = 0
        x12_acc = 0
        x21_try = 0
        x21_acc = 0

        for iswap in range(nswap):
            
            ri = np.random.rand(3) # Three uniform random numbers in range (0,1)
            ri = ri - 0.5          # Now in range (-0.5,+0.5) for box=1 units

            if np.random.rand() < 0.5: # Try swapping 1->2
                
                x12_try = x12_try + 1

                if n1>1: # Disallow n1->0
                    i  = np.random.randint(n1) # Choose atom at random in system 1
                    rj = np.delete(r1,i,0)     # Array of all the other atoms
                    partial_old = potential_1 ( r1[i,:], box1, r_cut, rj ) # Old atom potential, virial, etc
                    assert not partial_old.ovr, 'Overlap found on particle removal'
                    partial_new = potential_1 ( ri, box2, r_cut, r2 ) # New atom potential, virial, etc

                    if not partial_new.ovr: # Test for non-overlapping configuration
                        delta = ( partial_new.pot - partial_old.pot ) / temperature # Use cut (not shifted) potential
                        delta = delta - np.log ( box2**3 / ( n2+1 ) ) # Creation in 2
                        delta = delta + np.log ( box1**3 / n1 )       # Destruction in 1

                        if metropolis ( delta ): # Accept Metropolis test
                            r2      = np.append ( r2, ri[np.newaxis,:], 0 ) # Add new particle to r2 array
                            n2      = r2.shape[0]                           # New value of N2
                            r1      = np.copy(rj)                           # Delete particle from r1 array
                            n1      = r1.shape[0]                           # New value of N1
                            total1  = total1 - partial_old                  # Update total values
                            total2  = total2 + partial_new                  # Update total values
                            x12_acc = x12_acc + 1                           # Increment 1->2 move counter

            else: # Try swapping 2->1

                x21_try = x21_try + 1

                if n2>1: # Disallow n2->0
                    i  = np.random.randint(n2) # Choose atom at random in system 2
                    rj = np.delete(r2,i,0)     # Array of all the other atoms
                    partial_old = potential_1 ( r2[i,:], box2, r_cut, rj ) # Old atom potential, virial, etc
                    assert not partial_old.ovr, 'Overlap found on particle removal'
                    partial_new = potential_1 ( ri, box1, r_cut, r1 ) # New atom potential, virial, etc

                    if not partial_new.ovr: # Test for non-overlapping configuration
                        delta = ( partial_new.pot - partial_old.pot ) / temperature # Use cut (not shifted) potential
                        delta = delta - np.log ( box1**3 / ( n1+1 ) ) # Creation in 1
                        delta = delta + np.log ( box2**3 / n2 )       # Destruction in 2

                        if metropolis ( delta ): # Accept Metropolis test
                            r1      = np.append ( r1, ri[np.newaxis,:], 0 ) # Add new particle to r1 array
                            n1      = r1.shape[0]                           # New value of N1
                            r2      = np.copy(rj)                           # Delete particle from r2 array
                            n2      = r2.shape[0]                           # New value of N
                            total1  = total1 + partial_new                  # Update total values
                            total2  = total2 - partial_old                  # Update total values
                            x21_acc = x21_acc + 1                           # Increment 2->1 move counter

        x12_ratio = x12_acc/x12_try if x12_try>0 else 0.0
        x21_ratio = x21_acc/x21_try if x21_try>0 else 0.0

        # Volume move

        v_ratio = 0.0

        dv       = dv_max * ( 2.0*np.random.rand() - 1.0 ) # Uniform on (-dv_max,+dv_max)
        vol1_old = box1**3                                 # Old volume
        vol2_old = box2**3                                 # Old volume
        vol1_new = vol1_old - dv                           # New volume
        vol2_new = vol2_old + dv                           # New volume
        box1_new = vol1_new**(1.0/3.0)                     # New box length
        box2_new = vol2_new**(1.0/3.0)                     # New box length
        assert min(box1_new,box2_new)>2.0*r_cut, 'Box length too small'
        total1_new = potential ( box1_new, r_cut, r1 )
        total2_new = potential ( box2_new, r_cut, r2 )

        if not ( total1_new.ovr or total2_new.ovr ): # Test for non-overlapping configurations
            delta = total1_new.pot + total2_new.pot - total1.pot - total2.pot
            delta = delta / temperature
            delta = delta - n1*np.log(vol1_new/vol1_old) # Volume scaling in system 1
            delta = delta - n2*np.log(vol2_new/vol2_old) # Volume scaling in system 2
            if metropolis ( delta ): # Accept Metropolis test
                total1  = total1_new # Update total values
                total2  = total2_new # Update total values
                box1    = box1_new   # Update box lengths
                box2    = box2_new   # Update box lengths
                v_ratio = 1.0        # Set move counter

        blk_add ( calc_variables() )
        rho_vals[stp,:] = [n1/box1**3,n2/box2**3]
        eng_vals[stp,:] = [1.5*temperature+total1.pot/n1,1.5*temperature+total2.pot/n2]

    blk_end(blk)                                               # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'         # Number configuration by block
    write_cnf_atoms ( cnf1_prefix+sav_tag, n1, box1, r1*box1 ) # Save configuration
    write_cnf_atoms ( cnf2_prefix+sav_tag, n2, box2, r2*box2 ) # Save configuration

    # Increment histograms
    rho_h, rho_bins = np.histogram ( rho_vals, bins=nh, range=(rho_min,rho_max) )
    eng_h, eng_bins = np.histogram ( eng_vals, bins=nh, range=(eng_min,eng_max) )
    rho_hist = rho_hist + rho_h
    eng_hist = eng_hist + eng_h

run_end ( calc_variables() )

# Write out histograms
norm = 2*nstep*nblock*(rho_bins[1]-rho_bins[0])
rho_hist = rho_hist / norm
norm = 2*nstep*nblock*(eng_bins[1]-eng_bins[0])
eng_hist = eng_hist / norm
with open("his.out","w") as f:
    for k in range(nh):
        rho = (rho_bins[k]+rho_bins[k+1])/2.0
        eng = (eng_bins[k]+eng_bins[k+1])/2.0
        print("{:15.6f}{:15.6f}{:15.6f}{:15.6f}".format(rho,rho_hist[k],eng,eng_hist[k]),file=f)

write_cnf_atoms ( cnf1_prefix+out_tag, n1, box1, r1*box1 ) # Save configuration
write_cnf_atoms ( cnf2_prefix+out_tag, n2, box2, r2*box2 ) # Save configuration
conclusion()
