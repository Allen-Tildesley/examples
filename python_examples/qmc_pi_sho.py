#!/usr/bin/env python3
# qmc_pi_sho.py

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

"""Quantum Monte Carlo, path-integral, harmonic oscillator."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    from averages_module import VariableType

    # Preliminary calculations
    kin = 0.5 * p * temperature  # Kinetic energy for P-bead system

    # Move ratio
    m_r = VariableType ( nam = 'Move ratio', val = m_ratio, instant = False )

    # Classical potential energy
    pe_cl = VariableType ( nam = 'PE classical', val = pot_cl )

    # Quantum potential energy
    pe_qu = VariableType ( nam = 'PE quantum', val = pot_qu )

    # Energy
    energy = VariableType ( nam = 'Energy', val = kin+pot_cl-pot_qu )

    # Collect together into a list for averaging
    return [ m_r, pe_cl, pe_qu, energy ]

def e_pi_sho ( p, beta ):
    """Exact results for PI approximation of given order."""

    # Exact formulae given by
    # KS Schweizer, RM Stratt, D Chandler, and PG Wolynes, J Chem Phys, 75, 1347 (1981)
    # M Takahashi and M Imada, J Phys Soc Japan, 53, 3765 (1984)

    # For not-too-high P, we may express the results as a ratio of polynomials in alpha, 
    # with integer coefficients, most conveniently in partial-fraction form.
    # We give these up to P=8, and they are easy to obtain using a computer algebra package.

    # Otherwise, we use the floating-point formula, but this might become
    # inaccurate for certain values of the parameters

    import math
    import numpy as np
    
    assert p>0, 'Error in value of p'

    t = 1 / beta
    s = ( p*t ) **2

    if p==1:
        e = t
    elif p==2:
        e = 1.0 + 1.0 / np.polyval([4,1],s)
        e = e * t
    elif p==3:
        e = 1.0 + 2.0 / np.polyval([3,1],s)
        e = e * t
    elif p==4:
        e = 1.0 + 1.0 / np.polyval([4,1],s)
        e = e + 2.0 / np.polyval([2,1],s)
        e = e * t
    elif p==5:
        e = 1.0 + np.polyval([10,4],s)/np.polyval([5,5,1],s)
        e = e * t
    elif p==6:
        e = 1.0 + 1.0 / np.polyval([4,1],s)
        e = e + 2.0 / np.polyval([1,1],s)
        e = e + 2.0 / np.polyval([3,1],s)
        e = e * t
    elif p==7:
        e = 1.0 + np.polyval([28,28,6],s) / np.polyval([7,14,7,1],s)
        e = e * t
    elif p==8:
        e = 1.0 + 1.0 / np.polyval([4,1],s)
        e = e + 2.0 / np.polyval([2,1],s)
        e = e + np.polyval ([8,4],s) / np.polyval([2,4,1],s)
        e = e * t
    else:
        alpha = 0.5 * beta / p
        q1 = math.sqrt(1.0+alpha**2) + alpha
        q2 = math.sqrt(1.0+alpha**2) - alpha
        q1p = q1 ** p
        q2p = q2 ** p
        e = (q1p+q2p) / ( (q1p-q2p) * (q1+q2) )

    return e

# Program to calculate the average total energy E at temperature T
# for a particle in a harmonic potential, V=(x**2)/2,
# by simulating the discretized path integral ring polymer of P beads

# In atomic units, classical oscillation freqency omega=1, hbar=1, mass=1
# so T is equivalent to kT/hbar*omega and E is equivalent to E/hbar*omega
# Results are output as averages over the production period.
# The value of <E> may be compared with the exact result for given P for this simple problem
# as well as the exact quantum result for P=infinity.

# For this simple illustration we only use crude single-particle Metropolis moves
# It is possible to devise smarter sampling schemes for the ring polymer

# Reads several variables and options from standard input using JSON format
# Leave input line empty "{}" to accept supplied defaults

import json
import sys
import numpy as np
import math
from averages_module import run_begin, run_end, blk_begin, blk_end, blk_add
from maths_module import metropolis

print('qmc_pi_sho')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"p":8, "temperature":0.2, "nstep":50000, "nblock":20, "nequil":10, "dx_max":1.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
p           = nml["p"]           if "p"           in nml else defaults["p"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nequil      = nml["nequil"]      if "nequil"      in nml else defaults["nequil"]
dx_max      = nml["dx_max"]      if "dx_max"      in nml else defaults["dx_max"]

# Write out parameters
print( "{:40}{:15d}  ".format('Number of beads, P',                 p)           )
print( "{:40}{:15.6f}".format('Temperature',                        temperature) )
print( "{:40}{:15d}  ".format('Number of blocks for production',    nblock)      )
print( "{:40}{:15d}  ".format('Number of blocks for equilibration', nequil)      )
print( "{:40}{:15d}  ".format('Number of steps per block',          nstep)       )
print( "{:40}{:15.6f}".format('Max displacement',                   dx_max)      )
beta     = 1.0 / temperature
k_spring = p * temperature**2

np.random.seed()
x = np.zeros(p,dtype=np.float_) # Set up initial positions at origin

# Calculate initial values
pot_cl = 0.5 * np.sum ( x**2 ) / p                         # Classical potential energy
pot_qu = 0.5 * k_spring * np.sum ( ( x-np.roll(x,1) )**2 ) # Quantum potential energy

# Initialize arrays for averaging and write column headings
m_ratio = 0.0
run_begin ( calc_variables() )

for blk in range(-nequil,nblock): # Loop over blocks (including equilibration)

    blk_begin()

    for stp in range(nstep): # Loop over steps

        moves = 0
        for i in range(p): # Loop over beads

            # Identify neighbours
            ip1 = i+1 if i+1<p else 0
            im1 = i-1 if i>0   else p-1
            zeta = np.random.rand() # Uniform in range (0,1)
            zeta = 2.0*zeta - 1.0   # Now in range (-1,+1)

            xi = x[i]
            pot_cl_old = 0.5 * xi**2 / p
            pot_qu_old = 0.5 * k_spring * ( (xi-x[im1])**2 + (xi-x[ip1])**2 )
            xi         = xi + zeta * dx_max   # Trial move to new position
            pot_cl_new = 0.5 * xi**2 / p
            pot_qu_new = 0.5 * k_spring * ( (xi-x[im1])**2 + (xi-x[ip1])**2 )

            delta = ( pot_cl_new + pot_qu_new - pot_cl_old - pot_qu_old ) / temperature
            if metropolis ( delta ): # Accept Metropolis test
                pot_cl = pot_cl + pot_cl_new - pot_cl_old # Update classical potential energy
                pot_qu = pot_qu + pot_qu_new - pot_qu_old # Update quantum potential energy
                x[i]   = xi                               # Update position
                moves  = moves + 1                        # Increment move counter

        m_ratio = moves / p

        if blk >= 0:
            blk_add ( calc_variables() )

    if blk >= 0:
        blk_end(blk)

run_end ( calc_variables() )

e_qu = e_pi_sho ( p, beta )
print("{:8}{:8d}{:7}{:15.6f}".format('Exact P=',p,' energy',e_qu))
e_qu = 0.5 / math.tanh(0.5*beta)
print("{:23}{:15.6f}".format('Exact P=infinity energy',e_qu))
