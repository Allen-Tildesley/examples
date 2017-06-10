#!/usr/bin/env python3
# pair_distribution.py

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

"""Calculates pair distribution function g(r)."""

import json
import sys
import numpy as np
import math
from config_io_module import read_cnf_atoms
import os.path

print('pair_distribution')
# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"dr":0.02}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ', list(defaults.keys()))
    
# Set parameters to input values or defaults
dr = nml["dr"] if "dr" in nml else defaults["dr"]

# Write out parameters
print ( "{:40}{:15.6f}".format('g(r) spacing dr', dr)  )

# Read in configuration
cnf_prefix = 'cnf.'
if not os.path.isfile(cnf_prefix+'000'):
    print(cnf_prefix+'000 does not exist')
    sys.exit()
n, box, r = read_cnf_atoms(cnf_prefix+'000')
print("{:40}{:15d}  ".format('Number of particles', n)   )
print("{:40}{:15.6f}".format('Box (in sigma units)',box) )
dr = dr / box           # Convert to box=1 units
nk = math.floor(0.5/dr) # Accumulate out to half box length
r_max = nk*dr           # Actual r_max (box=1 units)
print( "{:40}{:15d}  ".format('Number of bins', nk)    )
print( "{:40}{:15.6f}".format('Maximum r/box',  r_max) )

h     = np.zeros(nk,dtype=np.int_) # Histogram bins initialized to zero
nstep = 0                          # Counts configurations

while True: # Loop until configurations or naming scheme exhausted
    if nstep >= 1000:
        break
    sav_tag   = str(nstep).zfill(3)
    file_name = cnf_prefix+sav_tag
    if not os.path.isfile(file_name):
        break
    n, box, r = read_cnf_atoms(file_name)
    print('Processing '+file_name)
    r = r / box # Convert to box=1 units

    # To make best use of NumPy, we loop over cyclic offset shifts and process N rij pairs at once.
    # factor=2 accounts for both ij and ji, but if N is even, the last shift 
    # includes both ij and ji pairs already, so factor=1
    # The idea dates back to S Brode and R Ahlrichs Comput Phys Commun 42, 51 (1986)
    nshift=n//2
    for shift in range(nshift):
        rij        = r - np.roll(r,shift+1,axis=0)                   # Difference between N pairs of coordinates
        rij        = rij - np.rint(rij)                              # Apply periodic boundaries
        rij_mag    = np.sqrt(np.sum(rij**2,axis=1))                  # Separation distances
        hist,edges = np.histogram(rij_mag,bins=nk,range=(0.0,r_max)) # Accumulate histogram of separations
        factor     = 1 if n%2==0 and shift==nshift-1 else 2          # Normally 2 except possibly on last shift
        h          = h + factor*hist                                 # Accumulate histogram
    nstep=nstep+1 # Increment configuration counter ready for next time

rho  = float(n) # Our calculation is done in box=1 units
h_id = ( 4.0 * np.pi * rho / 3.0) * ( edges[1:nk+1]**3 - edges[0:nk]**3 ) # Ideal number
g    = h / h_id / (n*nstep) # Average number

print('Output to pair_distribution.out')
edges = edges*box                       # Convert bin edges back to sigma=1 units
r_mid = 0.5*(edges[0:nk]+edges[1:nk+1]) # Mid points of bins
np.savetxt('pair_distribution.out',np.c_[r_mid,g],fmt="%15.8f")
