#!/usr/bin/env python3
# ewald.py

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

"""Test routines provided in ewald_module.py."""

import json
import sys
import numpy as np
from itertools        import product
from config_io_module import read_cnf_atoms
from ewald_module     import pot_r_ewald, pot_k_ewald, pot_k_pm_ewald

# Reads an atomic configuration with periodic boundary conditions from cnf.inp
# Calculates r-space and k-space contributions to potential energy
# for given screening parameter kappa and number of wave vectors determined by nk
# Adds the surface term.
# Compares with brute force summation in real space over shells of periodic boxes

print('ewald')
# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"kappa":6.0, "nk":8, "nbox":8}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ', list(defaults.keys()))
    
# Set parameters to input values or defaults
kappa = nml["kappa"] if "kappa" in nml else defaults["kappa"]
nk    = nml["nk"]    if "nk"    in nml else defaults["nk"]
nbox  = nml["nbox"]  if "nbox"  in nml else defaults["nbox"]

# Write out parameters
print ( "{:40}{:15.6f}".format('kappa*box',                           kappa) )
print ( "{:40}{:15d}"  .format('Wave-vector limit in each direction', nk)    )
print ( "{:40}{:15d}"  .format('Brute force box limit',               nbox)  )
nbox_sq = nbox**2

# Read in configuration
n, box, r = read_cnf_atoms('cnf.inp')
print("{:40}{:15d}  ".format('Number of particles', n))
print("{:40}{:15.6f}".format('Box (in sigma units)',box))
r = r / box        # Work throughout in unit box
r = r - np.rint(r) # Apply periodic boundaries
assert n%2 == 0, 'Error, we require even N for charge balance'
q = np.empty(n,dtype=np.float_)
q[::2] = 1.0
q[1::2] = -1.0
print ( "{:40}{:15.6f}".format('Net charge', np.sum(q)) )

# Compute potential energy by Ewald method
pot_r  = pot_r_ewald ( r, q, kappa )                # Real-space term involving screened Coulomb potential
pot_k  = pot_k_ewald ( nk, r, q, kappa )            # Reciprocal space term
dipole = np.sum ( q[:,np.newaxis]*r[:,:], axis=0 )  # Calculate overall box dipole
pot_s  = ( 2.0*np.pi / 3.0 ) * np.sum ( dipole**2 ) # Surface term
pot    = pot_r + pot_k + pot_s                      # Total potential
print ( "{:40}{:15.6f}".format('r-space potential energy', pot_r) )
print ( "{:40}{:15.6f}".format('k-space potential energy', pot_k) )
print ( "{:40}{:15.6f}".format('surface potential energy', pot_s) )
print ( "{:40}{:15.6f}".format('total   potential energy', pot)   )

# Compare with an illustrative (simplified) particle-mesh method
pot_k  = pot_k_pm_ewald ( nk, r, q, kappa ) # Reciprocal space term
pot    = pot_r + pot_k + pot_s              # Total potential
print ( "{:40}{:15.6f}".format('k-space potential energy (PME)', pot_k) )
print ( "{:40}{:15.6f}".format('total   potential energy (PME)', pot)   )

# Compare with brute force calculation
# Big multiple loop over all pairs and surrounding periodic boxes
# For clarity we count all pairs twice, as ij and ji
# Potentials are stored according to squared distance of neighbouring box
print('Brute force method')

pot_shell = np.zeros(nbox_sq+1, dtype=np.float_) # Zero array of potentials (not all the elements of this array will be filled)

for xbox, ybox, zbox in product ( range(-nbox,nbox+1), repeat=3 ): # Triple loop over surrounding periodic boxes
    rbox_sq = xbox**2 + ybox**2 + zbox**2
    if rbox_sq > nbox_sq: # Skip if outside maximum sphere of boxes
        continue
    rbox_vec = np.array([xbox,ybox,zbox], dtype=np.float_)
    start_shift = 0 if rbox_sq>0 else 1 # Skip all i==j if in central box
    for shift in range(start_shift,n):
        # Bare Coulomb term, including box vector, no periodic box correction
        rij     = r - np.roll(r,shift,axis=0) - rbox_vec
        rij_mag = np.sqrt(np.sum(rij**2,axis=1))
        pot_shell[rbox_sq] = pot_shell[rbox_sq] + np.sum(q*np.roll(q,shift) / rij_mag)

# Correct for double counting
pot_shell = pot_shell / 2.0

# Convert to cumulative sum
pot_shell = np.cumsum ( pot_shell )

# Write out results for increasing spherical cutoff
print('Shell      Potential')
for rbox in range(nbox+1):
    print ( "{:5d}{:15.6f}".format(rbox, pot_shell[rbox**2]) )

