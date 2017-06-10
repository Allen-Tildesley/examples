#!/usr/bin/env python3
# diffusion.py

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

"""Calculates vacf and msd from supplied sequence of configurations."""

def unfold ( r_old, r ):
    """Removes effects of periodic boundaries on particle trajectories.

    r_old is the configuration at the previous step 
    r is the current configuration
    box is accessed from the calling program.
    The function returns the unfolded version of r.
    """

    r_new = r - r_old                      #  Convert r to displacements relative to r_old
    r_new = r_new - np.rint(r_new/box)*box # Apply periodic boundaries to displacements
    r_new = r_new + r_old                  # Convert r back to absolute coordinates
    return r_new

# diffusion program

import json
import sys
import os
import numpy as np
from config_io_module import read_cnf_atoms

print('diffusion')

# Reads a trajectory from a sequence of configuration files
# Calculates velocity autocorrelation function, mean square displacement,
# and cross-correlation between initial velocity and displacement
# Results are written to a file 'diffusion.out' with diagnostics to standard output

# For illustration and simplicity, we adopt a scheme of formatted files of the same kind
# as those that are saved at the end of each block of our MD simulation examples
# We assume that the initial configuration of a run has been copied to cnf.000
# and subsequent configurations are called cnf.001 cnf.002 etc., up to (at most) cnf.999
# Obviously, in a practical application, a binary trajectory file would fulfil this role.

# Cubic periodic boundary conditions are assumed
# r and box are assumed to be in the same units (e.g. LJ sigma)
# box is assumed to be unchanged throughout

# Note that we never apply periodic boundary conditions to the atomic positions
# We unfold the trajectory by applying PBCs to the displacements between successive configurations
# This assumes that the atoms never move more than box/2 during that interval

# Values of basic parameters are read from standard input using JSON format
# Although a default value of delta=0.05 is supplied, it is really only a place-holder
# for the correct user-supplied value (time interval between configurations)

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nt":500, "origin_interval":10, "delta":0.05}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))
    
# Set parameters to input values or defaults
nt              = nml["nt"]              if "nt"              in nml else defaults["nt"]
origin_interval = nml["origin_interval"] if "origin_interval" in nml else defaults["origin_interval"]
delta           = nml["delta"]           if "delta"           in nml else defaults["delta"]

n0 = nt // origin_interval + 1 # Enough origins to span max correlation time

# Write out parameters
print( "{:40}{:15d}  ".format('Max correlation time nt',       nt)              )
print( "{:40}{:15d}  ".format('Origin interval',               origin_interval) )
print( "{:40}{:15d}  ".format('Number of time origins n0',     n0)              )
print( "{:40}{:15.6f}".format('Time interval between configs', delta)           )

# Check that the initial configuration exists
cnf_prefix = 'cnf.'
if not os.path.isfile(cnf_prefix+'000'):
    print(cnf_prefix+'000 does not exist')
    sys.exit()

n, box, r = read_cnf_atoms(cnf_prefix+'000') # Just to get N
print("{:40}{:15d}  ".format('Number of particles', n)   )
print("{:40}{:15.6f}".format('Box (in sigma units)',box) )

msd  = np.zeros(nt+1,dtype=np.float_)
rvcf = np.zeros(nt+1,dtype=np.float_)
vacf = np.zeros(nt+1,dtype=np.float_)
norm = np.zeros(nt+1,dtype=np.float_)
t0 = np.empty(n0,dtype=np.int_)
v0 = np.empty((n0,n,3),dtype=np.float_)
r0 = np.empty((n0,n,3),dtype=np.float_)
mk   = -1 # Storage location of time origin
full = False

t = 0
while True: # Loop until configurations or naming scheme exhausted
    if t > 999:
        break
    sav_tag   = str(t).zfill(3)
    file_name = cnf_prefix+sav_tag
    if not os.path.isfile(file_name):
        break
    n, box, r, v = read_cnf_atoms(file_name,with_v=True)
    print('Processing '+file_name)

    if t>0:
        r = unfold ( r_old, r )

    if t%origin_interval == 0: # Test to store as time origin
        mk = mk + 1
        if mk >= n0:
            full = True
            mk = mk - n0 # Overwrite older values
        t0[mk]     = t # Store time origin
        r0[mk,:,:] = r # Store position at time origin
        v0[mk,:,:] = v # Store velocity at time origin

    # Correlate with all time origins stored so far
    nk = n0 if full else mk+1
    for k in range(nk): # Loop over time origins
        dt = t - t0[k]
        assert dt>=0, "dt error {:5d}".format(dt)
        if dt<= nt : # Check that dt is in range
            msd[dt]  = msd[dt]  + np.sum( (r-r0[k,:,:])**2 )        # Increment msd
            rvcf[dt] = rvcf[dt] + np.sum( (r-r0[k,:,:])*v0[k,:,:] ) # Increment cross correlation function
            vacf[dt] = vacf[dt] + np.sum( v*v0[k,:,:] )             # Increment autocorrelation function
            norm[dt] = norm[dt] + 1.0                               # Increment normalizing factor
    r_old = r     # Ready to unfold next step
    t     = t + 1 # Number of next step

assert np.all(norm>0.5), 'Normalization array error' # Should never happen

# Normalize by N as well as time-origin normalizing factors
msd  = msd  / norm / n # 3D mean-squared displacement
rvcf = rvcf / norm / n # 3D cross-correlation function
vacf = vacf / norm / n # 3D autocorrelation function

print('Output to diffusion.out')
with open("diffusion.out","w") as f:
    for t in range(nt+1):
        print("{:15.6f}{:15.8f}{:15.8f}{:15.8f}".format(t*delta,vacf[t],rvcf[t],msd[t]), file=f)
