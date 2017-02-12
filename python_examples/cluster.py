#!/usr/bin/env python3
"""Identify atoms clusters in a configuration."""

def in_range ( j, k ):
    """Returns True if pair, whose indices are supplied, is in range.

    Uses r, box, and squared range criterion r_cl_sq from calling program.
    """
    
    rjk=r[j,:]-r[k,:]                # Separation vector
    rjk = rjk - np.rint(rjk/box)*box # Periodic boundary conditions
    rjk_sq = np.sum(rjk**2)          # Squared separation
    return rjk_sq < r_cl_sq # Returns True if pair is in range

import json
import sys
import numpy as np
from config_io_module import read_cnf_atoms

print('cluster')
# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"r_cl":1.1}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ', list(defaults.keys()))
    
# Set parameters to input values or defaults
r_cl = nml["r_cl"] if "r_cl" in nml else defaults["r_cl"]

# Write out parameters
print ( "{:40}{:15.6f}".format('Cluster separation distance', r_cl)  )

# Read in configuration
n, box, r = read_cnf_atoms('cnf.inp')
print("{:40}{:15d}  ".format('Number of particles', n))
print("{:40}{:15.6f}".format('Box (in sigma units)',box))

r       = r - np.rint(r/box)*box     # Apply periodic boundaries
r_cl_sq = r_cl**2                    # Used in in_range function
my_list = np.arange(n,dtype=np.int_) # Set up the list

for i in range(n-1): # Begin outer loop
    if i == my_list[i]:
        j=i
        while True: # Begin inner loop
            for k in range(i+1,n): # Begin innermost loop
                if my_list[k]==k:
                    if in_range(j,k):
                        my_list[k], my_list[j] = my_list[j], my_list[k] # Swap elements
            j = my_list[j]
            if j==i:
                break

# For diagnostic purposes, print out the cluster membership
# no particular sorting (e.g. by size)

done = np.zeros(n,dtype=np.int_)
cluster_id = 0

print('Cluster members .....')

while True: # Begin loop over remaining clusters
    if np.all(done>0): # Loop until all done
        break
    i=np.where(done==0)[0][0] # Find first zero
    cluster_id=cluster_id+1
    print("{}{:5d}{}".format('Cluster',cluster_id,' = '),end="")
    j = i
    done[j] = cluster_id
    print("{:5d}".format(j),end="")
    while True: # Begin loop to find other members of cluster
        j=my_list[j]
        if j==i: # link list has returned to start
            break
        done[j] = cluster_id
        print("{:5d}".format(j),end="")
    print()

# Count cluster members
print('Cluster Count')
for i in range(cluster_id):
    print("{:7d}{:5d}".format(i+1,np.count_nonzero(done==i+1)))
