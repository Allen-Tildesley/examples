#!/usr/bin/env python3
# cluster.py

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

"""Identify atom clusters in a configuration."""

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
defaults = {"r_cl":1.5}
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
n, box, r = read_cnf_atoms('cluster.inp')
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
