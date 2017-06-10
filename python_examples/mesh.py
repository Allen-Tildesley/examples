#!/usr/bin/env python3
# mesh.py

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

"""Assignment of charges to a 3-d mesh."""

# This program assigns a set of charges to a cubic mesh using the
# triangular shape cloud distribution described by Hockney and Eastwood (1988)
# The charges are positioned in a box of unit length.
# The charge mesh is indexed from 0 to sc-1 in each coordinate direction.

import json
import sys
import numpy as np
from mesh_module import mesh_function

print('mesh')
# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"n":4, "sc":10}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ', list(defaults.keys()))
    
# Set parameters to input values or defaults
n = nml["n"]   if "n"  in nml else defaults["n"]
sc = nml["sc"] if "sc" in nml else defaults["sc"]
h  = 1.0 / sc # mesh spacing

# Write out parameters
print ( "{:40}{:15d}"  .format('Number of charges', n ) )
print ( "{:40}{:15d}"  .format('Dimension of mesh', sc) )
print ( "{:40}{:15.6f}".format('Mesh spacing',      h ) )

# For illustration we choose random charge positions with coordinates in range (0,1)
# In a real application, we would convert positions into this range
np.random.seed()
r = np.random.random_sample( (n,3) )

# For illustration we choose +1 and -1 charges, alternately
q = np.empty(n,dtype=np.float_)
q[::2] = 1.0
q[1::2] = -1.0

rho = mesh_function ( r, q, sc )

# Output charge density
format_string="{:15.6f}"*sc
for n0 in range(sc):
    print ( "{}{:5d}".format('x-layer', n0))
    for n1 in range(sc):
        print ( format_string.format(*rho[n0,n1,:]) )

# Finally check integrated charge density
print( "{}{:10.6f}{:10.6f}".format('Total charge',np.sum(q),np.sum(rho)*h**3))
