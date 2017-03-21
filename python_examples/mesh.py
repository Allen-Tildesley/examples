#!/usr/bin/env python3
# mesh.py

#------------------------------------------------------------------------------------------------#
# This software was written in 2016/17                                                           #
# by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        #
# and Dominic J. Tildesley <dominic.tildesley@epfl.ch> ("the authors"),                          #
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
# There is no user input, some example parameters are set in the program.

import numpy as np

# Example values for illustration
n  = 4        # number of charges
sc = 10       # dimension of mesh
h  = 1.0 / sc # mesh spacing

# For illustration we choose random charge positions with coordinates in range (0,1)
# In a real application, we would convert positions into this range
np.random.seed()
r = np.random.random_sample( (n,3) )

# For illustration we choose +1 and -1 charges, alternately
q = np.empty(n,dtype=np.float_)
q[::2] = 1.0
q[1::2] = -1.0

rho = np.zeros( (sc,sc,sc), dtype=np.float_ ) # Zero mesh
v = np.empty( (3,3), dtype=np.float_ ) # array of weights

for i in range(n): # Loop over charges
    nr = np.rint(r[i,:]*sc).astype(np.int_) # nearest mesh point indices
    nr = np.mod(nr,sc)                      # with periodic boundaries
    dr = r[i,:] - nr * h    # vector to charge from mesh cell centre
    dr = dr - np.rint( dr ) # periodic boundaries
    dr = dr / h             # normalise by mesh cell size

    # weights for three point assignment scheme
    # first index of v is 0,1,2 for neighbour -1,0,1
    # second index of v is 0,1,2 for x,y,z
    v[0,:] = 0.5 * ( 0.5 - dr ) ** 2
    v[1,:] = 0.75 - dr**2
    v[2,:] = 0.5 * ( 0.5 + dr ) ** 2

    for i0 in range(3):            # Loop over x neighbours
        q0 = q[i]*v[i0,0]          # charge times x-weight
        n0 = np.mod(nr[0]+i0-1,sc) # x-neighbour index with periodic boundaries
        for i1 in range(3):            # Loop over y neighbours
            q1 = q0*v[i1,1]            # charge times xy-weight
            n1 = np.mod(nr[1]+i1-1,sc) # y-neighbour index with periodic boundaries
            for i2 in range(3):            # Loop over z neighbours
                q2 = q1*v[i2,2]            # charge times xyz-weight
                n2 = np.mod(nr[2]+i2-1,sc) # z-neighbour index with periodic boundaries
                rho[n0,n1,n2] = rho[n0,n1,n2] + q2 # mesh cell share of charge

rho = rho / h**3 # convert charges to charge densities

# Output charge density
format_string="{:15.6f}"*sc
for n0 in range(sc):
    print ( "{}{:5d}".format('x-layer', n0))
    for n1 in range(sc):
        print ( format_string.format(*rho[n0,n1,:]) )

# Finally check integrated charge density
print( "{}{:10.6f}{:10.6f}".format('Total charge',np.sum(q),np.sum(rho)*h**3))
