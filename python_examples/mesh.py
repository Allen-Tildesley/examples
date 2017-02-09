#!/usr/bin/env python3
# mesh program
"""Assignment of charges to a 3-d mesh.

This program assigns a set of charges to a cubic mesh using the
triangular shape cloud distribution described by Hockney and Eastwood (1988)
The charges are positioned in a box of unit length.
The charge mesh is indexed from 0 to sc-1 in each coordinate direction.
There is no user input, some example parameters are set in the program.
"""
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
