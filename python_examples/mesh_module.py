#!/usr/bin/env python3
# mesh_module.py

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

def mesh_function ( r, q, sc ):
    """Function to assign charges to a 3-d mesh."""

    # This function assigns a set of charges to a cubic mesh using the
    # triangular shape cloud distribution described by Hockney and Eastwood (1988)
    # It is assumed that the charges are in a box of unit length.
    # The charge mesh is indexed from 0 to sc-1 in each coordinate direction.

    import numpy as np

    n, d = r.shape
    assert d==3, 'r dimension error'
    assert n==q.size, 'q dimension error'

    h  = 1.0 / sc # mesh spacing

    rho = np.zeros( (sc,sc,sc), dtype=np.float_ ) # Zero mesh
    v   = np.empty( (3,3), dtype=np.float_ )      # Array of weights

    for i in range(n): # Loop over charges
        nr = np.rint(r[i,:]*sc).astype(np.int_) # Nearest mesh point indices
        nr = np.mod(nr,sc)                      # With periodic boundaries
        dr = r[i,:] - nr * h    # Vector to charge from mesh cell centre
        dr = dr - np.rint( dr ) # Periodic boundaries
        dr = dr / h             # Normalise by mesh cell size

        # Weights for three point assignment scheme
        # First index of v is 0,1,2 for neighbour -1,0,1
        # Second index of v is 0,1,2 for x,y,z
        v[0,:] = 0.5 * ( 0.5 - dr ) ** 2
        v[1,:] = 0.75 - dr**2
        v[2,:] = 0.5 * ( 0.5 + dr ) ** 2

        for i0 in range(3):            # Loop over x neighbours
            q0 = q[i]*v[i0,0]          # Charge times x-weight
            n0 = np.mod(nr[0]+i0-1,sc) # x-neighbour index with periodic boundaries
            for i1 in range(3):            # Loop over y neighbours
                q1 = q0*v[i1,1]            # Charge times xy-weight
                n1 = np.mod(nr[1]+i1-1,sc) # y-neighbour index with periodic boundaries
                for i2 in range(3):            # Loop over z neighbours
                    q2 = q1*v[i2,2]            # Charge times xyz-weight
                    n2 = np.mod(nr[2]+i2-1,sc) # z-neighbour index with periodic boundaries
                    rho[n0,n1,n2] = rho[n0,n1,n2] + q2 # Mesh cell share of charge

    rho = rho / h**3 # convert charges to charge densities
    return rho

def sharpen ( x ):
    """Returns sharpening factor for particle-mesh Ewald influence function."""

    import numpy as np

    tol = 1.0e-3

    if np.fabs(x) < tol:
        u = 1 - 0.5*x**2 # Taylor series
    else:
        u = ( np.sin(x)/x )**3

    return 1.0 / u**2
