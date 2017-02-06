#!/usr/bin/env python3
"""Angle-torsion potential and forces.

This module makes available the variable n and function force.
"""

import numpy as np
n = 4 # Four-atom potential
print('test_pot_twist module')
print('Returns potential and force for polymer angle-torsion')
print(n,'-atom potential',sep='')

def force ( r ):
    """Returns potential pot and numpy array f of shape (n,3), same as input argument.

    Demonstrates the calculation of forces from the polymer angle-twisting potential.
    We choose to make the polymer the minimum length needed for testing.
    Written for ease of comparison with the text rather than efficiency!
    """

    assert r.shape == (n,3), 'Incorrect shape of r'

    d = np.zeros_like(r)             # Create d vectors (bonds)
    d[1:n,:] = r[1:n,:] - r[0:n-1,:] # Compute d vectors (zero index not used)

    # Store C coefficients in a matrix
    # In the general case we would not need to calculate every pair
    # and also we would make use of the symmetry cc[a,b]=cc[b,a]
    cc = np.zeros((n,n),dtype=np.float_)  # Create C array (scalar products)
    for a in range(1,n):
        for b in range(1,n):
            cc[a,b]=np.dot(d[a,:],d[b,:]) # Compute C array (zero indices not used)

    # Store D coefficients in a matrix
    # In the general case we would not need to calculate every pair
    # and also we would make use of the symmetry dd[a,b]=dd[b,a]
    dd = np.zeros((n,n),dtype=np.float_)  # Create D array
    for a in range(1,n):
        for b in range(1,n):
            dd[a,b]=cc[a,a]*cc[b,b] - cc[a,b]**2 # Compute D array (zero indices not used)

    a = n-1 # For this test there is just one angle

    # Here is the potential as a function of cos(phi)
    # For testing we use the simplest form: v= -cos(phi)
    # The notation matches that used in the appendix

    prefac = 1.0 / np.sqrt(dd[a,a-1]*dd[a-1,a-2])
    fac = (cc[a,a-1]*cc[a-1,a-2]-cc[a,a-2]*cc[a-1,a-1])
    pot = prefac * fac # This is -cos(phi)

    # Here we include the derivative of the potential with respect to cos(phi) in the prefactor
    # For this simple case it is -1, so the forces are simply gradients of cos(phi) as in the text
    f = np.zeros_like(r) # create force array
    fac1 = fac/dd[a,a-1]
    fac2 = fac/dd[a-1,a-2]
    f[a,:] = -prefac * ( cc[a-1,a-2]*d[a-1,:] - cc[a-1,a-1]*d[a-2,:]
             -fac1 * ( cc[a-1,a-1]*d[a,:] - cc[a,a-1]*d[a-1,:] ) )
    f[a-1,:] = -prefac * ( cc[a-1,a-2]*d[a,:] - cc[a-1,a-2]*d[a-1,:]
               + cc[a,a-1]*d[a-2,:] + cc[a-1,a-1]*d[a-2,:] - 2.0*cc[a,a-2]*d[a-1,:]
               -fac2 * ( cc[a-2,a-2]*d[a-1,:] - cc[a-1,a-2]*d[a-2,:] )
               -fac1 * ( cc[a,a]*d[a-1,:] - cc[a-1,a-1]*d[a,:] - cc[a,a-1]*d[a,:] + cc[a,a-1]*d[a-1,:]) )
    f[a-2,:] = -prefac * ( -cc[a-1,a-2]*d[a,:] + cc[a,a-1]*d[a-1,:]
               -cc[a,a-1]*d[a-2,:] - cc[a-1,a-1]*d[a,:] + 2.0*cc[a,a-2]*d[a-1,:]
               -fac2 * ( cc[a-1,a-1]*d[a-2,:] - cc[a-2,a-2]*d[a-1,:] - cc[a-1,a-2]*d[a-1,:] + cc[a-1,a-2]*d[a-2,:] )
               -fac1 * ( -cc[a,a]*d[a-1,:] + cc[a,a-1]*d[a,:] ) )
    f[a-3,:] = -prefac * ( -cc[a,a-1]*d[a-1,:] + cc[a-1,a-1]*d[a,:]
               -fac2 * ( -cc[a-1,a-1]*d[a-2,:] + cc[a-1,a-2]*d[a-1,:] ) )

    return pot, f
