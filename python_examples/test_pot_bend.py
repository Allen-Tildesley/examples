#!/usr/bin/env python3
"""Angle-bending potential and forces.

This module makes available the variable n and function force.
"""

import numpy as np
n = 3 # Three-atom potential
print('test_pot_bend module')
print('Returns potential and force for polymer angle-bending')
print(n,'-atom potential',sep='')

def force ( r ):
    """Returns potential pot and numpy array f of shape (n,3), same as input argument.

    Demonstrates the calculation of forces from the polymer angle-bending potential.
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

    a = n-1 # For this test there is just one angle

    # Here is the potential as a function of cos(theta)
    # For testing we use the simplest form: v= -cos(theta)
    # The notation matches that used in the appendix

    prefac = 1.0 / np.sqrt(cc[a,a]*cc[a-1,a-1])
    fac    = cc[a,a-1]
    pot    = -prefac*fac # This is -cos(theta)

    # Here we include the derivative of the potential with respect to cos(theta) in the prefactor
    # For this simple case it is -1, so the forces are simply gradients of cos(theta) as in the text
    f    = np.empty_like(r) # Create force array
    fac1 = fac / cc[a,a]
    fac2 = fac / cc[a-1,a-1]
    f[a,:]   = -prefac * ( fac1*d[a,:] - d[a-1,:] )
    f[a-1,:] =  prefac * ( fac1*d[a,:] - fac2*d[a-1,:] + d[a,:] - d[a-1,:] )
    f[a-2,:] =  prefac * ( fac2*d[a-1,:] - d[a,:] )

    return pot, f
