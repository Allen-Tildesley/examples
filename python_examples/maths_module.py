#!/usr/bin/env python3
"""Routines for maths, random numbers, order parameters.
"""
def random_vector():
    import numpy as np

    zeta = np.random.rand(2) # Two uniformly sampled random numbers in range (0,1)
    c = 2.0*zeta[0] - 1.0    # Random cos(theta) uniformly sampled in range (-1,+1)
    if c >= 1.0:             # Guard against very small chance of roundoff error
        s = 0.0              # Set sin(theta) to zero
    else:
        s = np.sqrt(1.0-c**2) # Calculate sin(theta) from cos(theta), always positive

    phi = zeta[1] * 2.0*np.pi # Random angle uniformly sampled in range (0,2*pi)

    return np.array ( ( s*np.cos(phi), s*np.sin(phi), c ), dtype='f8' ) # Random unit vector
