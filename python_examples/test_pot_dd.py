#!/usr/bin/env python3
"""Dipole-dipole potential and forces.

This module makes available the variable n and function force.
"""

import numpy as np
n = 2 # Two-molecule potential
print('test_pot_dd module')
print('Returns potential and force for dipole-dipole')
print(n,'-molecule potential',sep='')

def force ( r, e ):
    """Returns potential pot and numpy arrays f, t of shape (n,3), same as input arguments.

    Demonstrates the calculation of forces from the dipole-dipole potential.
    Written for ease of comparison with the text rather than efficiency!
    """
    from math import isclose
    
    assert r.shape == (n,3), 'Incorrect shape of r'
    assert e.shape == (n,3), 'Incorrect shape of e'

    # Notation to match appendix
    i = 0
    j = 1

    ei = e[i,:]
    ej = e[j,:]
    assert isclose(np.sum(ei**2),1.0), 'Non-unit vector {} {} {}'.format(*ei)
    assert isclose(np.sum(ej**2),1.0), 'Non-unit vector {} {} {}'.format(*ej)

    rij = r[i,:] - r[j,:]
    rij_mag = np.sqrt( np.sum(rij**2) ) # Magnitude of separation vector
    sij = rij / rij_mag                 # Unit vector
    ci  = np.dot( ei, sij )
    cj  = np.dot( ej, sij )
    cij = np.dot( ei, ej  )

    # The dipole-dipole potential with mu_i = mu_j = 1
    pot = (cij-3.0*ci*cj)/rij_mag**3

    # Force vector
    fij = 3.0 * ( (cij-5.0*ci*cj)*sij + cj*ei + ci*ej ) / rij_mag**4

    # Torque gradients
    gi = ( ej - 3.0*cj*sij ) / rij_mag**3
    gj = ( ei - 3.0*ci*sij ) / rij_mag**3

    # Final forces and torques
    f = np.empty_like(r)
    t = np.empty_like(r)
    f[i,:] = fij
    f[j,:] = -fij
    t[i,:] = -np.cross(ei,gi)
    t[j,:] = -np.cross(ej,gj)

    return pot, f, t
