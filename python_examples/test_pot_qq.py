#!/usr/bin/env python3
"""Quadrupole-quadrupole potential and forces.

This module makes available the variable n and function force.
"""

import numpy as np
n = 2 # Two-molecule potential
print('test_pot_qq module')
print('Returns potential and force for quadrupole-quadrupole')
print(n,'-molecule potential',sep='')

def force ( r, e ):
    """Returns potential pot and numpy arrays f, t of shape (n,3), same as input arguments.

    Demonstrates the calculation of forces from the quadrupole-quadrupole potential.
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
    cij = np.dot( ei, ej )

    # The quadrupole-quadrupole potential with Q_i = 1, Q_j = 1
    vij = 0.75 * (1.0 - 5.0*ci**2 - 5.0*cj**2 + 2.0*cij**2
          + 35.0*(ci*cj)**2 - 20.0*ci*cj*cij) / rij_mag**5

    # Forces and torque gradients for quadrupole-quadrupole potential with Q_i = 1, Q_j = 1
    dvdrij = -5.0 * vij / rij_mag
    dvdci  =  7.5 * (ci*(7.0*cj**2-1.0)-2.0*cj*cij) / rij_mag**5
    dvdcj  =  7.5 * (cj*(7.0*ci**2-1.0)-2.0*ci*cij) / rij_mag**5
    dvdcij = -3.0 * (5.0*ci*cj-cij) / rij_mag**5
    fij    = - dvdrij*sij - dvdci*(ei-ci*sij)/rij_mag - dvdcj*(ej-cj*sij)/rij_mag
    gi     = dvdci*sij + dvdcij*ej
    gj     = dvdcj*sij + dvdcij*ei

    # Final potential, forces and torques
    pot    = vij
    f      = np.zeros_like(r)
    t      = np.zeros_like(r)
    f[i,:] = fij
    f[j,:] = -fij
    t[i,:] = -np.cross(ei,gi)
    t[j,:] = -np.cross(ej,gj)

    return pot, f, t
