#!/usr/bin/env python3
"""Dipole-quadrupole potential and forces.

This module makes available the variable n and function force.
"""

import numpy as np
n = 2 # Two-molecule potential
print('test_pot_dq module')
print('Returns potential and force for dipole-quadrupole')
print(n,'-molecule potential',sep='')

def force ( r, e ):
    """Returns potential pot and numpy arrays f, t of shape (n,3), same as input arguments.

    Demonstrates the calculation of forces from the dipole-quadrupole potential.
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

    # The dipole-quadrupole potential with mu_i = 1, Q_j = 1
    vij_dq = 1.5 * (ci*(1.0-5.0*cj**2)+2*cj*cij)/rij_mag**4

    # The quadrupole-dipole potential with Q_i = 1, mu_j = 1
    vij_qd = -1.5 * (cj*(1.0-5.0*ci**2)+2*ci*cij)/rij_mag**4

    # Forces and torque gradients for dipole-quadrupole potential with mu_i = 1, Q_j = 1
    dvdrij = -4.0*vij_dq/rij_mag
    dvdci  =  1.5 * (1-5.0*cj**2) / rij_mag**4
    dvdcj  =  3.0 * (cij-5.0*ci*cj) / rij_mag**4
    dvdcij =  3.0 * cj / rij_mag**4
    fij_dq = - dvdrij*sij - dvdci*(ei-ci*sij)/rij_mag - dvdcj*(ej-cj*sij)/rij_mag
    gi_dq  = dvdci*sij + dvdcij*ej
    gj_dq  = dvdcj*sij + dvdcij*ei

    # Forces and torque gradients for quadrupole-dipole potential with Q_i = 1, mu_j = 1
    dvdrij = -4.0*vij_qd/rij_mag
    dvdci  = -3.0 * (cij-5.0*ci*cj) / rij_mag**4
    dvdcj  = -1.5 * (1-5.0*ci**2) / rij_mag**4
    dvdcij = -3.0 * ci / rij_mag**4
    fij_qd = - dvdrij*sij - dvdci*(ei-ci*sij)/rij_mag - dvdcj*(ej-cj*sij)/rij_mag
    gi_qd  = dvdci*sij + dvdcij*ej
    gj_qd  = dvdcj*sij + dvdcij*ei

    # Final potential, forces and torques
    pot    = vij_dq + vij_qd
    f      = np.zeros_like(r)
    t      = np.zeros_like(r)
    f[i,:] =  fij_dq + fij_qd
    f[j,:] = -fij_dq - fij_qd
    t[i,:] = -np.cross(ei,gi_dq+gi_qd)
    t[j,:] = -np.cross(ej,gj_dq+gj_qd)

    return pot, f, t
