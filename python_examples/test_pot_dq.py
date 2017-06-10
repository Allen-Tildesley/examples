#!/usr/bin/env python3
# test_pot_dq.py

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

"""Dipole-quadrupole potential and forces."""

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
    cij = np.dot( ei, ej  )

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
    f      = np.empty_like(r)
    t      = np.empty_like(r)
    f[i,:] =  fij_dq + fij_qd
    f[j,:] = -fij_dq - fij_qd
    t[i,:] = -np.cross(ei,gi_dq+gi_qd)
    t[j,:] = -np.cross(ej,gj_dq+gj_qd)

    return pot, f, t
