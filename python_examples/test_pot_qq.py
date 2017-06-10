#!/usr/bin/env python3
# test_pot_qq.py

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

"""Quadrupole-quadrupole potential and forces."""

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
    cij = np.dot( ei, ej  )

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
    f      = np.empty_like(r)
    t      = np.empty_like(r)
    f[i,:] = fij
    f[j,:] = -fij
    t[i,:] = -np.cross(ei,gi)
    t[j,:] = -np.cross(ej,gj)

    return pot, f, t
