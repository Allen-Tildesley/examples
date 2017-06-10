#!/usr/bin/env python3
# test_pot_at.py

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

"""Axilrod-Teller triple-dipole potential and forces."""

import numpy as np
n = 3 # Three-atom potential
print('test_pot_at module')
print('Returns potential and force for Axilrod-Teller triple-dipole')
print(n,'-atom potential',sep='')

def force ( r ):
    """Returns potential pot and numpy array f of shape (n,3), same as input argument.

    Demonstrates the calculation of forces from the Axilrod-Teller triple-dipole potential.
    Written for ease of comparison with the text rather than efficiency!
    """
    
    assert r.shape == (n,3), 'Incorrect shape of r'

    # Notation to match appendix
    i = 0
    j = 1
    k = 2

    # Note that we define the separation vectors in a cyclic way
    rij = r[i,:] - r[j,:]
    rjk = r[j,:] - r[k,:]
    rki = r[k,:] - r[i,:]
    rij_sq = np.sum(rij**2)
    rjk_sq = np.sum(rjk**2)
    rki_sq = np.sum(rki**2)
    rij2 = 1.0/rij_sq
    rjk2 = 1.0/rjk_sq
    rki2 = 1.0/rki_sq
    rij_mag = np.sqrt( rij_sq )
    rjk_mag = np.sqrt( rjk_sq )
    rki_mag = np.sqrt( rki_sq )
    ci = np.dot( rki, rij )
    cj = np.dot( rij, rjk )
    ck = np.dot( rjk, rki )
    prefac = 1.0/(rij_mag*rjk_mag*rki_mag)**5

    pot = prefac * ( rij_sq*rjk_sq*rki_sq - 3.0*ci*cj*ck ) # The triple-dipole potential with strength=nu=1

    f   = np.empty_like(r) # Create force array
    fac = 5.0*(rij_sq*rjk_sq*rki_sq-3.0*ci*cj*ck)
    f[i,:] = prefac * ( fac*(rij2*rij-rki2*rki) 
             + 3.0*ci*(ck-cj)*rjk + 3.0*cj*ck*(rki-rij) 
             + 2.0*(rij_sq*rjk_sq*rki-rjk_sq*rki_sq*rij) )
    f[j,:] = prefac * ( fac*(rjk2*rjk-rij2*rij) 
             + 3.0*cj*(ci-ck)*rki + 3.0*ck*ci*(rij-rjk) 
             + 2.0*(rjk_sq*rki_sq*rij-rki_sq*rij_sq*rjk) )
    f[k,:] = prefac * ( fac*(rki2*rki-rjk2*rjk) 
             + 3.0*ck*(cj-ci)*rij + 3.0*ci*cj*(rjk-rki) 
             + 2.0*(rki_sq*rij_sq*rjk-rij_sq*rjk_sq*rki) )

    return pot, f
