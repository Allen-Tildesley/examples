#!/usr/bin/env python3
# test_pot_bend.py

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

"""Angle-bending potential and forces."""

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
