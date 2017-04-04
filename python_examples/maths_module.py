#!/usr/bin/env python3
# maths_module.py

#------------------------------------------------------------------------------------------------#
# This software was written in 2016/17                                                           #
# by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        #
# and Dominic J. Tildesley <dominic.tildesley@epfl.ch> ("the authors"),                          #
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

"""Routines for maths, random numbers, order parameters."""

# Routines associated with random number generation

def random_vector():
    """Returns a random unit vector as a numpy array of 3 elements."""

    import numpy as np

    zeta = np.random.rand(2) # Two uniformly sampled random numbers in range (0,1)
    c = 2.0*zeta[0] - 1.0    # Random cos(theta) uniformly sampled in range (-1,+1)
    if c >= 1.0:             # Guard against very small chance of roundoff error
        s = 0.0              # Set sin(theta) to zero
    else:
        s = np.sqrt(1.0-c**2) # Calculate sin(theta) from cos(theta), always positive

    phi = zeta[1] * 2.0*np.pi # Random angle uniformly sampled in range (0,2*pi)

    return np.array ( ( s*np.cos(phi), s*np.sin(phi), c ), dtype=np.float_ ) # Random unit vector

def random_perpendicular_vector ( old ):
    """Returns a uniformly sampled unit vector perpendicular to the old vector.""" 

    import numpy as np
    
    # Note that we do not require the reference vector to be of unit length
    # However we do require its length to be greater than a small tolerance!

    tol = 1.e-6

    assert old.size==3, 'Error in old vector dimension'
    norm = np.sum ( old**2 ) # Old squared length
    assert norm>tol, "{}{:15.3e}{:15.3e}".format('old normalization error', norm, tol)
    n = old / np.sqrt(norm) # Normalized old vector

    while True: # Loop until generated vector is not too small
       e    = random_vector () # Randomly oriented unit vector
       proj = np.dot ( e, n )  # Projection along old
       e    = e - proj * n     # Make e perpendicular to old
       norm = np.sum ( e**2 )  # Squared length
       if norm > tol:          # Accept, unless e is too small (which is unlikely)
           break

    e = e / np.sqrt ( norm ) # Normalize
    return e

def random_quaternion():
    """Returns a random unit quaternion as a numpy array of 4 elements."""

    import numpy as np
    
    while True: # Loop until within unit disk
        zeta = 2.0*np.random.rand(2) - 1.0 # Two uniform random numbers between -1 and 1
        norm1 = np.sum ( zeta**2 )         # Squared magnitude
        if norm1 < 1.0:                    # Test for within unit disk
            break

    while True: # Loop until within unit disk
        beta = 2.0*np.random.rand(2) - 1.0 # Two uniform random numbers between -1 and 1
        norm2 = np.sum ( beta**2 )         # Squared magnitude
        if norm2 < 1.0:                    # Test for within unit disk
            break

    f = np.sqrt ( ( 1.0 - norm1 ) / norm2 )
    return np.array ( ( zeta[0], zeta[1], beta[0]*f, beta[1]*f ), dtype=np.float_ ) # Random quaternion

def random_translate_vector ( dr_max, old ):
    """Returns a vector translated by a random amount."""

    import numpy as np

    # A randomly chosen vector is added to the old one

    zeta = np.random.rand(3)   # Three uniform random numbers in range (0,1)
    zeta = 2.0*zeta - 1.0      # Now in range (-1,+1)
    return old + zeta * dr_max # Move to new position

def random_rotate_vector ( angle_max, old ):
    """Returns a vector rotated by a small amount relative to the old one."""

    import numpy as np
    
    # A small randomly chosen vector is added to the old one, and the result renormalized
    # Provided angle_max is << 1, it is approximately the maximum rotation angle (in radians)
    # The magnitude of the rotation is not uniformly sampled, but this should not matter

    # Note that the old vector should be normalized and we test for this

    tol = 1.e-6

    norm = np.sum ( old**2 ) # Old squared length
    assert np.fabs ( norm - 1.0 ) < tol, "{}{:20.8e}{:20.8e}".format('old normalization error', norm, tol)

    # Choose new orientation by adding random small vector
    e    = old + angle_max * random_vector ()
    norm = np.sum ( e**2 )
    return e / np.sqrt(norm) # Normalize

def metropolis ( delta ):
    """Conduct Metropolis test, with safeguards."""

    import numpy as np
    
    exponent_guard = 75.0

    if delta > exponent_guard: # Too high, reject without evaluating
        return False
    elif delta < 0.0: # Downhill, accept without evaluating
        return True
    else:
        zeta = np.random.rand() # Uniform random number in range (0,1)
        return np.exp(-delta) > zeta # Metropolis test

# Low-level mathematical operations

def rotate_vector ( angle, axis, old ):

    """Returns a vector rotated from the old one by specified angle about specified axis."""

    import numpy as np
    from math import isclose
    
    # Note that the axis vector should be normalized and we test for this
    # In general, the old vector need not be normalized, and the same goes for the result
    # although quite often in our applications they will be

    assert old.size == 3, 'Incorrect size of old'
    assert axis.size == 3, 'Incorrect size of axis'
    assert isclose(np.sum(axis**2),1.0), 'Non-unit vector {} {} {}'.format(*axis)

    c    = np.cos ( angle )
    s    = np.sin ( angle )
    proj = np.dot ( axis, old ) # The two vectors need not be perpendicular

    # Standard (Goldstein) rotation formula
    e = c * old + ( 1.0 - c ) * proj * axis + s * np.cross ( axis, old )

    return e

def nematic_order ( e ):
    """Returns a nematic orientational order parameter."""

    import numpy as np
    
    # Calculate the nematic order parameter <P2(cos(theta))>
    # where theta is the angle between a molecular axis and the director
    # which is the direction that maximises the order parameter
    # This is obtained by finding the largest eigenvalue of
    # the 3x3 second-rank traceless order tensor

    # Note that this is not the same as the order parameter characterizing a crystal

    n, d = e.shape
    assert d==3, 'Error in e dimension '

    # Order tensor: outer product of each orientation vector, summed over n molecules
    q = np.sum ( e[:,:,np.newaxis]*e[:,np.newaxis,:], axis=0)
    q = 1.5 * q / n                # Normalize
    for i in range(3):
        q[i,i] = q[i,i] - 0.5 # Make traceless

    evals = np.linalg.eigvalsh(q)
    return evals[2]
