#!/usr/bin/env python3
# mc_chain_sw_module.py

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

"""Monte Carlo, single chain, square wells."""

fast = True # Change this to replace NumPy potential evaluation with slower Python

def introduction():
    """Prints out introductory statements at start of run."""

    print('Hard-sphere chain with fixed bond length')
    print('Square-well attractive potential')
    print('Diameter, sigma = 1')
    if fast:
        print('Fast NumPy potential routine')
    else:
        print('Slow Python potential routine')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def regrow ( s, m_max, k_max, bond, q_range, r, q ):
    """Carries out single regrowth move, returning new r, q and indicator of success."""

    # A short sequence of m atoms (m<=m_max) is deleted and regrown in the CBMC manner
    # We randomly select which end of the chain to apply each of these operations to
    # At each stage, k_max different atom positions are tried
    # The bond length is fixed throughout
    # Weights used in the regrowth are athermal, computed only on the basis of the
    # hard-core overlap part of the non-bonded interactions: essentially they count non-overlaps
    # Hence they are suitable for use in both NVT and Wang-Landau simulations

    # r_old and r_new are used as working arrays

    import numpy as np
    from maths_module import random_vector

    w_tol = 1.e-10 # Min weight tolerance
    n, d = r.shape
    assert d==3, 'Dimension error in regrow'
    r_try = np.empty((k_max,3),dtype=np.float_)
    w     = np.empty(k_max,dtype=np.float_)

    r_old = np.copy(r) # Store copy of r
    q_old = q          # Store old q

    if m_max <= 0:
        return r_old, q_old, False
    
    m = 1+np.random.randint ( m_max ) # Number of atoms to regrow
    c = np.random.randint ( 4 )       # Growth option

    # PART 1: CONSTRUCT NEW CONFIGURATION WITH NEW WEIGHT

    if c==0: # Remove from end and add to end
        r[:n-m,:] = r_old[:n-m,:] # Copy first n-m atoms

    elif c==1: # Remove from end and add to start
        r[:n-m,:] = r_old[n-m-1::-1,:] # Copy and reverse first n-m atoms

    elif c==2: # Remove from start and add to start
        r[:n-m,:] = r_old[:m-1:-1,:] # Copy and reverse last n-m atoms

    else: # Remove from start and add to end
        r[:n-m,:] = r_old[m:,:] # Copy last n-m atoms

    # Take the opportunity to place atom 0 at the origin
    r0        = np.copy(r[0,:])
    r[:n-m,:] = r[:n-m,:] - r0

    w_new = np.float_(1.0)
    for i in range(n-m,n): # Loop to regrow last m atoms, computing new weight
        for k in range(k_max): # Loop over k_max tries
            r_try[k,:] = r[i-1,:] + bond * random_vector()  # Trial position in random direction from i-1
            w[k]       = weight_1 ( r_try[k,:], r[:i-1,:] ) # Store overlap weight for this try
        w_sum = np.sum(w)
        if w_sum<w_tol: # Early exit if this happens at any stage
            return r_old, q_old, False
        w = w / w_sum
        k = np.random.choice(k_max,p=w) # Pick winning try according to weights
        r[i,:] = r_try[k,:]             # Store winning position
        w_new = w_new * w_sum           # Accumulate total weight

    if w_new<w_tol: # Exit if this happens
        return r_old, q_old, False

    q_new = qcount ( r, q_range ) # Compute full new nonbonded energy
    r_new = np.copy(r)            # Store new configuration

    # END OF PART 1: NEW CONFIGURATION AND WEIGHT ARE COMPLETE

    # PART 2: RECONSTRUCT OLD CONFIGURATION WITH OLD WEIGHT

    if c==0 or c==1: # Remove and hence reconstruct at end
        r[:,:] = r_old[:,:] # Copy all n atoms
    else: # Remove and reconstruct at start
        r[:,:] = r_old[::-1,:] # Copy and reverse all n atoms

    w_old = np.float_(1.0)
    for i in range(n-m,n):

        # Old position and weight are stored as try 0
        r_try[0,:] = r[i,:]
        w[0]       = 1 # Current weight must be 1

        # Remaining tries only required to compute weight
        for k in range(1,k_max): # Loop over k_max-1 other tries
            r_try[k,:] = r[i-1,:] + bond * random_vector()  # Trial position in random direction from i-1
            w[k]       = weight_1 ( r_try[k,:], r[:i-1,:] ) # Store overlap weight for this try

        w_sum  = np.sum(w)
        r[i,:] = r_try[0,:]    # Restore winning position (always the original one)
        w_old  = w_old * w_sum # Accumulate total weight

    assert w_old>w_tol, 'Old weight error'

    # END OF PART 2: OLD CONFIGURATION AND WEIGHT ARE COMPLETE

    # Choose either old or new configuration

    if accept ( s, q_old, q_new, w_old, w_new ):
        return r_new, q_new, True
    else:
        return r_old, q_old, False

def pivot ( s, phi_max, q_range, r, q ):
    """Carries out a pivot move, returning new r, q and indicator of success."""

    # An atom is picked at random, and the part of the chain lying to one side of it
    # is rotated as a whole, by a random angle, about a randomly oriented axis
    # There are no weights to take into account in the acceptance/rejection decision
    # (the function weight_1 is used simply to indicate overlap / no overlap)

    # r_old and r_new are used as working arrays

    import numpy as np
    from maths_module import random_vector, rotate_vector

    n, d = r.shape
    assert d==3, 'Dimension error in regrow'
    if n<3:
        return r, q, False # Makes no sense to pivot such a short chain

    r_old = np.copy(r) # Store copy of r
    q_old = q          # Store old q
        
    c = np.random.randint ( 2 ) # Which part to pivot (actually redundant here)

    if c==0: # Copy atoms
        r[:,:] = r_old[:,:]

    else: # Copy and reverse atoms
        r[:,:] = r_old[::-1,:]

    # Take the opportunity to place atom 0 at the origin
    r0 = np.copy(r[0,:])
    r  = r - r0

    j = np.random.randint ( 1, n-1 ) # Pivot position (not at either end)
    rj = r[j,:] # Pivot point

    for i in range(1,j+1): # Loop over static atoms, redundant but for roundoff
        if weight_1 ( r[i,:], r[:i-1,:] )==0: # Check for overlaps
            return r_old, q_old, False

    # Pivot, and check for overlap in new configuration
    # NB include overlap checks within rotated segment, because of roundoff

    axis = random_vector ()                         # Pivot rotation axis
    phi  = phi_max * ( 2.0*np.random.rand() - 1.0 ) # Pivot rotation angle in desired range

    for i in range(j+1,n): # Loop over moving atoms

       rij    = r[i,:] - rj                      # Relative vector of atom
       rij    = rotate_vector ( phi, axis, rij ) # Rotate relative vector
       r[i,:] = rj + rij                         # New absolute position

       if weight_1 ( r[i,:], r[:i-1,:] )==0: # Check for overlaps
           return r_old, q_old, False

    q_new = qcount ( r, q_range ) # Compute full new nonbonded energy
    r_new = np.copy(r)          # Store new configuration

    # Choose either old or new configuration

    if accept ( s, q_old, q_new ):
        return r_new, q_new, True
    else:
        return r_old, q_old, False

def crank ( s, phi_max, q_range, r, q ):
    """Carries out a crankshaft move, returning new r, q and indicator of success."""

    # An atom is picked at random. Unless it is an end-atom, a rotation axis
    # is defined as the line joining the two atoms on either side, and the
    # chosen atom is rotated about that axis by a random angle.
    # In the case of end-atoms, the axis is defined by the line joining its
    # nearest and next-nearest neighbours.
    # There are no weights to take into account in the acceptance/rejection decision
    # (the function weight_1 is used simply to indicate overlap / no overlap)

    # r_old and r_new are used as working arrays

    import numpy as np
    from maths_module import rotate_vector

    n, d = r.shape
    assert d==3, 'Dimension error in regrow'
    if n<3:
        return r, q, False # Makes no sense to crank such a short chain

    r_old = np.copy(r) # Store copy of r
    q_old = q          # Store old q

    r[:,:] = r_old[:,:] # Copy old position (somewhat redundant here)

    i = np.random.randint ( n ) # Pick random atom to move

    if i == 0:                           # Rotate about 1--2 bond
       axis = r[1,:] - r[2,:]            #   axis of rotation
       rj   = r[1,:]                     #   reference position on axis
       rnot = r[2:,:]                    #   other atoms (not 0 or 1)
    elif i == n-1:                       # Rotate about (n-3)--(n-2) bond
       axis = r[n-2,:] - r[n-3,:]        #   axis of rotation
       rj   = r[n-2,:]                   #   reference position on axis
       rnot = r[:n-2,:]                  #   other atoms (not n-1 or n-2)
    else:                                # Rotate about (i-1)--(i+1) bond
       axis = r[i+1,:] - r[i-1,:]        #   axis of rotation
       rj   = r[i-1,:]                   #   reference position on axis
       rnot = np.delete(r,[i-1,i,i+1],0) #   other atoms (not i-1, i, i+1)

    q_new = q_old - qcount_1 ( r[i,:], rnot, q_range ) # Subtract old energies for moving atom

    norm = np.sqrt(np.sum(axis**2)) # Squared length of rotation axis
    axis = axis / norm              # Unit vector along rotation axis

    phi  = phi_max * ( 2.0*np.random.rand() - 1.0 ) # Rotation angle in desired range

    rij    = r[i,:] - rj                      # Relative vector of atom
    rij    = rotate_vector ( phi, axis, rij ) # Rotate relative vector
    r[i,:] = rj + rij                         # New absolute position

    if weight_1 ( r[i,:], rnot )==0 : # Check for overlaps
        return r_old, q_old, False

    q_new = q_new + qcount_1 ( r[i,:], rnot, q_range ) # Add new energies for moving atom
    r_new = np.copy(r)                                 # Store new configuration
    
    # Choose either old or new configuration

    if accept ( s, q_old, q_new ):
        return r_new, q_new, True
    else:
        return r_old, q_old, False

def weight ( r ):
    """Takes in coordinate array, and returns weight = 0 (overlap) or 1 (no overlap)."""

    # Arithmetically w is the product of all the pair Boltzmann factors
    # Here, each is 0 or 1, which allows us to treat them as Boolean variables
    # but in the more general case they would be real and we would multiply them all together

    n, d = r.shape
    assert d==3, 'Dimension error for r in potential'

    for i in range(n-2):
        if weight_1 ( r[i,:], r[i+2:,:] )==0: # not i+1
            return 0

    return 1 # No overlaps detected, so weight is one

def weight_1 ( ri, r ):
    """Takes in coordinates of an atom and returns weight = 0 (overlap) or 1 (no overlap).

    Partner coordinate array is supplied.
    """

    import numpy as np

    # This is athermal, based on hard core nonbonded overlaps only
    # The result is either 0 (overlap) or 1 (non overlap)
    # Arithmetically w is the product of all the pair Boltzmann factors
    # Here, each is 0 or 1, which allows us to treat them as Boolean variables
    # but in the more general case they would be real and we would multiply them all together

    # In general, r will be a subset of the complete set of simulation coordinates
    # and none of its rows should be identical to ri

    nj, d = r.shape
    assert d==3, 'Dimension error for r in weight_1'
    assert ri.size==3, 'Dimension error for ri in weight_1'

    if fast:
        rij    = ri - r                # Separation vectors
        rij_sq = np.sum(rij**2,axis=1) # Squared separations
        if np.any(rij_sq<1.0):
            return 0 # Overlap detected
        
        return 1 # No overlap detected

    else:
        for rj in r:
            rij    = ri - rj         # Separation vector
            rij_sq = np.sum(rij**2)  # Squared separation
            if rij_sq<1.0:
                return 0
        return 1

def qcount ( r, q_range ):
    """Returns a count of all square-well interactions within q_range."""

    n, d = r.shape
    assert d==3, 'Dimension error for r in potential'

    q = 0
    
    for i in range(n-2):
        q = q + qcount_1 ( r[i,:], r[i+2:,:], q_range ) # not i+1

    return q

def qcount_1 ( ri, r, q_range ):
    """Returns a count of all square-well interactions within q_range of given atom.

    Partner coordinate array is supplied.
    """

    import numpy as np

    # In general, r will be a subset of the complete set of simulation coordinates
    # and none of its rows should be identical to ri

    nj, d = r.shape
    assert d==3, 'Dimension error for r in weight_1'
    assert ri.size==3, 'Dimension error for ri in weight_1'

    q_range_sq = q_range**2
    
    if fast:
        rij    = ri - r                # Separation vectors
        rij_sq = np.sum(rij**2,axis=1) # Squared separations
        return np.count_nonzero ( rij_sq<q_range_sq )

    else:
        q = 0
        for rj in r:
            rij    = ri - rj         # Separation vector
            rij_sq = np.sum(rij**2)  # Squared separation
            if rij_sq<q_range_sq:
                q = q + 1
        return q

def accept ( s, q_old, q_new, w_old=1.0, w_new=1.0 ):
    """Returns accept/reject decision based on supplied entropy table."""

    import numpy as np
    from maths_module import metropolis

    # This routine is essentially the Metropolis-Hastings formula, generalized
    # to use a specified tabulated function of the energy in the exponent
    # For NVT MC, s(q) is just E(q)/kT = -q/kT since energy is -q
    # For Wang-Landau, s(q) is the entropy function, which is updated through the run
    # In either case, weights (real, but taking positive integer values here) may optionally appear

    w_tol = np.float_(0.5)
    
    # Check that we are within bounds for the entropy table
    assert q_new>=0 and q_new<s.size, 'q_new out of bounds'
    assert q_old>=0 and q_old<s.size, 'q_old out of bounds'

    # Check values of weights
    assert w_old>w_tol, 'Impossible weight error'
    if w_new<w_tol: # Should really have been detected before
        return False
    
    delta = s[q_new] - s[q_old] # Change in entropy function

    # Equivalent to exp(-delta) -> (w_new/w_old)*exp(-delta)
    delta = delta - np.log(w_new/w_old)

    result = metropolis ( delta )
    return result
