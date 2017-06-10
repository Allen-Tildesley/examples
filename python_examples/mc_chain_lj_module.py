#!/usr/bin/env python3
# mc_chain_lj_module.py

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

"""Monte Carlo, single chain, LJ atoms."""

fast = True # Change this to replace NumPy potential evaluation with slower Python

class PotentialType:
    """A composite variable for interactions."""

    def __init__(self, pot, ovr):
        self.pot = pot # the potential energy
        self.ovr = ovr # a flag indicating overlap (i.e. pot too high to use)

    def __add__(self, other):
        pot = self.pot +  other.pot
        ovr = self.ovr or other.ovr
        return PotentialType(pot,ovr)

    def __sub__(self, other):
        pot = self.pot -  other.pot
        ovr = self.ovr or other.ovr # This is meaningless, but inconsequential
        return PotentialType(pot,ovr)

def introduction():
    """Prints out introductory statements at start of run."""

    # This model, specifically its collapse behaviour, is discussed in detail by
    # F Calvo, JPK Doye, DJ Wales, J Chem Phys 116, 2642 (2002)
    # A similar model, with rigid bond lengths, is discussed by
    # A Irback, E Sandelin, J Chem Phys 110, 12256 (1999)
    # Both types of model are also studied by
    # P Grassberger, R Hegger, J Phys Cond Matt 7, 3089 (1995)

    print('LJ chain, no cutoff, no shift, no periodic box')
    print('Diameter, sigma = 1')
    print('Well depth, epsilon = 1')
    print('Harmonic spring bond potential')
    if fast:
        print('Fast NumPy potential routine')
    else:
        print('Slow Python potential routine')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def regrow ( temperature, m_max, k_max, bond, k_spring, r ):
    """Carries out single regrowth move, returning new r and indicator of success."""

    # A short sequence of m atoms (m<=m_max) is deleted and regrown in the CBMC manner
    # We randomly select which end of the chain to apply each of these operations to
    # At each stage, k_max different atom positions are tried
    # Function random_bond selects bond lengths according to the internal (harmonic) potential
    # Rosenbluth weights are computed using the external (nonbonded) potential
    # Acceptance/rejection is determined using these weights

    # r_old and r_new are used as working arrays

    import numpy as np
    from maths_module import random_vector

    w_tol = 1.e-10 # Min weight tolerance
    n, d = r.shape
    assert d==3, 'Dimension error in regrow'
    r_try = np.empty((k_max,3),dtype=np.float_)
    w     = np.empty(k_max,dtype=np.float_)

    std   = np.sqrt(temperature/k_spring) # Spring bond standard deviation
    d_max = 3.0*std                       # Impose a limit on variation, say 3*std
    assert d_max<0.75*bond, 'Spring bond strength error'
    d_max = d_max + bond # This is the actual max d allowed

    r_old = np.copy(r) # Store copy of r

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

    w_new = 1.0
    for i in range(n-m,n): # Loop to regrow last m atoms, computing new weight
        for k in range(k_max): # Loop over k_max tries
            d          = random_bond ( bond, std, d_max )      # Generate random bond length around d=bond
            r_try[k,:] = r[i-1,:] + d * random_vector()        # Trial position in random direction from i-1
            partial    = potential_1 ( r_try[k,:], r[:i-1,:] ) # Nonbonded interactions with earlier atoms (not i-1)
            w[k]       = 0.0 if partial.ovr else np.exp(-partial.pot/temperature) # Weight for this try
        w_sum = np.sum(w)
        if w_sum<w_tol: # Early exit if this happens at any stage
            return r_old, False
        w = w / w_sum
        k = np.random.choice(k_max,p=w) # Pick winning try according to weights
        r[i,:] = r_try[k,:]             # Store winning position
        w_new = w_new * w_sum           # Accumulate total weight

    if w_new<w_tol: # Exit if this happens
        return r_old, False

    r_new = np.copy(r) # Store new configuration

    # END OF PART 1: NEW CONFIGURATION AND WEIGHT ARE COMPLETE

    # PART 2: RECONSTRUCT OLD CONFIGURATION WITH OLD WEIGHT

    if c==0 or c==1: # Remove and hence reconstruct at end
        r[:,:] = r_old[:,:] # Copy all n atoms
    else: # Remove and reconstruct at start
        r[:,:] = r_old[::-1,:] # Copy and reverse all n atoms

    w_old = 1.0
    for i in range(n-m,n):

        # Old position and weight are stored as try 0
        r_try[0,:] = r[i,:]
        partial    = potential_1 ( r_try[0,:], r[:i-1,:] ) # Nonbonded energy with earlier atoms (not i-1)
        w[0]       = 0.0 if partial.ovr else np.exp(-partial.pot/temperature) # Weight for this try

        # Remaining tries only required to compute weight
        for k in range(1,k_max): # Loop over k_max-1 other tries
            d          = random_bond ( bond, std, d_max )      # Generate random bond length around d=bond
            r_try[k,:] = r[i-1,:] + d * random_vector()        # Trial position in random direction from i-1
            partial    = potential_1 ( r_try[k,:], r[:i-1,:] ) # Nonbonded interactions with earlier atoms (not i-1)
            w[k]       = 0.0 if partial.ovr else np.exp(-partial.pot/temperature) # Weight for this try

        w_sum = np.sum(w)
        r[i,:] = r_try[0,:]   # Restore winning position (always the original one)
        w_old = w_old * w_sum # Accumulate total weight

    assert w_old>w_tol, 'Old weight error'

    # END OF PART 2: OLD CONFIGURATION AND WEIGHT ARE COMPLETE

    # Choose either old or new configuration according to weight
    # All non-bonded Boltzmann factors are incorporated into the weights
    # All spring-bond Boltzmann factors are included in the selection of d

    zeta = np.random.rand()
    if zeta < (w_new/w_old):
        return r_new, True
    else:
        return r_old, False

def potential ( r ):
    """Takes in coordinate array, and calculates total potential etc.

    The results are returned as total, a PotentialType variable."""

    # total.pot is the nonbonded potential energy for whole system
    # total.ovr is a flag indicating overlap (potential too high) to avoid overflow
    # If this flag is True, the value of total.pot should not be used
    # Actual calculation is performed by function potential_1

    n, d = r.shape
    assert d==3, 'Dimension error for r in potential'

    total = PotentialType ( pot=0.0, ovr=False ) # Initialize

    for i in range(n-2):
        partial = potential_1 ( r[i,:], r[i+2:,:] ) # not i+1
        if partial.ovr:
            total.ovr = True
            break
        total = total + partial

    return total

def potential_1 ( ri, r ):
    """Takes in coordinates of an atom and calculates its interactions.

    Partner coordinate array is supplied.
    """

    import numpy as np

    # partial.pot is the potential energy of atom ri with a set of other atoms
    # partial.ovr is a flag indicating overlap (potential too high) to avoid overflow
    # If this is True, the value of partial.pot should not be used
    # In general, r will be a subset of the complete set of simulation coordinates
    # and none of its rows should be identical to ri

    # Coordinates are assumed to be in LJ units where sigma = 1, no box, no periodic boundaries
    # Results are in LJ units where sigma = 1, epsilon = 1

    nj, d = r.shape
    assert d==3, 'Dimension error for r in potential_1'
    assert ri.size==3, 'Dimension error for ri in potential_1'

    sr2_ovr      = 1.77 # Overlap threshold (pot > 100)

    # Initialize
    partial = PotentialType ( pot=0.0, ovr=False )

    if fast:
        rij    = ri - r                # Separation vectors
        rij_sq = np.sum(rij**2,axis=1) # Squared separations
        sr2    = 1.0 / rij_sq          # (sigma/rij)**2
        ovr    = sr2 > sr2_ovr         # Overlap if too close

        if any(ovr):
            partial.ovr=True
            return partial
        
        sr6  = sr2 ** 3
        sr12 = sr6 ** 2
        pot  = sr12 - sr6
        
        partial = PotentialType ( pot=np.sum(pot), ovr=False )

    else:
        for rj in r:
            rij    = ri - rj         # Separation vector
            rij_sq = np.sum(rij**2)  # Squared separation
            sr2    = 1.0 / rij_sq    # (sigma/rij)**2
            ovr    = sr2 > sr2_ovr   # Overlap if too close

            if ovr:
                partial.ovr=True
                return partial

            sr6  = sr2 ** 3
            sr12 = sr6 ** 2
            pot  = sr12 - sr6

            partial = partial + PotentialType ( pot=pot, ovr=ovr )

    # Multiply results by numerical factors
    partial.pot = partial.pot * 4.0 # 4*epsilon

    return partial

def random_bond ( b, std, d_max):
    """Returns random bond length sampled from weighted Gaussian distribution."""

    # Uses von Neumann's rejection method to sample (d**2)*exp(-0.5*(d-b)**2/std**2)
    # The sampled distribution is the same but with d in the prefactor (which arises from
    # the Jacobian in 3D) replaced by the constant d_max, which makes it a Gaussian function
    # Hence, the range must be restricted to d<d_max, for the rejection method to work
    # It will be reasonably efficient provided std is small compared with b
    # This is essentially the same method as an example in
    # Understanding Molecular Simulation by D Frenkel and B Smit

    import numpy as np

    while True:
        d = b + std*np.random.randn()
        if d<0.0 or d>d_max:
            continue
        zeta = np.random.rand()
        if zeta <= (d/d_max)**2: # Compare with ratio of distributions to accept
            break

    return d

def spring_pot ( b, k, r ):
    """Returns internal spring potential energy."""

    import numpy as np

    rij = r[:-1,:]-r[1:,:] # Get all nearest neighbour vectors
    d = np.sqrt(np.sum(rij**2,axis=1)) # All nearest neighbour distances
    pot = 0.5*k*np.sum ( (d-b)**2 )

    return pot
