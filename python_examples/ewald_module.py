#!/usr/bin/env python3
# ewald_module.f90

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

"""r-space and k-space parts of Ewald sum for ions."""

  # References:
  # Woodcock and Singer, Trans. Faraday Soc. 67, 12, 1971.
  # de Leeuw et al., Proc. Roy. Soc. A 373, 27, 1980.
  # Heyes, J. Chem. Phys. 74, 1924, 1981.
  # see also Fincham, mdions, CCP5 program library.

  # The self term is subtracted from the k-space contribution
  # The surface term for simulations in vacuum is not included
  # A cubic box and unit box length are assumed throughout
  # No special lists are used

def pot_r_ewald ( r, q, kappa ):
    """Returns r-space part of potential energy."""
   
    import numpy as np
    from itertools import combinations
    from scipy.special import erfc
    
    n,d = r.shape
    assert d==3, 'r dimension error in pot_r_ewald'
    assert n==q.size, 'q dimension error in pot_r_ewald'

    pot = 0.0

    for i,j in combinations(range(n),2):  # Double loop over pairs of atoms
        rij  = r[i,:] - r[j,:]            # Separation vector (box=1 units)
        rij  = rij - np.rint ( rij )      # Periodic boundaries
        rij_mag = np.sqrt(np.sum(rij**2)) # Magnitude of separation
        vij = q[i] * q[j] * erfc ( kappa * rij_mag ) / rij_mag # Screened Coulomb term
        pot = pot + vij

    return pot
 
def pot_k_ewald ( nk, r, q, kappa ):
    """Returns k-space part of potential energy."""
   
    import numpy as np
    from itertools import combinations

    twopi = 2.0*np.pi
    twopi_sq = twopi**2
    
    if pot_k_ewald.first_call: # Precalculation of expressions at the start

        b = 1.0 / 4.0 / kappa / kappa
        k_sq_max = nk**2            # Store this in module data
        kfac = np.zeros(k_sq_max+1,dtype=np.float_)
       
        # Lazy triple loop, which over-writes the same values of ksq several times
        # and leaves some locations empty (those which are not the sum of three squared integers)
        # We are only interested in the value of k**2, so skip all negative components
        for kx,ky,kz in combinations(range(nk+1),3):  # Triple loop over k vector components

            k_sq = kx**2 + ky**2 + kz**2

            if k_sq <= k_sq_max and k_sq>0: # Test to ensure within range
                kr_sq      = twopi_sq * k_sq           # k**2 in real units
                kfac[k_sq] = twopi * np.exp ( -b * kr_sq ) / kr_sq # Stored expression for later use

        pot_k_ewald.first_call = False

    # Double-check value on later calls
    assert k_sq_max == nk**2, 'nk error'

    eikx = np.zeros((n,nk+1),dtype=np.complex_)   # omits negative k indices
    eiky = np.zeros((n,2*nk+1),dtype=np.complex_) # includes negative k indices
    eiky = np.zeros((n,2*nk+1),dtype=np.complex_) # includes negative k indices

    # Calculate kx, ky, kz = 0, 1 explicitly
    eikx[:,   0] = 1.0 + 0.0j
    eiky[:,nk+0] = 1.0 + 0.0j
    eikz[:,nk+0] = 1.0 + 0.0j

    eikx[:,   1] = np.cos(twopi*r[:,1]) + np.sin(twopi*r[:,1])*1j
    eiky[:,nk+1] = np.cos(twopi*r[:,2]) + np.sin(twopi*r[:,2])*1j
    eikz[:,nk+1] = np.cos(twopi*r[:,3]) + np.sin(twopi*r[:,3])*1j

    # Calculate remaining positive kx, ky and kz by recurrence
    for k in range(2,nk+1):
        eikx[:,   k] = eikx[:,   k-1] * eikx[:,   1]
        eiky[:,nk+k] = eiky[:,nk+k-1] * eiky[:,nk+1]
        eikz[:,nk+k] = eikz[:,nk+k-1] * eikz[:,nk+1]

    # Negative k values are complex conjugates of positive ones
    # We do not need negative values of kx
    eiky[:,0:nk] = np.conj ( eiky[:,2*nk:nk:-1] )
    eikz[:,0:nk] = np.conj ( eikz[:,2*nk:nk:-1] )

    pot = 0.0

    for kx in range(nk+1): # Outer loop over non-negative kx

        factor = 1.0 if kx==0 else 2.0 # Accounts for skipping negative kx

        for ky,kz in combinations(range(2*nk+1),2):  # Double loop over ky, kz vector components

            k_sq = kx**2 + (ky-nk)**2 + (kz-nk)**2

            if k_sq <= k_sq_max and k_sq > 0: # Test to ensure within range

                term = np.sum ( q[:] * eikx[:,kx] * eiky[:,ky] * eikz[:,kz] ) # Sum over all ions
                pot  = pot + factor * kfac[k_sq] * np.real ( np.conj(term)*term )

    # Subtract self part of k-space sum
    pot = pot - kappa * np.sum ( q**2 ) / np.sqrt(np.pi)

    return pot

pot_k_ewald.firstcall=True
