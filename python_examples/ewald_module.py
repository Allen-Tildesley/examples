#!/usr/bin/env python3
# ewald_module.f90

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
    assert d==3, 'r dimension error'
    assert n==q.size, 'q dimension error'

    pot = 0.0

    for i,j in combinations(range(n),2):  # Double loop over pairs of atoms
        rij  = r[i,:] - r[j,:]            # Separation vector (box=1 units)
        rij  = rij - np.rint ( rij )      # Periodic boundaries
        rij_mag = np.sqrt(np.sum(rij**2)) # Magnitude of separation
        vij = q[i] * q[j] * erfc ( kappa * rij_mag ) / rij_mag # Screened Coulomb term
        pot = pot + vij

    return pot

first_call = True

def pot_k_ewald ( nk, r, q, kappa ):
    """Returns k-space part of potential energy."""
   
    import numpy as np
    from itertools import product
    global first_call, k_sq_max, kfac

    n,d = r.shape
    assert d==3, 'r dimension error'
    assert n==q.size, 'q dimension error'
  
    twopi = 2.0*np.pi
    twopi_sq = twopi**2
    
    if first_call: # Precalculation of expressions at the start

        b = 1.0 / 4.0 / kappa / kappa
        k_sq_max = nk**2            # Store this in module data
        kfac = np.zeros(k_sq_max+1,dtype=np.float_)
       
        # Lazy triple loop, which over-writes the same values of ksq several times
        # and leaves some locations empty (those which are not the sum of three squared integers)
        # We are only interested in the value of k**2, so skip all negative components
        for kx,ky,kz in product(range(nk+1),repeat=3):  # Triple loop over k vector components

            k_sq = kx**2 + ky**2 + kz**2

            if k_sq <= k_sq_max and k_sq>0: # Test to ensure within range
                kr_sq      = twopi_sq * k_sq           # k**2 in real units
                kfac[k_sq] = twopi * np.exp ( -b * kr_sq ) / kr_sq # Stored expression for later use

        first_call = False

    # Double-check value on later calls
    assert k_sq_max == nk**2, 'nk error'

    eikx = np.zeros((n,nk+1),  dtype=np.complex_) # omits negative k indices
    eiky = np.zeros((n,2*nk+1),dtype=np.complex_) # includes negative k indices
    eikz = np.zeros((n,2*nk+1),dtype=np.complex_) # includes negative k indices

    # Calculate kx, ky, kz = 0, 1 explicitly
    eikx[:,   0] = 1.0 + 0.0j
    eiky[:,nk+0] = 1.0 + 0.0j
    eikz[:,nk+0] = 1.0 + 0.0j

    eikx[:,   1] = np.cos(twopi*r[:,0]) + np.sin(twopi*r[:,0])*1j
    eiky[:,nk+1] = np.cos(twopi*r[:,1]) + np.sin(twopi*r[:,1])*1j
    eikz[:,nk+1] = np.cos(twopi*r[:,2]) + np.sin(twopi*r[:,2])*1j

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

        for ky,kz in product(range(-nk,nk+1),repeat=2):  # Double loop over ky, kz vector components

            k_sq = kx**2 + ky**2 + kz**2

            if k_sq <= k_sq_max and k_sq > 0: # Test to ensure within range

                term = np.sum ( q[:] * eikx[:,kx] * eiky[:,nk+ky] * eikz[:,nk+kz] ) # Sum over all ions
                pot  = pot + factor * kfac[k_sq] * np.real ( np.conj(term)*term )

    # Subtract self part of k-space sum
    pot = pot - kappa * np.sum ( q**2 ) / np.sqrt(np.pi)

    return pot

def pot_k_pm_ewald ( nk, r, q, kappa ):
    """Returns k-space part of potential energy (particle-mesh method)."""
   
    import numpy as np
    from mesh_module import mesh_function, sharpen

    n,d = r.shape
    assert d==3, 'r dimension error'
    assert n==q.size, 'q dimension error'
  
    twopi  = 2.0*np.pi
    fourpi = 4.0*np.pi
    dk     = twopi # k-spacing for box=1
    dr     = 1.0 / (2*nk) # r-spacing
    sc     = 2*nk # Use nk to determine mesh dimension
    
    # Set up indices in wraparound convention
    imesh      = np.arange(sc)
    imesh[nk:] = imesh[nk:] - sc
    kmesh      = (dk*imesh).tolist()

    # Assign charge density to complex array ready for FFT
    # Assume r in unit box with range (-0.5,0.5)
    # Mesh function expects coordinates in range 0..1 so we add 0.5
    fft_inp = mesh_function ( r+0.5, q, 2*nk ).astype(np.complex_)

    # Forward FFT incorporating scaling by number of grid points
    fft_out = np.fft.fftn(fft_inp) / (sc**3)

    # Square modulus of charge density
    rho_sq = np.real ( fft_out*np.conjugate(fft_out) )

    pot = 0.0

    # Triple loop over xyz grid points in k-space
    for ix, kx in enumerate(kmesh):
        fx = sharpen( 0.5*kx*dr )
        for iy, ky in enumerate(kmesh):
            fxy = fx*sharpen( 0.5*ky*dr )
            for iz, kz in enumerate(kmesh):
                fxyz = fxy*sharpen( 0.5*kz*dr )

                if ix==0 and iy==0 and iz==0: # Skip zero wave vector
                    continue
                k_sq = kx**2 + ky**2 + kz**2  # Squared magnitude of wave vector
                g    = (fourpi/k_sq) * np.exp ( -k_sq / (4.0*kappa**2) ) # Uncorrected influence function
                g    = g * fxyz # Apply simple correction
                pot  = pot + g * rho_sq[ix,iy,iz]

    # Divide by 2; box volume is 1
    pot = pot / 2.0

    # Subtract self part of k-space sum
    pot = pot - kappa * np.sum ( q**2 ) / np.sqrt(np.pi)

    return pot
