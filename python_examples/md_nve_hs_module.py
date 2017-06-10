#!/usr/bin/env python3
# md_nve_hs_module.py

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

"""Collisions and overlap for MD of hard spheres."""

def introduction():
    """Prints out introductory statements at start of run."""
    
    print('Hard sphere potential')
    print('Diameter, sigma = 1')
    print('Energy, kT = 1')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def update ( i, box, r, v ):
    """Looks for collisions with i, for particles j>i."""

    import numpy as np
    
    # This routine loops over atoms in r seeking shortest collision time with i
    # We are only interested in j>i, and collision info will be stored under index i
    # r and v just contain slices j>=i so r[0,:] is position of i, r[1:,:] is all j>i
    # We do all this with NumPy routines

    huge = 1.e9

    rij   = r[0,:] - r[1:,:]           # All the separation vectors
    rij   = rij - np.rint ( rij )        # Periodic boundaries, box=1 units
    rij   = rij * box                  # Now in sigma=1 units
    vij   = v[0,:] - v[1:,:]           # All the relative velocities
    bij   = np.sum ( rij*vij, axis=1 ) # All the dot products
    rijsq = np.sum ( rij**2, axis=1 )  # All the squared distances
    vijsq = np.sum ( vij**2, axis=1 )  # All the squared speeds

    discr = bij ** 2 - vijsq * ( rijsq - 1.0 ) # All the discriminants, sigma**2 = 1.0

    # Collision times calculated only where collisions are possible
    mask = np.logical_and (bij<0.0, discr>0.0)
    tij = np.where ( mask,  (-bij-np.sqrt(np.fabs(discr)))/vijsq, huge )

    k = np.argmin(tij)      # Location of minimum time in array tij
    return tij[k], k+i+1 # Return partner, remembering that tij[0] corresponds to i+1

def dndate ( i, box, r, v, coltime, partner ):
    """Looks for collisions with r[i], for particles j<i."""

    import numpy as np
    
    # This routine loops over atoms in r seeking shortest collision time with i
    # We are only interested in j<i, and collision info will be stored under index j
    # r and v contain slices up to and including i
    # coltime and partner correspond to j<i and
    # we only update the value corresponding to any earlier collision times 
    # We do all this with NumPy routines

    huge = 1.e9
    
    rij   = r[i] - r[:i,:]               # All the separation vectors
    rij   = rij - np.rint ( rij )        # Periodic boundaries, box=1 units
    rij   = rij * box                  # Now in sigma=1 units
    vij   = v[i] - v[:i,:]               # All the relative velocities
    bij   = np.sum ( rij*vij, axis=1 ) # All the dot products
    rijsq = np.sum ( rij**2, axis=1 ) # All the squared distances
    vijsq = np.sum ( vij**2, axis=1 ) # All the squared speeds

    discr = bij ** 2 - vijsq * ( rijsq - 1.0 ) # All the discriminants, sigma**2 = 1.0

    # Collision times calculated only where collisions are possible
    mask    = np.logical_and ( bij<0.0, discr>0.0 )
    tij     = np.where ( mask,  (-bij-np.sqrt(np.fabs(discr)))/vijsq, huge )
    mask    = tij < coltime
    coltime = np.where ( mask, tij, coltime )
    partner = np.where ( mask, i, partner )

    return coltime, partner

def overlap ( box, r ):
    import numpy as np
    from itertools import combinations
    
    box_sq  = box**2

    for ri, rj in combinations(r,2):
        rij = ri - rj
        rij = rij - np.rint ( rij )
        rij_sq = np.sum ( rij**2 )  # squared distance
        rij_sq = rij_sq * box_sq # now in sigma=1 units

        if rij_sq < 1.0:
            rij_mag = np.sqrt(rij_sq)
            print('Warning, overlap detected')
            return True

    return False
  
def collide ( ri, vi, rj, vj, box ):
    """Implements collision dynamics, updating the velocities."""

    import numpy as np
    
    # The colliding pair (i,j) is assumed to be in contact already
    
    rij = ri - rj
    rij = rij - np.rint ( rij ) # Separation vector
    rij = rij * box             # Now in sigma=1 units
    vij = vi - vj               # Relative velocity

    factor = np.dot ( rij, vij )
    vij    = -factor * rij

    vi = vi + vij
    vj = vj - vij
    virial = np.dot ( vij, rij ) / 3.0
    return vi, vj, virial
