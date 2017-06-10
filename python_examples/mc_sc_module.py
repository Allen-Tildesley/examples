#!/usr/bin/env python3
# mc_sc_module.py

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

"""Overlap and move routines for MC simulation, hard spherocylinders."""

fast   = True # Change this to replace NumPy overlap evaluation with slower Python
length = 5.0 # Cylinder length L (in units where D=1) used throughout the module

def introduction():
    """Prints out introductory statements at start of run."""

    import numpy as np
    
    vmol   = np.pi * ( 0.25*length + 1/6 ) # Spherocylinder volume

    print('Hard spherocylinder potential')
    print("{:40}{:15.6f}".format('Spherocylinder L/D ratio',   length))
    print("{:40}{:15.6f}".format('Spherocylinder volume/D**3', vmol))
    print('Diameter, D = 1')
    print('Energy, kT = 1')
    if fast:
        print('Fast NumPy overlap routine')
    else:
        print('Slow Python overlap routine')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def overlap ( box, r, e ):
    """Takes in box and coordinate & orientation arrays, and signals any overlap."""

    # Actual calculation is performed by function overlap_1

    n, d = r.shape
    assert d==3, 'Dimension error for r in overlap'
    assert d==e.shape[1], 'Dimension error for e in overlap'
    assert n==e.shape[0], "{}{:d}{:d}".format('Dimension error for e in overlap',n,e.shape[0])

    for i in range(n-1):
        if overlap_1 ( r[i,:], e[i,:], box, r[i+1:,:], e[i+1:,:] ):
            return True # Immediate return on detection of overlap

    return False

def overlap_1 ( ri, ei, box, r, e ):
    """Takes in coordinates and orientations of a molecules and signals any overlap.

    Values of box and partner coordinate array are supplied.
    """

    import numpy as np

    # In general, r will be a subset of the complete set of simulation coordinates
    # and none of its rows should be identical to ri

    # It is assumed that positions are in units where box = 1

    nj, d = r.shape
    assert d==3, 'Dimension error for r in overlap_1'
    assert d==e.shape[1], 'Dimension error for e in overlap_1'
    assert nj==e.shape[0], "{}{:d}{:d}".format('Dimension error for e in overlap_1',nj,e.shape[0])
    assert ri.size==3, 'Dimension error for ri in overlap_1'
    assert ei.size==3, 'Dimension error for ei in overlap_1'

    range = 1.0 + length
    assert range <= box/2.0, "{}{:15.6f}{:15.6f}".format('Box too small', box, range)
    range_box_sq = ( range / box ) ** 2 # Squared range in box=1 units
    box_sq       = box**2               # Squared box length
    
    if fast:
        rij = ri - r                    # Get all separation vectors from partners
        rij = rij - np.rint(rij)        # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2,axis=1)  # Squared separations
        rij_sq = rij_sq * box_sq        # Now in D=1 units
        rij    = rij    * box           # Now in D=1 units
        rei    = np.sum ( rij*ei, axis=1 ) # All dot products
        rej    = np.sum ( rij*e,  axis=1 ) # All dot products
        eij    = np.sum ( ei *e,  axis=1 ) # All dot products
        sij_sq = all_dist_sq ( rij_sq, rei, rej, eij ) # Squared distance between line segments
        return np.any(sij_sq<1.0)

    # Otherwise use slow method
    for j,rj in enumerate(r):
        rij = ri - rj             # Separation vector
        rij = rij - np.rint(rij)  # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2)   # Squared separation
        if rij_sq > range_box_sq: # Check no possibility of overlap
            continue
        rij_sq = rij_sq * box_sq # Now in D=1 units
        rij    = rij    * box    # Now in D=1 units
        rei    = np.dot ( rij, ei     )
        rej    = np.dot ( rij, e[j,:] )
        eij    = np.dot ( ei,  e[j,:] )

        sij_sq = dist_sq ( rij_sq, rei, rej, eij ) # Squared distance between line segments
        if sij_sq<1.0:
            return True # Overlap detected, return immediately

    return False

def n_overlap ( box, r, e ):
    """Takes in box and coordinate and orientation arrays, and counts overlaps."""

    # This routine is used in the calculation of pressure
    # Actual calculation is performed by function n_overlap_1

    n, d = r.shape
    assert d==3, 'Dimension error for r in n_overlap'
    assert d==e.shape[1], 'Dimension error for e in n_overlap'
    assert n==e.shape[0], "{}{:d}{:d}".format('Dimension error for e in n_overlap',n,e.shape[0])

    n_ovr = 0
    for i in range(n-1):
        n_ovr = n_ovr + n_overlap_1 ( r[i,:], e[i,:], box, r[i+1:,:], e[i+1:,:] )

    return n_ovr

def n_overlap_1 ( ri, ei, box, r, e ):
    """Takes in coordinates and orientations of a molecule and counts overlaps.

    Values of box and partner coordinate array are supplied.
    Fast or slow algorithm selected.
    """

    import numpy as np

    # In general, r will be a subset of the complete set of simulation coordinates
    # and none of its rows should be identical to ri

    # It is assumed that positions are in units where box = 1

    nj, d = r.shape
    assert d==3, 'Dimension error for r in n_overlap_1'
    assert d==e.shape[1], 'Dimension error for e in n_overlap_1'
    assert nj==e.shape[0], "{}{:d}{:d}".format('Dimension error for e in n_overlap_1',nj,e.shape[0])
    assert ri.size==3, 'Dimension error for ri in n_overlap_1'
    assert ei.size==3, 'Dimension error for ei in n_overlap_1'

    range = 1.0 + length
    assert range <= box/2.0, "{}{:15.6f}{:15.6f}".format('Box too small', box, range)
    range_box_sq = ( range / box ) ** 2 # Squared range in box=1 units
    box_sq       = box**2               # Squared box length

    if fast:
        rij = ri - r                    # Get all separation vectors from partners
        rij = rij - np.rint(rij)        # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2,axis=1)  # Squared separations
        rij_sq = rij_sq * box_sq        # Now in D=1 units
        rij    = rij    * box           # Now in D=1 units
        rei    = np.sum ( rij*ei, axis=1 ) # All dot products
        rej    = np.sum ( rij*e,  axis=1 ) # All dot products
        eij    = np.sum ( ei *e,  axis=1 ) # All dot products
        sij_sq = all_dist_sq ( rij_sq, rei, rej, eij ) # Squared distance between line segments
        return np.count_nonzero(sij_sq<1.0)

    # Otherwise use slow method
    n_ovr = 0
    for j,rj in enumerate(r):
        rij = ri - rj             # Separation vector
        rij = rij - np.rint(rij)  # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2)   # Squared separation
        if rij_sq > range_box_sq: # Check no possibility of overlap
            continue
        rij_sq = rij_sq * box_sq # Now in D=1 units
        rij    = rij    * box    # Now in D=1 units
        rei    = np.dot ( rij, ei     )
        rej    = np.dot ( rij, e[j,:] )
        eij    = np.dot ( ei,  e[j,:] )

        sij_sq = dist_sq ( rij_sq, rei, rej, eij ) # Squared distance between line segments
        if sij_sq<1.0:
            n_ovr = n_ovr + 1

    return n_ovr

def dist_sq ( rij_sq, rei, rej, eij ):
    import numpy as np

    tol = 1.0e-6
    ell2   = length/2.0   # Half length
    sin_sq = 1.0 - eij**2 # Squared sine of angle between line segments

    if sin_sq < tol: 
        ci = -rei
        cj =  rej
    else:
        ci = ( - rei + eij * rej ) / sin_sq
        cj = (   rej - eij * rei ) / sin_sq

    ai = np.fabs ( ci )
    aj = np.fabs ( cj )
    if ai > ell2:
        ci = ell2*np.sign(ci)
    if aj > ell2:
        cj = ell2*np.sign(cj)

    if ai > aj:
        cj =  rej + ci * eij
    else:
        ci = -rei + cj * eij

    ai = np.fabs ( ci )
    aj = np.fabs ( cj )
    if ai > ell2:
        ci = ell2*np.sign(ci)
    if aj > ell2:
        cj = ell2*np.sign(cj)

    di =  2.0 * rei + ci - cj * eij
    dj = -2.0 * rej + cj - ci * eij

    return rij_sq + ci * di + cj * dj # Squared distance between line segments
    
def all_dist_sq ( rij_sq, rei, rej, eij ):
    import numpy as np
    
    tol = 1.0e-6
    ell2   = length/2.0   # Half length
    sin_sq = 1.0 - eij**2 # Squared sines of angles between line segments
    mask = sin_sq>tol
    ci = np.where ( mask, (-rei+eij*rej)/sin_sq, -rei )
    cj = np.where ( mask, ( rej-eij*rei)/sin_sq,  rej )

    ai = np.fabs ( ci )
    aj = np.fabs ( cj )
    ci = np.where ( ai>ell2, ell2*np.sign(ci), ci )
    cj = np.where ( aj>ell2, ell2*np.sign(cj), cj )
    mask = ai>aj
    cj = np.where ( mask, rej+ci*eij, cj )
    mask = np.logical_not(mask)
    ci = np.where ( mask, -rei+cj*eij, ci )

    ai = np.fabs ( ci )
    aj = np.fabs ( cj )
    ci = np.where ( ai>ell2, ell2*np.sign(ci), ci )
    cj = np.where ( aj>ell2, ell2*np.sign(cj), cj )

    di =  2.0 * rei + ci - cj * eij
    dj = -2.0 * rej + cj - ci * eij

    return rij_sq + ci * di + cj * dj # Squared distances between line segments
