#!/usr/bin/env python3
# mc_poly_lj_module.py

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

"""Routines for MC simulation, polyatomic molecule, LJ atoms."""

import numpy as np

fast = True # Change this to replace NumPy potential evaluation with slower Python

# Bond vectors in body-fixed frame
# Isosceles triangle, 3 sites, with unit bond length and bond angle alpha
# which we set to 75 degrees here
alpha  = 75.0 * np.pi / 180.0
alpha2 = alpha / 2.0
na = 3
db = np.array([[-np.sin(alpha2), 0.0,   -np.cos(alpha2)/3.0],
               [0.0,            0.0, 2.0*np.cos(alpha2)/3.0],
               [np.sin(alpha2), 0.0,    -np.cos(alpha2)/3.0]], dtype=np.float_)
diameter = 2.0 * np.sqrt ( np.max ( np.sum(db**2,axis=1) ) ) # Molecular diameter

# Cutoff distance and force-shift parameters (all private) chosen as per the reference:
# S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)
r_cut    = 2.612 # in sigma=1 units, where r_cut = 1.2616 nm, sigma = 0.483 nm
sr_cut   = 1.0/r_cut
sr_cut6  = sr_cut**6
sr_cut12 = sr_cut6**2
lambda1  = 4.0*(7.0*sr_cut6-13.0*sr_cut12)
lambda2  = -24.0*(sr_cut6-2.0*sr_cut12)*sr_cut

class PotentialType:
    """A composite variable for interactions."""

    def __init__(self, pot, vir, ovr):
        self.pot = pot # the potential energy
        self.vir = vir # the virial
        self.ovr = ovr # a flag indicating overlap (i.e. pot too high to use)

    def __add__(self, other):
        pot = self.pot +  other.pot
        vir = self.vir +  other.vir
        ovr = self.ovr or other.ovr
        return PotentialType(pot,vir,ovr)

    def __sub__(self, other):
        pot = self.pot -  other.pot
        vir = self.vir -  other.vir
        ovr = self.ovr or other.ovr # This is meaningless, but inconsequential
        return PotentialType(pot,vir,ovr)

def introduction():
    """Prints out introductory statements at start of run."""

    print('Lennard-Jones potential')
    print('Cut-and-force-shifted')
    print('Diameter, sigma = 1')
    print('Well depth, epsilon = 1')
    if fast:
        print('Fast NumPy potential routine')
    else:
        print('Slow Python potential routine')

    print( "{:40}{:15d}".format('Number of atoms per molecule', na) )
    for i, b in enumerate(db):
        print( "{}{:2d}{:15.6f}{:15.6f}{:15.6f}".format('Body-fixed atom vector',i,*b))
    print( "{:40}{:15.6f}".format('Molecular diameter', diameter) )

    print( "{:40}{:15.6f}".format('r_cut', r_cut) )
    print( "{:40}{:15.6f}".format('Force-shift lambda1', lambda1) )
    print( "{:40}{:15.6f}".format('Force-shift lambda2', lambda2) )

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def potential ( box, r, d ):
    """Takes in box and r & d arrays, and calculates total potential etc.

    The results are returned as total, a PotentialType variable.
    """
    # Actual calculation performed by function potential_1

    n, ndim = r.shape
    assert ndim==3, 'Dimension error for r'
    nn, nna, ndim = d.shape
    assert nna==na and ndim==3, 'Dimension error for d'
    assert n==nn, 'Dimension mismatch for r and d'

    total = PotentialType ( pot=0.0, vir=0.0, ovr=False )

    for i in range(n-1):
        partial = potential_1 ( r[i,:], d[i,:,:], box, r[i+1:,:], d[i+1:,:,:] )
        if partial.ovr:
            total.ovr = True
            break
        total = total + partial

    return total

def potential_1 ( ri, di, box, r, d ):
    """Takes in r & d of a molecule and calculates its interactions.

    Values of box and partner coordinate arrays are supplied.
    The results are returned as partial, a PotentialType variable.
    """

    import numpy as np

    # partial.pot is the nonbonded cut (not shifted) potential energy of atom ri with a set of other atoms
    # partial.vir is the corresponding virial of atom ri
    # partial.ovr is a flag indicating overlap (potential too high) to avoid overflow
    # If this is True, the values of partial.pot etc should not be used
    # In general, r & d will be subsets of the complete set of simulation coordinates & bond vectors
    # and none of their rows should be identical to ri, di

    # It is assumed that r has been divided by box
    # Results are in LJ units where sigma = 1, epsilon = 1
    # Note that this is the force-shifted LJ potential with a linear smoothing term
    # S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)

    nj, ndim = r.shape
    assert ndim==3, 'Dimension error for r'
    nnj, nna, ndim = d.shape
    assert nna==na and ndim==3, 'Dimension error for d'
    assert nj==nnj, 'Dimension mismatch between r and d'
    assert ri.size==3, 'Dimension error for ri'
    nna, ndim = di.shape
    assert nna==na and ndim==3, 'Dimension error for di'

    sr2_ovr       = 1.77 # Overlap threshold (pot > 100)
    rm_cut_box    = ( r_cut + diameter ) / box # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box**2              # squared
    assert rm_cut_box<0.5, 'rm_cut/box too large'
    r_cut_sq = r_cut ** 2

    partial = PotentialType ( pot=0.0, vir=0.0, ovr=False ) # Initialize

    if fast:
        rij = ri - r             # Get all separation vectors from partners
        rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
        rij = rij * box          # Now in sigma=1 units
        for a in range(na):
            for b in range(na):
                rab      = rij + di[a,:] - d[:,b,:] # All atom-atom vectors for given a and b
                rab_sq   = np.sum(rab**2,axis=1)    # Squared separations
                in_range = rab_sq < r_cut_sq        # Set flags for within cutoff
                sr2      = 1.0 / rab_sq             # (sigma/rab)**2
                ovr      = sr2 > sr2_ovr            # Set flags for any overlaps
                if np.any(ovr):
                    partial.ovr = True
                    return partial
                rmag  = np.sqrt(rab_sq)
                sr6  = sr2 ** 3
                sr12 = sr6 ** 2
                pot   = np.where ( in_range,
                        4.0*(sr12-sr6) + lambda1 + lambda2*rmag, 0.0 ) # LJ atom-atom pair potentials (force-shifted)
                virab = np.where ( in_range,
                        24.0*(2.0*sr12-sr6) - lambda2*rmag, 0.0 ) # LJ atom-atom pair virials
                fab   = virab * sr2
                fab   = rab * fab[:,np.newaxis] # LJ atom-atom pair forces
                partial = partial + PotentialType ( pot=np.sum(pot), vir=np.sum(rij*fab), ovr=False )

    else:
        for j, rj in enumerate(r):
            rij = ri - rj            # Separation vector
            rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
            rij_sq = np.sum(rij**2)  # Squared separation
            if rij_sq < rm_cut_box_sq: # Check within cutoff
                rij = rij * box # Now in sigma=1 units
                for a in range(na):
                    for b in range(na):
                        rab    = rij + di[a,:] - d[j,b,:] # Atom-atom vector, sigma=1 units
                        rab_sq = np.sum ( rab**2 )        # Squared atom-atom separation, sigma=1 units

                        if rab_sq < r_cut_sq: #     Test within potential cutoff 
                            sr2    = 1.0 / rab_sq    # (sigma/rab)**2
                            ovr    = sr2 > sr2_ovr   # Overlap if too close
                            if ovr:
                                partial.ovr=True
                                return partial
                            rmag  = np.sqrt(rab_sq)
                            sr6   = sr2 ** 3
                            sr12  = sr6 ** 2
                            pot   = 4.0*(sr12-sr6) + lambda1 + lambda2*rmag # LJ atom-atom pair potential (force-shifted)
                            virab = 24.0*(2.0*sr12-sr6) - lambda2*rmag      # LJ atom-atom pair virial
                            fab   = rab * virab * sr2                       # LJ atom-atom pair force
                            partial = partial + PotentialType ( pot=pot, vir=np.sum(rij*fab), ovr=ovr )

    # Include numerical factors
    partial.vir = partial.vir / 3.0 # Divide virial by 3

    return partial
