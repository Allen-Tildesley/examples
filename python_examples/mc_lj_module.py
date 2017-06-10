#!/usr/bin/env python3
# mc_lj_module.py

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

"""Energy and move routines for MC simulation, LJ potential."""

fast = True # Change this to replace NumPy potential evaluation with slower Python

class PotentialType:
    """A composite variable for interactions."""

    def __init__(self, pot, vir, lap, ovr):
        self.pot = pot # the potential energy cut (but not shifted) at r_cut
        self.vir = vir # the virial
        self.lap = lap # the Laplacian
        self.ovr = ovr # a flag indicating overlap (i.e. pot too high to use)

    def __add__(self, other):
        pot = self.pot +  other.pot
        vir = self.vir +  other.vir
        lap = self.lap +  other.lap
        ovr = self.ovr or other.ovr
        return PotentialType(pot,vir,lap,ovr)

    def __sub__(self, other):
        pot = self.pot -  other.pot
        vir = self.vir -  other.vir
        lap = self.lap -  other.lap
        ovr = self.ovr or other.ovr # This is meaningless, but inconsequential
        return PotentialType(pot,vir,lap,ovr)

def introduction():
    """Prints out introductory statements at start of run."""

    print('Lennard-Jones potential')
    print('Cut (but not shifted)')
    print('Diameter, sigma = 1')
    print('Well depth, epsilon = 1')
    if fast:
        print('Fast NumPy potential routine')
    else:
        print('Slow Python potential routine')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def potential ( box, r_cut, r ):
    """Takes in box, cutoff range, and coordinate array, and calculates total potential etc.

    The results are returned as total, a PotentialType variable.
    """
    # Actual calculation performed by function potential_1

    n, d = r.shape
    assert d==3, 'Dimension error for r in potential'

    total = PotentialType ( pot=0.0, vir=0.0, lap=0.0, ovr=False )

    for i in range(n-1):
        partial = potential_1 ( r[i,:], box, r_cut, r[i+1:,:] )
        if partial.ovr:
            total.ovr = True
            break
        total = total + partial

    return total

def potential_1 ( ri, box, r_cut, r ):
    """Takes in coordinates of an atom and calculates its interactions.

    Values of box, cutoff range, and partner coordinate array are supplied.
    The results are returned as partial, a PotentialType variable.
    """

    import numpy as np

    # partial.pot is the nonbonded cut (not shifted) potential energy of atom ri with a set of other atoms
    # partial.vir is the corresponding virial of atom ri
    # partial.lap is the corresponding Laplacian of atom ri
    # partial.ovr is a flag indicating overlap (potential too high) to avoid overflow
    # If this is True, the values of partial.pot etc should not be used
    # In general, r will be a subset of the complete set of simulation coordinates
    # and none of its rows should be identical to ri

    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1

    nj, d = r.shape
    assert d==3, 'Dimension error for r in potential_1'
    assert ri.size==3, 'Dimension error for ri in potential_1'

    sr2_ovr      = 1.77 # Overlap threshold (pot > 100)
    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    if fast:
        rij = ri - r                    # Get all separation vectors from partners
        rij = rij - np.rint(rij)        # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2,axis=1)  # Squared separations
        in_range = rij_sq < r_cut_box_sq                    # Set flags for within cutoff
        rij_sq   = rij_sq * box_sq                          # Now in sigma=1 units
        sr2      = np.where ( in_range, 1.0 / rij_sq, 0.0 ) # (sigma/rij)**2, only if in range
        ovr      = sr2 > sr2_ovr                            # Set flags for any overlaps
        if np.any(ovr):
            partial = PotentialType ( pot=0.0, vir=0.0, lap=0.0, ovr=True )
            return partial
        sr6  = sr2 ** 3
        sr12 = sr6 ** 2
        pot  = sr12 - sr6                    # LJ pair potentials (cut but not shifted)
        vir  = pot + sr12                    # LJ pair virials
        lap  = ( 22.0*sr12 - 5.0*sr6 ) * sr2 # LJ pair Laplacians
        partial = PotentialType ( pot=np.sum(pot), vir=np.sum(vir), lap=np.sum(lap), ovr=False )

    else:
        partial = PotentialType ( pot=0.0, vir=0.0, lap=0.0, ovr=False )
        for rj in r:
            rij = ri - rj            # Separation vector
            rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
            rij_sq = np.sum(rij**2)  # Squared separation
            if rij_sq < r_cut_box_sq: # Check within cutoff
                rij_sq = rij_sq * box_sq # Now in sigma=1 units
                sr2    = 1.0 / rij_sq    # (sigma/rij)**2
                ovr    = sr2 > sr2_ovr   # Overlap if too close
                if ovr:
                    partial.ovr=True
                    return partial
                sr6  = sr2 ** 3
                sr12 = sr6 ** 2
                pot  = sr12 - sr6                    # LJ pair potential (cut but not shifted)
                vir  = pot + sr12                    # LJ pair virial
                lap  = ( 22.0*sr12 - 5.0*sr6 ) * sr2 # LJ pair Laplacian
                partial = partial + PotentialType ( pot=pot, vir=vir, lap=lap, ovr=ovr )

    # Multiply results by numerical factors
    partial.pot = partial.pot * 4.0        # 4*epsilon
    partial.vir = partial.vir * 24.0 / 3.0 # 24*epsilon and divide virial by 3
    partial.lap = partial.lap * 24.0 * 2.0 # 24*epsilon and factor 2 for ij and ji

    return partial

def force_sq ( box, r_cut, r ):
    """Calculates total squared force."""

    import numpy as np

    n, d = r.shape

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    f = np.zeros_like(r) # Initialize

    if fast:
        for i in range(n-1):
            rij = r[i,:]-r[i+1:,:]            # Separation vectors for j>i
            rij = rij - np.rint(rij)          # Periodic boundary conditions in box=1 units
            rij_sq    = np.sum(rij**2,axis=1) # Squared separations for j>1
            in_range  = rij_sq < r_cut_box_sq # Set flags for within cutoff
            rij_sq    = rij_sq * box_sq       # Now in sigma=1 units
            rij       = rij * box             # Now in sigma=1 units
            sr2       = np.where ( in_range, 1.0 / rij_sq, 0.0 ) # (sigma/rij)**2, only if in range
            sr6       = sr2 ** 3
            sr12      = sr6 ** 2
            fij       = (2*sr12 - sr6)*sr2
            fij       = rij * fij[:,np.newaxis] # LJ pair forces
            f[i,:]    = f[i,:] + np.sum(fij,axis=0)
            f[i+1:,:] = f[i+1:,:] - fij
    else:
        for i in range(n-1): # Outer loop over atoms
            for j in range(i+1,n): # Inner loop over atoms
                rij = r[i,:] - r[j,:]       # Separation vector
                rij = rij - np.rint ( rij ) # Periodic boundary conditions in box=1 units
                rij_sq = np.sum ( rij**2 )  # Squared separation

                if rij_sq < r_cut_box_sq: # Check within cutoff

                    rij_sq = rij_sq * box_sq # Now in sigma=1 units
                    rij = rij * box          # Now in sigma=1 units
                    sr2    = 1.0 / rij_sq
                    sr6    = sr2 ** 3
                    sr12   = sr6 ** 2
                    fij    = rij * (2.0*sr12 - sr6) *sr2 # LJ pair forces
                    f[i,:] = f[i,:] + fij
                    f[j,:] = f[j,:] - fij

    f = f * 24 # Numerical factor 24*epsilon

    return np.sum(f**2)
