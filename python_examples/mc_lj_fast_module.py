#!/usr/bin/env python3
"""Energy and move routines for MC simulation, LJ potential.

Fast version, using NumPy.
"""

lt, gt = -1, 1 # Options for j_range

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
    print('Fast version built around NumPy routines')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def potential ( box, r_cut, r ):
    """Takes in box, cutoff range, and coordinate array, and calculates total potential etc.

    The results are returned as total, a PotentialType variable.
    """

    n, d = r.shape
    assert d==3, 'Dimension error for r in potential'

    total = PotentialType ( pot=0.0, vir=0.0, lap=0.0, ovr=False ) # Initialize

    for i in range(n-1):
        partial = potential_1 ( r[i,:], i, box, r_cut, r, gt )
        if partial.ovr:
            total.ovr = True
            break
        total = total + partial

    return total

def potential_1 ( ri, i, box, r_cut, r, j_range=None ):
    """Takes in coordinates and index of an atom and calculates its interactions.

    Values of box, cutoff range, and coordinate array are supplied,
    as is the (optional) range of j-neighbours to be considered.
    """

    import numpy as np
    import sys

    # partial.pot is the nonbonded cut (not shifted) potential energy of atom ri with a set of other atoms
    # partial.vir is the corresponding virial of atom ri
    # partial.lap is the corresponding Laplacian of atom ri
    # partial.ovr is a flag indicating overlap (potential too high) to avoid overflow
    # If this is True, the values of partial.pot etc should not be used
    # The coordinates in ri are not necessarily identical with those in r[i,:]
    # The argument j_range may restrict partner indices to j>i, or j<i
    # Usually i takes a value in range(n) indicating the atom to be skipped
    # but the routine should also handle the case j_range=None, i=n, when no atoms are skipped

    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1

    n, d = r.shape
    assert d==3, 'Dimension error for r in potential_1'
    assert ri.size==3, 'Dimension error for ri in potential_1'

    sr2_ovr      = 1.77 # Overlap threshold (pot > 100)
    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    # Get all separation vectors from partners excluding i
    if j_range is None:
        if i in range(n):
            rij = ri-np.delete(r,i,0)
        else:
            rij = ri-r
    else:
        assert j_range==gt or j_range==lt, 'j-range error in potential_1'
        rij = ri-r[:i,:] if j_range==lt else ri-r[i+1:,:]

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

    for i in range(n-1):
        rij = r[i,:]-r[i+1:,:]           # Separation vectors for j>i
        rij = rij - np.rint(rij)         # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2,axis=1)   # Squared separations for j>1
        in_range = rij_sq < r_cut_box_sq # Set flags for within cutoff

        rij_sq = rij_sq * box_sq # Now in sigma=1 units
        rij    = rij * box       # Now in sigma=1 units
        sr2    = np.where ( in_range, 1.0 / rij_sq, 0.0 ) # (sigma/rij)**2, only if in range

        sr6  = sr2 ** 3
        sr12 = sr6 ** 2
        fac  = (2*sr12 - sr6)*sr2
        fij  = rij * fac[:,np.newaxis] # LJ pair forces

        f[i,:]    = f[i,:] + np.sum(fij,axis=0)
        f[i+1:,:] = f[i+1:,:] - fij

    f = f * 24 # Numerical factor 24*epsilon

    return np.sum(f**2)
