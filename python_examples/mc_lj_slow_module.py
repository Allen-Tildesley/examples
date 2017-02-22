#!/usr/bin/env python3
"""Energy and move routines for MC simulation, LJ potential.

Slow version using Python loop.
"""

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
    print('Slow version built around Python loops')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def potential ( box, r_cut, r ):
    """Takes in box, cutoff range, and coordinate array, and calculates total potential etc.

    The results are returned as total, a PotentialType variable."""

    # total.pot is the nonbonded cut (not shifted) potential energy for whole system
    # total.vir is the corresponding virial for whole system
    # total.lap is the corresponding Laplacian for whole system
    # total.ovr is a flag indicating overlap (potential too high) to avoid overflow
    # If this flag is True, the values of total.pot etc should not be used
    # Actual calculation is performed by function potential_1

    n, d = r.shape
    assert d==3, 'Dimension error for r in potential'

    total = PotentialType ( pot=0.0, vir=0.0, lap=0.0, ovr=False ) # Initialize

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

    # Initialize
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
