#!/usr/bin/env python3
"""Force routine for MD simulation, Lennard-Jones atoms."""

class PotentialType:
    """A composite variable for interactions."""

    def __init__(self, cut, pot, vir, lap, ovr):
        self.cut = cut # the potential energy cut (but not shifted) at r_cut
        self.pot = pot # the potential energy cut-and-shifted at r_cut
        self.vir = vir # the virial
        self.lap = lap # the Laplacian
        self.ovr = ovr # a flag indicating overlap (i.e. pot too high to use)

    def __add__(self, other):
        cut = self.cut +  other.cut
        pot = self.pot +  other.pot
        vir = self.vir +  other.vir
        lap = self.lap +  other.lap
        ovr = self.ovr or other.ovr

        return PotentialType(cut,pot,vir,lap,ovr)

def introduction():
    """Prints out introductory statements at start of run."""
    
    print('Lennard-Jones potential')
    print('Cut-and-shifted version for dynamics')
    print('Cut (but not shifted) version also calculated')
    print('Diameter, sigma = 1')
    print('Well depth, epsilon = 1')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def force ( box, r_cut, r ):
    """Takes in box, cutoff range, and coordinate array, and calculates forces and potentials etc."""

    import numpy as np
    
    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1

    n, d = r.shape
    assert d==3, 'Dimension error in force'
    
    sr2_ovr      = 1.77 # Overlap threshold (pot > 100)
    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    # Calculate potential at cutoff
    sr2     = 1.0 / r_cut**2 # in sigma=1 units
    sr6     = sr2 ** 3
    sr12    = sr6 **2
    pot_cut = sr12 - sr6 # Without numerical factor 4

    # Initialize
    f = np.zeros_like(r)
    total = PotentialType ( cut=0.0, pot=0.0, vir=0.0, lap=0.0, ovr=False )

    for i in range(n-1):
        for j in range(i+1,n):
            rij = r[i,:]-r[j,:]      # Separation vector
            rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
            rij_sq = np.sum(rij**2)  # Squared separation

            if rij_sq < r_cut_box_sq: # Check within cutoff
                rij_sq = rij_sq * box_sq # Now in sigma=1 units
                rij    = rij * box       # Now in sigma=1 units
                sr2    = 1.0 / rij_sq    # (sigma/rij)**2
                ovr    = sr2 > sr2_ovr   # Overlap if too close

                sr6  = sr2 ** 3
                sr12 = sr6 ** 2
                cut  = sr12 - sr6                    # LJ pair potential (cut but not shifted)
                vir  = cut + sr12                    # LJ pair virial
                pot  = cut - pot_cut                 # LJ pair potential (cut-and-shifted)
                lap  = ( 22.0*sr12 - 5.0*sr6 ) * sr2 # LJ pair Laplacian
                fij  = rij * vir * sr2               # LJ pair forces

                total  = total + PotentialType ( cut=cut, pot=pot, vir=vir, lap=lap, ovr=ovr )
                f[i,:] = f[i,:] + fij
                f[j,:] = f[j,:] - fij

    # Multiply results by numerical factors
    f         = f         * 24.0       # 24*epsilon
    total.cut = total.cut * 4.0        # 4*epsilon
    total.pot = total.pot * 4.0        # 4*epsilon
    total.vir = total.vir * 24.0 / 3.0 # 24*epsilon and divide virial by 3
    total.lap = total.lap * 24.0 * 2.0 # 24*epsilon and factor 2 for ij and ji
    
    return total, f

def force_faster ( box, r_cut, r ):
    """Takes in box, cutoff range, and coordinate array, and calculates forces and potentials etc."""

    import numpy as np
    
    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1

    n, d = r.shape
    assert d==3, 'Dimension error in force_faster'
    
    sr2_ovr      = 1.77 # Overlap threshold (pot > 100)
    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    # Calculate potential at cutoff
    sr2     = 1.0 / r_cut**2 # in sigma=1 units
    sr6     = sr2 ** 3
    sr12    = sr6 **2
    pot_cut = sr12 - sr6 # Without numerical factor 4

    # Initialize
    f = np.zeros_like(r)
    total = PotentialType ( cut=0.0, pot=0.0, vir=0.0, lap=0.0, ovr=False )

    for i in range(n-1):
        rij = r[i,:]-r[i+1:,:]           # Separation vectors for j>i
        rij = rij - np.rint(rij)         # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2,axis=1)   # Squared separations for j>1
        in_range = rij_sq < r_cut_box_sq # Set flags for within cutoff

        rij_sq = rij_sq * box_sq # Now in sigma=1 units
        rij    = rij * box       # Now in sigma=1 units
        sr2    = np.where ( in_range, 1.0 / rij_sq, 0.0 ) # (sigma/rij)**2, only if in range
        ovr    = sr2 > sr2_ovr   # Overlap if too close

        sr6  = sr2 ** 3
        sr12 = sr6 ** 2
        cut  = sr12 - sr6                                  # LJ pair potential (cut but not shifted)
        vir  = cut + sr12                                  # LJ pair virial
        pot  = np.where ( in_range, cut - pot_cut, 0.0 )   # LJ pair potential (cut-and-shifted)
        lap  = ( 22.0*sr12 - 5.0*sr6 ) * sr2               # LJ pair Laplacian
        fij  = rij * vir[:,np.newaxis] * sr2[:,np.newaxis] # LJ pair forces

        total     = total + PotentialType ( cut=np.sum(cut), pot=np.sum(pot), vir=np.sum(vir), lap=np.sum(lap), ovr=np.any(ovr) )
        f[i,:]    = f[i,:] + np.sum(fij,axis=0)
        f[i+1:,:] = f[i+1:,:] - fij

    # Multiply results by numerical factors
    f         = f         * 24.0       # 24*epsilon
    total.cut = total.cut * 4.0        # 4*epsilon
    total.pot = total.pot * 4.0        # 4*epsilon
    total.vir = total.vir * 24.0 / 3.0 # 24*epsilon and divide virial by 3
    total.lap = total.lap * 24.0 * 2.0 # 24*epsilon and factor 2 for ij and ji

    return total, f

def hessian ( box, r_cut, r, f ):
    """Calculates Hessian function (for 1/N correction to config temp)."""

    import numpy as np

    # This routine is only needed in a constant-energy ensemble
    # It is assumed that positions are in units where box = 1
    # but the result is given in units where sigma = 1 and epsilon = 1
    # It is assumed that forces have already been calculated in array f

    n, d = r.shape
    assert d==3, 'Dimension error in hessian'
    assert np.all ( r.shape==f.shape ), 'Dimension mismatch in hessian'

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    hes = 0.0

    for i in range(n-1):
        for j in range(i+1,n):
            rij = r[i,:] - r[j,:]       # Separation vector
            rij = rij - np.rint ( rij ) # Periodic boundary conditions in box=1 units
            rij_sq = np.sum ( rij**2 )  # Squared separation

            if rij_sq < r_cut_box_sq:
                rij_sq = rij_sq * box_sq # Now in sigma=1 units
                rij = rij * box    # Now in sigma=1 units
                fij = f[i,:] - f[j,:] # Difference in forces

                ff   = np.dot(fij,fij)
                rf   = np.dot(rij,fij)
                sr2  = 1.0 / rij_sq
                sr6  = sr2 ** 3
                sr8  = sr6 * sr2
                sr10 = sr8 * sr2
                v1   = 24.0 * ( 1.0 - 2.0 * sr6 ) * sr8
                v2   = 96.0 * ( 7.0 * sr6 - 2.0 ) * sr10
                hes  = hes + v1 * ff + v2 * rf**2

    return hes

def hessian_faster ( box, r_cut, r, f ):
    """Calculates Hessian function (for 1/N correction to config temp)."""

    import numpy as np

    # This routine is only needed in a constant-energy ensemble
    # It is assumed that positions are in units where box = 1
    # but the result is given in units where sigma = 1 and epsilon = 1
    # It is assumed that forces have already been calculated in array f

    n, d = r.shape
    assert d==3, 'Dimension error in hessian_faster'
    assert np.all ( r.shape==f.shape ), 'Dimension mismatch in hessian_faster'

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    hes = 0.0

    for i in range(n-1):
        rij = r[i,:] - r[i+1:,:]         # Separation vectors
        rij = rij - np.rint ( rij )      # Periodic boundary conditions in box=1 units
        rij_sq = np.sum ( rij**2 )       # Squared separation
        in_range = rij_sq < r_cut_box_sq # Set flags for within cutoff

        rij_sq = rij_sq * box_sq # Now in sigma=1 units
        rij = rij * box          # Now in sigma=1 units
        fij = f[i,:] - f[i+1:,:] # Differences in forces

        ff   = np.sum(fij*fij,axis=1)
        rf   = np.sum(rij*fij,axis=1)
        sr2  = np.where ( in_range, 1.0 / rij_sq, 0.0 ) # Only where in range
        sr6  = sr2 ** 3
        sr8  = sr6 * sr2
        sr10 = sr8 * sr2
        v1   = 24.0 * ( 1.0 - 2.0 * sr6 ) * sr8
        v2   = 96.0 * ( 7.0 * sr6 - 2.0 ) * sr10
        hes  = hes + np.sum(v1 * ff) + np.sum(v2 * rf**2)

    return hes

