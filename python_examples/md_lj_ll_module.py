#!/usr/bin/env python3
# md_lj_ll_module.py

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

"""Force routine for MD simulation, LJ atoms, using neighbour lists."""

fast = True # Change this to replace NumPy force evaluation with slower Python

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
    if fast:
        print('Fast NumPy force routine')
    else:
        print('Slow Python force routine')
    print('Uses neighbour lists')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def force ( box, r_cut, r ):
    """Takes in box, cutoff range, and coordinate array, and calculates forces and potentials etc."""

    import numpy as np
    from itertools import product
    import math
    
    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1
    # Uses neighbour lists

    n = r.shape[0]

    # Set up vectors to half the cells in neighbourhood of 3x3x3 cells in cubic lattice
    # The cells are chosen so that if (d0,d1,d2) appears, then (-d0,-d1,-d2) does not.
    d = np.array ( [ [ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [-1, 1, 0],
                     [ 0, 1, 0], [ 0, 0, 1], [-1, 0, 1], [ 1, 0, 1], [-1,-1, 1],
                     [ 0,-1, 1], [ 1,-1, 1], [-1, 1, 1], [ 0, 1, 1], [ 1, 1, 1] ] )

    r = r - np.rint(r) # Ensure all atoms in periodic box
    
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
    f     = np.zeros_like(r)
    total = PotentialType ( cut=0.0, pot=0.0, vir=0.0, lap=0.0, ovr=False )

    # Calculate cell index triplets
    sc = math.floor(box/r_cut)                          # Number of cells along box edge
    c  = np.floor((r+0.5)*sc).astype(np.int_)           # N*3 array of cell indices for all atoms
    assert np.all(c>=0) and np.all(c<sc), 'Index error' # Simplistic "guard" against roundoff

    if fast:
        
        # Build list of arrays, each array holding positions of atoms in a cell
        # At the same time, define a matching set of force arrays in each cell
        # i and j number the atoms in each cell; we do not refer explicitly to indices in r
        rc, fc = [], []                        # Initially empty lists of positions and forces
        for ci in product(range(sc),repeat=3): # Triple loop over cells
            mask = np.all(c==ci,axis=1)        # Mask identifies atoms in this cell
            rc.append(r[mask,:])               # Copy atom coordinates into array, add to list
            fc.append(np.zeros_like(rc[-1]))   # Zero corresponding forces, add to list

        for ci1, rci in enumerate(rc):            # Loop over i-cells, getting all atoms in each i-cell as an array
            ci = np.unravel_index(ci1,(sc,sc,sc)) # Get i-cell triple-indices
            if rci.size==0:                       # Handle empty cell
                continue

            for dj in d:                                              # Loop over neighbouring j-cells
                cj  = ci + dj                                         # Compute neighbour j-cell triple-indices
                cj1 = np.ravel_multi_index(cj,(sc,sc,sc),mode='wrap') # Convert j-cell to single-index
                rcj = rc[cj1]                                         # Get atoms in j-cell as an array
                if rcj.size==0:                                       # Handle empty cell
                    continue

                rij      = rci[:,np.newaxis,:]-rcj[np.newaxis,:,:] # Separation vectors for all i and j
                rij      = rij - np.rint(rij)                      # PBCs in box=1 units
                rij_sq   = np.sum(rij**2,axis=2)                   # Squared separations
                in_range = rij_sq < r_cut_box_sq                   # Set flags for within cutoff

                if ci1==cj1:
                    np.fill_diagonal(in_range,False) # Eliminate i==j when i-cell==j-cell
                    np.fill_diagonal(rij_sq,1.0)     # Avoid divide-by-zero below

                rij_sq = rij_sq * box_sq                         # Now in sigma=1 units
                rij    = rij * box                               # Now in sigma=1 units
                sr2    = np.where ( in_range, 1.0/rij_sq, 0.0 )  # (sigma/rij)**2, only if in range
                ovr    = sr2 > sr2_ovr                           # Overlap if too close
                sr6    = sr2 ** 3
                sr12   = sr6 ** 2
                cut    = sr12 - sr6                              # LJ potential (cut but not shifted)
                vir    = cut + sr12                              # LJ virial
                pot    = np.where ( in_range, cut-pot_cut, 0.0 ) # LJ potential (cut-and-shifted)
                lap    = ( 22.0*sr12 - 5.0*sr6 ) * sr2           # LJ Laplacian
                fij    = vir * sr2                               # LJ scalar part of forces
                fij    = rij * fij[:,:,np.newaxis]               # LJ pair forces

                if ci1==cj1: # Correct for double-counting ij and ji when i-cell==j-cell
                    fij = fij / 2
                    total = total + PotentialType ( cut=np.sum(cut)/2, pot=np.sum(pot)/2,
                                                    vir=np.sum(vir)/2, lap=np.sum(lap)/2, ovr=np.any(ovr) )
                else:
                    total = total + PotentialType ( cut=np.sum(cut), pot=np.sum(pot),
                                                    vir=np.sum(vir), lap=np.sum(lap), ovr=np.any(ovr) )

                fc[ci1][:,:] = fc[ci1][:,:] + np.sum(fij,axis=1) # Aggregate force on atoms in i-cell
                fc[cj1][:,:] = fc[cj1][:,:] - np.sum(fij,axis=0) # Aggregate force on atoms in j-cell

        # Copy forces from list of cell arrays to main force array
        for ci in product(range(sc),repeat=3):                          # Triple loop over cells
            mask      = np.all(c==ci,axis=1)                            # Mask identifies atoms in this cell
            ci1       = np.ravel_multi_index(ci,(sc,sc,sc),mode='wrap') # Single-index
            f[mask,:] = fc[ci1]                                         # Copy atom forces from correct cell

    else:
                
        # Build list of arrays, each array holding indices of atoms in a cell
        # ki and kj are atom indices in the r array; i and j number the atoms in each cell
        k_array = np.arange(n)                 # Atom indices 0..N-1
        kc = []                                # Initially empty list of indices
        for ci in product(range(sc),repeat=3): # Triple loop over cells
            mask = np.all(c==ci,axis=1)        # Mask identifies atoms in this cell
            kc.append(k_array[mask])           # Copy atom indices into array, add to list

        for ci1, kci in enumerate(kc): # Loop over i-cells, getting atom indices as an array
            ci = np.unravel_index(ci1,(sc,sc,sc)) # Get i-cell triple-indices

            for dj in d: # Loop over neighbouring j-cells
                cj  = ci + dj                                         # Compute neighbour j-cell triple-indices
                cj1 = np.ravel_multi_index(cj,(sc,sc,sc),mode='wrap') # Convert to single-index
                kcj = kc[cj1]                                         # Get indices of atoms in j-cell as an array

                for i, ki in enumerate(kci):     # Loop over individual atoms in i-cell
                    j0 = i+1 if cj1==ci1 else 0  # Only look upwards if i-cell==j-cell
                    if j0 >= kcj.size:           # Handles (redundantly) empty j-cell and the case 
                        continue                 # where j-cell==i-cell and i is last atom

                    for kj in kcj[j0:]:          # Loop over individual atoms in j-cell
                        rij = r[ki,:]-r[kj,:]    # Separation vector
                        rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
                        rij_sq = np.sum(rij**2)  # Squared separation

                        if rij_sq < r_cut_box_sq: # Check within cutoff
                            rij_sq  = rij_sq * box_sq # Now in sigma=1 units
                            rij     = rij * box       # Now in sigma=1 units
                            sr2     = 1.0 / rij_sq    # (sigma/rij)**2
                            ovr     = sr2 > sr2_ovr   # Overlap if too close
                            sr6     = sr2 ** 3
                            sr12    = sr6 ** 2
                            cut     = sr12 - sr6                    # LJ potential (cut but not shifted)
                            vir     = cut + sr12                    # LJ virial
                            pot     = cut - pot_cut                 # LJ potential (cut-and-shifted)
                            lap     = ( 22.0*sr12 - 5.0*sr6 ) * sr2 # LJ Laplacian
                            fij     = rij * vir * sr2               # LJ forces
                            total   = total + PotentialType ( cut=cut, pot=pot, vir=vir, lap=lap, ovr=ovr )
                            f[ki,:] = f[ki,:] + fij
                            f[kj,:] = f[kj,:] - fij

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
    from itertools import product
    import math

    # This routine is only needed in a constant-energy ensemble
    # It is assumed that positions are in units where box = 1
    # but the result is given in units where sigma = 1 and epsilon = 1
    # It is assumed that forces have already been calculated in array f
    # Uses neighbour lists

    n = r.shape[0]
    assert np.all ( r.shape==f.shape ), 'Dimension mismatch in hessian'

    # Set up vectors to half the cells in neighbourhood of 3x3x3 cells in cubic lattice
    # The cells are chosen so that if (d1,d2,d3) appears, then (-d1,-d2,-d3) does not.
    d = np.array ( [ [ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [-1, 1, 0],
                     [ 0, 1, 0], [ 0, 0, 1], [-1, 0, 1], [ 1, 0, 1], [-1,-1, 1],
                     [ 0,-1, 1], [ 1,-1, 1], [-1, 1, 1], [ 0, 1, 1], [ 1, 1, 1] ] )

    r = r - np.rint(r) # Ensure all atoms in periodic box

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    hes = 0.0

    # Calculate cell index triplets
    sc = math.floor(box/r_cut)                          # Number of cells along box edge
    c  = np.floor((r+0.5)*sc).astype(np.int_)           # N*3 array of cell indices for all atoms
    assert np.all(c>=0) and np.all(c<sc), 'Index error' # Simplistic "guard" against roundoff

    if fast:
        
        # Build list of arrays, each array holding positions of atoms in a cell
        # At the same time, build a matching set of force arrays in each cell
        # i and j number the atoms in each cell; we do not refer explicitly to indices in r
        rc, fc = [], []                        # Initially empty lists of positions and forces
        for ci in product(range(sc),repeat=3): # Triple loop over cells
            mask = np.all(c==ci,axis=1)        # Mask identifies atoms in this cell
            rc.append(r[mask,:])               # Copy atom coordinates into array, add to list
            fc.append(f[mask,:])               # Copy corresponding forces, add to list
        
        for ci1, rci in enumerate(rc):            # Loop over i-cells, getting all atoms in each i-cell as an array
            ci = np.unravel_index(ci1,(sc,sc,sc)) # Get i-cell triple-indices
            fci = fc[ci1]                         # Get i-cell atom forces
            if rci.size==0:                       # Handle empty cell
                continue

            for dj in d:                                              # Loop over neighbouring j-cells
                cj = ci + dj                                          # Compute neighbour j-cell triple-indices
                cj1 = np.ravel_multi_index(cj,(sc,sc,sc),mode='wrap') # Convert j-cell to single-index
                rcj = rc[cj1]                                         # Get atoms in j-cell as an array
                fcj = fc[cj1]                                         # Get j-cell atom forces
                if rcj.size==0:                                       # Handle empty cell
                    continue

                rij      = rci[:,np.newaxis,:]-rcj[np.newaxis,:,:] # Separation vectors for all i and j
                rij      = rij - np.rint(rij)                      # PBCs in box=1 units
                rij_sq   = np.sum(rij**2,axis=2)                   # Squared separations
                in_range = rij_sq < r_cut_box_sq                   # Set flags for within cutoff

                if ci1==cj1:
                    np.fill_diagonal(in_range,False) # Eliminate i=j when i-cell is j-cell
                    np.fill_diagonal(rij_sq,1.0)     # Avoid divide-by-zero below

                rij_sq   = rij_sq * box_sq                         # Now in sigma=1 units
                rij      = rij * box                               # Now in sigma=1 units
                fij      = fci[:,np.newaxis,:]-fcj[np.newaxis,:,:] # Differences in forces for all i and j

                ff   = np.sum(fij*fij,axis=2)
                rf   = np.sum(rij*fij,axis=2)
                sr2  = np.where ( in_range, 1.0 / rij_sq, 0.0 ) # Only where in range
                sr6  = sr2 ** 3
                sr8  = sr6 * sr2
                sr10 = sr8 * sr2
                v1   = 24.0 * ( 1.0 - 2.0 * sr6 ) * sr8
                v2   = 96.0 * ( 7.0 * sr6 - 2.0 ) * sr10
                if ci1==cj1: # Correct for double-counting ij and ji
                    hes  = hes + np.sum(v1 * ff)/2 + np.sum(v2 * rf**2)/2
                else:
                    hes  = hes + np.sum(v1 * ff) + np.sum(v2 * rf**2)

    else:
                
        # Build list of arrays, each array holding indices of atoms in a cell
        # ki and kj are atom indices in the r array; i and j number the atoms in each cell
        k_array = np.arange(n)                 # Atom indices 0..N-1
        kc = []                                # Initially empty list of indices covering each cell
        for ci in product(range(sc),repeat=3): # Triple loop over cells
            mask = np.all(c==ci,axis=1)        # Mask identifies atoms in this cell
            kc.append(k_array[mask])           # Copy atom indices into array, add to list

        for ci1, kci in enumerate(kc): # Loop over i-cells, getting atom indices as an array
            ci = np.unravel_index(ci1,(sc,sc,sc)) # Get i-cell triple-indices

            for dj in d: # Loop over neighbouring j-cells
                cj  = ci + dj                                         # Compute neighbour j-cell triple-indices
                cj1 = np.ravel_multi_index(cj,(sc,sc,sc),mode='wrap') # Convert to single-index
                kcj = kc[cj1]                                         # Get indices of atoms in j-cell as an array

                for i, ki in enumerate(kci):    # Loop over individual atoms in i-cell
                    j0 = i+1 if cj1==ci1 else 0 # Only look upwards if i-cell==j-cell
                    if j0 >= kcj.size:          # Handles (redundantly) empty j-cell and the case 
                        continue                # where j-cell==i-cell and i is last atom

                    for kj in kcj[j0:]:          # Loop over individual atoms in j-cell
                        rij = r[ki,:]-r[kj,:]    # Separation vector
                        rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
                        rij_sq = np.sum(rij**2)  # Squared separation

                        if rij_sq < r_cut_box_sq:
                            rij_sq = rij_sq * box_sq # Now in sigma=1 units
                            rij = rij * box          # Now in sigma=1 units
                            fij = f[ki,:] - f[kj,:]  # Difference in forces
                            
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
