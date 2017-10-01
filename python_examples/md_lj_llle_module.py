#!/usr/bin/env python3
# md_lj_llle_module.py

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

"""Force routine for MD simulation, LJ, Lees-Edwards, using neighbour lists."""

fast = True # Change this to replace NumPy force evaluation with slower Python
r_cut_sq = 2.0**(1.0/3.0)
r_cut    = r_cut_sq**(1.0/2.0) # Minimum and cutoff of WCA LJ potential

class PotentialType:
    """A composite variable for interactions."""

    def __init__(self, pot, vir, lap, pyx, ovr):
        self.pot = pot # the potential energy cut-and-shifted at r_cut
        self.vir = vir # the virial
        self.lap = lap # the Laplacian
        self.pyx = pyx # the off-diagonal virial part of the pressure tensor
        self.ovr = ovr # a flag indicating overlap (i.e. pot too high to use)

    def __add__(self, other):
        pot = self.pot +  other.pot
        vir = self.vir +  other.vir
        lap = self.lap +  other.lap
        pyx = self.pyx +  other.pyx
        ovr = self.ovr or other.ovr

        return PotentialType(pot,vir,lap,pyx,ovr)

def introduction():
    """Prints out introductory statements at start of run."""
    
    print('WCA shifted Lennard-Jones potential')
    print('Diameter, sigma = 1')
    print('Well depth, epsilon = 1')
    print('Lees-Edwards boundaries')
    if fast:
        print('Fast NumPy force routine')
    else:
        print('Slow Python force routine')
    print('Uses neighbour lists')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def force ( box, strain, r ):
    """Takes in box, strain, and coordinate array, and calculates forces and potentials etc."""

    import numpy as np
    from itertools import product
    import math
    
    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1
    # Lees-Edwards boundaries, in sliding brick arrangement
    # Flow/gradient/vorticity directions are x/y/z == 0/1/2
    # Uses neighbour lists

    n = r.shape[0]

    # Set up vectors to half the cells in neighbourhood of 3x3x3 cells in cubic lattice
    # The cells are chosen so that if (d0,d1,d2) appears, then (-d0,-d1,-d2) does not.
    # The last three cells are extra ones, to cope with the sheared system
    d = np.array ( [ [ 0, 0,  0], [ 1, 0, 0], [ 1, 0, 1], [-1, 0, 1], [ 0, 0, 1], # 5 cells with d1=0
                     [ 1, 1, -1], [ 1, 1, 0], [ 1, 1, 1],    # 3 cells with d0= 1, d1=1
                     [ 0, 1, -1], [ 0, 1, 0], [ 0, 1, 1],    # 3 cells with d0= 0, d1=1
                     [-1, 1, -1], [-1, 1, 0], [-1, 1, 1],    # 3 cells with d0=-1, d1=1
                     [-2, 1, -1], [-2, 1, 0], [-2, 1, 1] ] ) # 3 cells with d0=-2, d1=1

    r[:,0] = r[:,0] - np.rint(r[:,1])*strain # Extra correction in box=1 units
    r      = r      - np.rint(r)             # Ensure all atoms in periodic box
    
    sr2_ovr      = 1.77 # Overlap threshold (pot > 100)
    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    # Initialize
    f     = np.zeros_like(r)
    total = PotentialType ( pot=0.0, vir=0.0, pyx=0.0, lap=0.0, ovr=False )

    # Calculate cell index triplets
    sc = math.floor(box/r_cut)                          # Number of cells along box edge
    c  = np.floor((r+0.5)*sc).astype(np.int_)           # N*3 array of cell indices for all atoms
    assert np.all(c>=0) and np.all(c<sc), 'Index error' # Simplistic "guard" against roundoff

    shift = math.floor(strain*sc) # Strain measured in cell lengths

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

            # Set up correct neighbour cell indices
            if ci[1]==sc-1:                # i-cell is in the top layer
                dd       = d.copy()        # Standard list copied, including extra 3 cells
                dd[5:,0] = d[5:,0] - shift # All those looking up need adjustment in the x direction
            else:                          # i-cell is not in top layer
                dd = d[:-3,:].copy()       # Last three extra cells are not needed; shift is not needed
            
            for dj in dd:                                             # Loop over neighbouring j-cells
                cj  = ci + dj                                         # Compute neighbour j-cell triple-indices
                cj1 = np.ravel_multi_index(cj,(sc,sc,sc),mode='wrap') # Convert j-cell to single-index
                rcj = rc[cj1]                                         # Get atoms in j-cell as an array
                if rcj.size==0:                                       # Handle empty cell
                    continue

                rij        = rci[:,np.newaxis,:]-rcj[np.newaxis,:,:] # Separation vectors for all i and j
                rij[:,:,0] = rij[:,:,0] - np.rint(rij[:,:,1])*strain # Extra correction in box=1 units
                rij        = rij - np.rint(rij)                      # PBCs in box=1 units
                rij_sq     = np.sum(rij**2,axis=2)                   # Squared separations
                in_range   = rij_sq < r_cut_box_sq                   # Set flags for within cutoff

                if ci1==cj1:
                    np.fill_diagonal(in_range,False) # Eliminate i==j when i-cell==j-cell
                    np.fill_diagonal(rij_sq,1.0)     # Avoid divide-by-zero below

                rij_sq = rij_sq * box_sq                         # Now in sigma=1 units
                rij    = rij * box                               # Now in sigma=1 units
                sr2    = np.where ( in_range, 1.0/rij_sq, 0.0 )  # (sigma/rij)**2, only if in range
                ovr    = sr2 > sr2_ovr                           # Overlap if too close
                sr6    = sr2 ** 3
                sr12   = sr6 ** 2
                pot    = sr12 - sr6                              # LJ potential (cut but not shifted)
                vir    = pot + sr12                              # LJ virial
                pot    = np.where ( in_range, pot+0.25, 0.0 )    # WCA LJ pair potential (cut-and-shifted)
                lap    = ( 22.0*sr12 - 5.0*sr6 ) * sr2           # LJ Laplacian
                fij    = vir * sr2                               # LJ scalar part of forces
                fij    = rij * fij[:,:,np.newaxis]               # LJ pair forces
                pyx    = rij[:,:,1]*fij[:,:,0]                   # Off-diagonal element of pressure tensor

                if ci1==cj1: # Correct for double-counting ij and ji when i-cell==j-cell
                    fij = fij / 2
                    total = total + PotentialType ( pot=np.sum(pot)/2, vir=np.sum(vir)/2, 
                                                    pyx=np.sum(pyx)/2, lap=np.sum(lap)/2, ovr=np.any(ovr) )
                else:
                    total = total + PotentialType ( pot=np.sum(pot), vir=np.sum(vir), 
                                                    pyx=np.sum(pyx), lap=np.sum(lap), ovr=np.any(ovr) )

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

            # Set up correct neighbour cell indices
            if ci[1]==sc-1: # i-cell is in the top layer
                dd       = d                # Standard list copied, including extra 3 cells
                dd[5:,0] = dd[5:,0] - shift # All those looking up need adjustment in the x direction
            else:
                dd = d[:-3,:] # Last three extra cells are not needed; shift is not needed

            for dj in dd: # Loop over neighbouring j-cells
                cj  = ci + dj                                         # Compute neighbour j-cell triple-indices
                cj1 = np.ravel_multi_index(cj,(sc,sc,sc),mode='wrap') # Convert to single-index
                kcj = kc[cj1]                                         # Get indices of atoms in j-cell as an array

                for i, ki in enumerate(kci):     # Loop over individual atoms in i-cell
                    j0 = i+1 if cj1==ci1 else 0  # Only look upwards if i-cell==j-cell
                    if j0 >= kcj.size:           # Handles (redundantly) empty j-cell and the case 
                        continue                 # where j-cell==i-cell and i is last atom

                    for kj in kcj[j0:]:                          # Loop over individual atoms in j-cell
                        rij    = r[ki,:]-r[kj,:]                 # Separation vector
                        rij[0] = rij[0] - np.rint(rij[1])*strain # Extra correction in box=1 units
                        rij    = rij - np.rint(rij)              # Periodic boundary conditions in box=1 units
                        rij_sq = np.sum(rij**2)                  # Squared separation

                        if rij_sq < r_cut_box_sq: # Check within cutoff
                            rij_sq  = rij_sq * box_sq # Now in sigma=1 units
                            rij     = rij * box       # Now in sigma=1 units
                            sr2     = 1.0 / rij_sq    # (sigma/rij)**2
                            ovr     = sr2 > sr2_ovr   # Overlap if too close
                            sr6     = sr2 ** 3
                            sr12    = sr6 ** 2
                            pot     = sr12 - sr6                    # LJ potential (cut but not shifted)
                            vir     = pot + sr12                    # LJ virial
                            pot     = pot + 0.25                    # WCA LJ potential (cut-and-shifted)
                            lap     = ( 22.0*sr12 - 5.0*sr6 ) * sr2 # LJ Laplacian
                            fij     = rij * vir * sr2               # LJ forces
                            pyx     = rij[1]*fij[0]                 # Off-diagonal element of pressure tensor
                            total   = total + PotentialType ( pot=pot, vir=vir, pyx=pyx, lap=lap, ovr=ovr )
                            f[ki,:] = f[ki,:] + fij
                            f[kj,:] = f[kj,:] - fij

    # Multiply results by numerical factors
    f         = f         * 24.0       # 24*epsilon
    total.pot = total.pot * 4.0        # 4*epsilon
    total.vir = total.vir * 24.0 / 3.0 # 24*epsilon and divide virial by 3
    total.pyx = total.pyx * 24.0       # 24*epsilon
    total.lap = total.lap * 24.0 * 2.0 # 24*epsilon and factor 2 for ij and ji
    
    return total, f

