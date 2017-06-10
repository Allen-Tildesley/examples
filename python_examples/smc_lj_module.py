#!/usr/bin/env python3
# smc_lj_module.py

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

"""Energy, force, and move routines for SMC, LJ potential."""

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

    def __sub__(self, other):
        cut = self.cut -  other.cut
        pot = self.pot -  other.pot
        vir = self.vir -  other.vir
        lap = self.lap -  other.lap
        ovr = self.ovr or other.ovr # This is meaningless, but inconsequential
        return PotentialType(cut,pot,vir,lap,ovr)

def introduction():
    """Prints out introductory statements at start of run."""
    
    print('Lennard-Jones potential')
    print('Cut-and-shifted version for SMC dynamics')
    print('Cut (but not shifted) version also calculated')
    print('Diameter, sigma = 1')
    print('Well depth, epsilon = 1')
    if fast:
        print('Fast NumPy force routine')
    else:
        print('Slow Python force routine')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def force ( box, r_cut, r ):
    """Takes in box, cutoff range, and coordinate array, and calculates forces and potentials etc."""

    import numpy as np

    # Actual calculation is performed by function force_1

    n, d = r.shape
    assert d==3, 'Dimension error for r in force'

    total = PotentialType ( pot=0.0, cut=0.0, vir=0.0, lap=0.0, ovr=False )
    f = np.zeros_like(r)

    for i in range(n-1):
        partial, f_partial = force_1 ( r[i,:], box, r_cut, r[i+1:,:] )
        if partial.ovr:
            total.ovr = True
            break
        total = total + partial
        f[i]      = f[i]      + np.sum(f_partial,axis=0)
        f[i+1:,:] = f[i+1:,:] - f_partial

    return total, f
  
def force_1 ( ri, box, r_cut, r ):
    """Takes in coordinates of an atom and calculates its interactions.

    Values of box, cutoff range, and partner coordinate array are supplied.
    The results are returned as partial, a PotentialType variable, and the forces f_partial.
    """

    import numpy as np

    # partial.pot is the cut-and-shifted potential energy of atom i with a set of other atoms
    # partial.cut is the cut (but not shifted) version of the above
    # partial.vir is the corresponding virial of atom i
    # partial.lap is the corresponding Laplacian of atom i
    # partial.ovr is a flag indicating overlap (potential too high) to avoid overflow
    # If this is True, the values of partial.pot etc should not be used
    # f_partial contains the force on ri due to each atom in r
    # It is assumed that the calling routine knows what to do with these

    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1
    # Note that we use a shifted LJ potential here

    nj, d = r.shape
    assert d==3, 'Dimension error for r in force_1'
    assert ri.size==3, 'Dimension error for ri in force_1'

    sr2_ovr      = 1.77 # Overlap threshold (pot > 100)
    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    # Calculate potential at cutoff
    sr2     = 1.0 / r_cut**2 # in sigma=1 units
    sr6     = sr2 ** 3
    sr12    = sr6 **2
    pot_cut = sr12 - sr6 # Without numerical factor 4

    if fast:
        rij = ri - r                    # Get all separation vectors from partners
        rij = rij - np.rint(rij)        # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2,axis=1)  # Squared separations
        in_range = rij_sq < r_cut_box_sq                    # Set mask for within cutoff
        rij_sq   = rij_sq * box_sq                          # Now in sigma=1 units
        rij      = rij * box                                # Now in sigma=1 units
        sr2      = np.where ( in_range, 1.0 / rij_sq, 0.0 ) # (sigma/rij)**2, only if in range
        ovr      = sr2 > sr2_ovr                            # Set flags for any overlaps
        if np.any(ovr):
            partial = PotentialType ( pot=0.0, cut=0.0, vir=0.0, lap=0.0, ovr=True )
            f_partial = np.zeros_like(r)
            return partial, f_partial
        sr6  = sr2 ** 3
        sr12 = sr6 ** 2
        cut  = sr12 - sr6                                # LJ pair potentials (cut but not shifted)
        vir  = cut + sr12                                # LJ pair virials
        pot  = np.where ( in_range, cut - pot_cut, 0.0 ) # LJ pair potential (cut-and-shifted)
        lap  = ( 22.0*sr12 - 5.0*sr6 ) * sr2             # LJ pair Laplacians
        fij  = vir*sr2                                   # LJ scalar part of forces
        f_partial = rij * fij[:,np.newaxis]              # LJ pair forces on i due to each j
        partial = PotentialType ( cut=np.sum(cut), pot=np.sum(pot), vir=np.sum(vir), lap=np.sum(lap), ovr=np.any(ovr) )

    else:
        partial = PotentialType ( pot=0.0, cut=0.0, vir=0.0, lap=0.0, ovr=False )
        f_partial = np.zeros_like(r)
        for j, rj in enumerate(r):
            rij = ri - rj            # Separation vector
            rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
            rij_sq = np.sum(rij**2)  # Squared separation
            if rij_sq < r_cut_box_sq:    # Check within cutoff
                rij_sq = rij_sq * box_sq # Now in sigma=1 units
                rij    = rij * box       # Now in sigma=1 units
                sr2    = 1.0 / rij_sq    # (sigma/rij)**2
                ovr    = sr2 > sr2_ovr   # Overlap if too close
                if ovr:
                    partial.ovr=True
                    return partial, f_partial
                sr6  = sr2 ** 3
                sr12 = sr6 ** 2
                cut  = sr12 - sr6                    # LJ pair potential (cut but not shifted)
                vir  = cut + sr12                    # LJ pair virial
                pot  = cut - pot_cut                 # LJ pair potential (cut-and-shifted)
                lap  = ( 22.0*sr12 - 5.0*sr6 ) * sr2 # LJ pair Laplacian
                f_partial[j,:] = rij * vir * sr2     # LJ pair force on i due to j
                partial = partial + PotentialType ( pot=pot, cut=cut, vir=vir, lap=lap, ovr=ovr )

    # Multiply results by numerical factors
    partial.pot = partial.pot * 4.0
    partial.cut = partial.cut * 4.0
    partial.vir = partial.vir * 24.0 / 3.0
    partial.lap = partial.lap * 24.0 * 2.0
    partial.ovr = False # No overlaps detected (redundant but for clarity)
    f_partial   = f_partial   * 24.0

    return partial, f_partial
