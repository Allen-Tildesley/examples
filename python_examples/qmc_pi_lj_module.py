#!/usr/bin/env python3
# qmc_pi_lj_module.py

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

"""Energy and move routines for PIMC simulation, LJ potential."""

fast = True # Change this to replace NumPy potential evaluation with slower Python

class PotentialType:
    """A composite variable for interactions."""

    def __init__(self, pot, vir, ovr):
        self.pot = pot # the potential energy cut (but not shifted) at r_cut
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

    p, n, d = r.shape
    assert d==3, 'Dimension error for r'

    total = PotentialType ( pot=0.0, vir=0.0, ovr=False )

    for k in range(p):
        for i in range(n-1):
            partial = potential_1 ( r[k,i,:], box, r_cut, r[k,i+1:,:], p )
            if partial.ovr:
                total.ovr = True
                return total
            total = total + partial

    return total

def potential_1 ( rki, box, r_cut, r, p ):
    """Takes in coordinates of an atom and calculates its interactions.

    Values of box, cutoff range, and partner coordinate array are supplied.
    The results are returned as partial, a PotentialType variable.
    """

    import numpy as np

    # partial.pot is the nonbonded cut (not shifted) potential energy of atom rki with a set of other atoms
    # partial.vir is the corresponding virial of atom rki
    # partial.lap is the corresponding Laplacian of atom rki
    # partial.ovr is a flag indicating overlap (potential too high) to avoid overflow
    # If this is True, the values of partial.pot etc should not be used
    # In general, r will be a subset of the complete set of simulation coordinates
    # and none of its rows should be identical to rki

    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1

    nj, d = r.shape
    assert d==3, 'Dimension error for r in potential_1'
    assert rki.size==3, 'Dimension error for rki in potential_1'

    sr2_ovr      = 1.77 # Overlap threshold (pot > 100)
    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    if fast:
        r_ki_kj    = rki - r                    # Get all separation vectors from partners
        r_ki_kj    = r_ki_kj - np.rint(r_ki_kj) # Periodic boundary conditions in box=1 units
        r_ki_kj_sq = np.sum(r_ki_kj**2,axis=1)  # Squared separations
        in_range   = r_ki_kj_sq < r_cut_box_sq  # Set flags for within cutoff
        r_ki_kj_sq = r_ki_kj_sq * box_sq        # Now in sigma=1 units
        sr2        = np.where ( in_range, 1.0 / r_ki_kj_sq, 0.0 ) # (sigma/r_ki_kj)**2, only if in range
        ovr        = sr2 > sr2_ovr                                # Set flags for any overlaps
        if np.any(ovr):
            partial = PotentialType ( pot=0.0, vir=0.0, ovr=True )
            return partial
        sr6  = sr2 ** 3
        sr12 = sr6 ** 2
        pot  = sr12 - sr6                    # LJ pair potentials (cut but not shifted)
        vir  = pot + sr12                    # LJ pair virials
        partial = PotentialType ( pot=np.sum(pot), vir=np.sum(vir), ovr=False )

    else:
        partial = PotentialType ( pot=0.0, vir=0.0, ovr=False )
        for rkj in r:
            r_ki_kj = rki - rkj                  # Separation vector
            r_ki_kj = r_ki_kj - np.rint(r_ki_kj) # Periodic boundary conditions in box=1 units
            r_ki_kj_sq = np.sum(r_ki_kj**2)      # Squared separation
            if r_ki_kj_sq < r_cut_box_sq:        # Check within cutoff
                r_ki_kj_sq = r_ki_kj_sq * box_sq # Now in sigma=1 units
                sr2    = 1.0 / r_ki_kj_sq        # (sigma/r_ki_kj)**2
                ovr    = sr2 > sr2_ovr           # Overlap if too close
                if ovr:
                    partial.ovr=True
                    return partial
                sr6  = sr2 ** 3
                sr12 = sr6 ** 2
                pot  = sr12 - sr6                    # LJ pair potential (cut but not shifted)
                vir  = pot + sr12                    # LJ pair virial
                partial = partial + PotentialType ( pot=pot, vir=vir, ovr=ovr )

    # Include numerical factors
    partial.pot = partial.pot * 4.0 / p      # Classical potentials are weaker by a factor p
    partial.vir = partial.vir * 24.0 / (3*p) # 24*epsilon and divide virial by 3 & by p

    return partial

def spring ( box, k_spring, r ):
    """Takes in box, spring strength, and coordinate array, and calculates quantum spring potential."""

    import numpy as np

    # Actual calculation performed by function spring_1

    p, n, d = r.shape
    assert d==3, 'Dimension error for r in potential'

    total = 0.0
    
    for i in range(n): # Loop over ring polymers
        for k in range(p): # Loop over atoms within polymer
            kp = (k+1)%p # Neighbour index
            partial = spring_1 ( r[k,i,:], r[kp,i,:], box, k_spring )
            total = total + partial

    return total

def spring_1 ( rki, rli, box, k_spring ):
    """Returns quantum potential for given atom."""

    import numpy as np

    # Coordinate array of neighbour is supplied
    
    r_ki_li    = rki - rli                   # Separation vector
    r_ki_li    = r_ki_li - np.rint(r_ki_li)  # Periodic boundary conditions in box=1 units
    r_ki_li_sq = np.sum(r_ki_li**2) * box**2 # Squared separation in sigma=1 units

    return 0.5 * k_spring * r_ki_li_sq
