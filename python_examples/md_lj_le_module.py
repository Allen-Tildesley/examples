#!/usr/bin/env python3
# md_lj_le_module.py

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

"""Force routine for MD simulation, Lennard-Jones atoms, Lees-Edwards boundaries."""

fast     = True # Change this to replace NumPy force evaluation with slower Python
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

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def force ( box, strain, r ):
    """Takes in box, strain, and coordinate array, and calculates forces and potentials etc."""

    import numpy as np
    
    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1
    # Lees-Edwards boundaries, in sliding brick arrangement
    # Flow/gradient/vorticity directions are x/y/z == 0/1/2

    n, d = r.shape
    assert d==3, 'Dimension error in force'
    
    sr2_ovr      = 1.77 # Overlap threshold (pot > 100)
    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    # Initialize
    f = np.zeros_like(r)
    total = PotentialType ( pot=0.0, vir=0.0, pyx=0.0, lap=0.0, ovr=False )

    if fast:
        for i in range(n-1):
            rij       = r[i,:]-r[i+1:,:]                         # Separation vectors for j>i
            rij[:,0]  = rij[:,0] - np.rint ( rij[:,1] ) * strain # Extra correction in box=1 units
            rij       = rij - np.rint(rij)                       # Periodic boundary conditions in box=1 units
            rij_sq    = np.sum(rij**2,axis=1)                    # Squared separations for j>1
            in_range  = rij_sq < r_cut_box_sq                    # Set flags for within cutoff
            rij_sq    = rij_sq * box_sq                          # Now in sigma=1 units
            rij       = rij * box                                # Now in sigma=1 units
            sr2       = np.where ( in_range, 1.0/rij_sq, 0.0 )   # (sigma/rij)**2, only if in range
            ovr       = sr2 > sr2_ovr                            # Overlap if too close
            sr6       = sr2 ** 3
            sr12      = sr6 ** 2
            pot       = sr12 - sr6                               # LJ pair potential (cut but not shifted)
            vir       = pot + sr12                               # LJ pair virial
            pot       = np.where ( in_range, pot+0.25, 0.0 )     # WCA LJ pair potential (cut-and-shifted)
            lap       = ( 22.0*sr12 - 5.0*sr6 ) * sr2            # LJ pair Laplacian
            fij       = vir * sr2                                # LJ scalar part of forces
            fij       = rij * fij[:,np.newaxis]                  # LJ pair forces
            pyx       = rij[:,1]*fij[:,0]                        # Off-diagonal element of pressure tensor
            total     = total + PotentialType ( pot=np.sum(pot), vir=np.sum(vir), 
                                                pyx=np.sum(pyx), lap=np.sum(lap), ovr=np.any(ovr) )
            f[i,:]    = f[i,:] + np.sum(fij,axis=0)
            f[i+1:,:] = f[i+1:,:] - fij

    else:
        for i in range(n-1): # Outer loop
            for j in range(i+1,n): # Inner loop
                rij    = r[i,:]-r[j,:]                        # Separation vector
                rij[0] = rij[0] - np.rint ( rij[1] ) * strain # Extra correction in box=1 units
                rij    = rij - np.rint(rij)                   # Periodic boundary conditions in box=1 units
                rij_sq = np.sum(rij**2)                       # Squared separation

                if rij_sq < r_cut_box_sq: # Check within cutoff
                    rij_sq = rij_sq * box_sq # Now in sigma=1 units
                    rij    = rij * box       # Now in sigma=1 units
                    sr2    = 1.0 / rij_sq    # (sigma/rij)**2
                    ovr    = sr2 > sr2_ovr   # Overlap if too close

                    sr6  = sr2 ** 3
                    sr12 = sr6 ** 2
                    pot  = sr12 - sr6                    # LJ pair potential (cut but not shifted)
                    vir  = pot + sr12                    # LJ pair virial
                    pot  = pot + 0.25                    # WCA LJ pair potential (cut-and-shifted)
                    lap  = ( 22.0*sr12 - 5.0*sr6 ) * sr2 # LJ pair Laplacian
                    fij  = rij * vir * sr2               # LJ pair forces
                    pyx  = rij[1]*fij[0]                 # Off-diagonal element of pressure tensor

                    total  = total + PotentialType ( pot=pot, vir=vir, pyx=pyx, lap=lap, ovr=ovr )
                    f[i,:] = f[i,:] + fij
                    f[j,:] = f[j,:] - fij

    # Multiply results by numerical factors
    f         = f         * 24.0       # 24*epsilon
    total.pot = total.pot * 4.0        # 4*epsilon
    total.vir = total.vir * 24.0 / 3.0 # 24*epsilon and divide virial by 3
    total.pyx = total.pyx * 24.0       # 24*epsilon
    total.lap = total.lap * 24.0 * 2.0 # 24*epsilon and factor 2 for ij and ji
    
    return total, f
