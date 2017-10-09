#!/usr/bin/env python3
# md_poly_lj_module.py

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

"""Force routine for MD simulation, polyatomic molecules, LJ atoms."""

import numpy as np

fast = True # Change this to replace NumPy force evaluation with slower Python

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

m = np.array([ 1.0/3.0, 1.0/3.0, 1.0/3.0 ],dtype=np.float_) # Masses add up to 1.0

# The following section sets the diagonal moments of inertia "realistically"
# based on the values of atomic masses and bond vectors (above), with some checking
# that the total mass is 1, the COM is at the origin, and the inertia tensor is diagonal.
# However, there is nothing to stop the user replacing this section with a statement just
# setting the values of inertia[:]. The masses m are not used by the calling program.
# It might be advantageous, for instance, to artificially increase the values in inertia.
    
# Ensure that the db bonds, xyz molecular axes, and masses are chosen such that
# the total mass is 1 and the centre-of-mass is at the origin
assert np.isclose(np.sum(m),1.0), 'M is not 1.0 {}'.format(np.sum(m))
com = np.sum ( m[:,np.newaxis]*db, axis = 0 )
assert np.all ( np.isclose(com,0.0) ), 'COM error {} {} {}'.format(*com)

# Ensure that the db bonds, xyz molecular axes, and masses are chosen such that
# the off-diagonal elements of the inertia tensor are zero
inertia = -np.einsum('ij,ik->jk',m[:,np.newaxis]*db,db)
offdiag = inertia[np.triu_indices(3,1)]
assert np.all ( np.isclose(offdiag,0.0) ), 'Inertia not diagonal {} {} {}'.format(*offdiag)

# Calculate the diagonal elements of the inertia tensor
inertia = np.sum(m[:,np.newaxis]*db**2) + np.diagonal ( inertia )

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

def introduction():
    """Prints out introductory statements at start of run."""
    
    print('Lennard-Jones potential')
    print('Cut-and-force-shifted')
    print('Diameter, sigma = 1')
    print('Well depth, epsilon = 1')
    if fast:
        print('Fast NumPy force routine')
    else:
        print('Slow Python force routine')

    print( "{:40}{:15d}".format('Number of atoms per molecule', na) )
    for i, b in enumerate(db):
        print( "{}{:2d}{:15.6f}{:15.6f}{:15.6f}".format('Body-fixed atom vector',i,*b))
    print( "{:40}{:15.6f}".format('Molecular diameter', diameter) )

    print( "{:40}{:15.6f}".format('r_cut', r_cut) )
    print( "{:40}{:15.6f}".format('Force-shift lambda1', lambda1) )
    print( "{:40}{:15.6f}".format('Force-shift lambda2', lambda2) )

    print( "{:40}{:15.6f}{:15.6f}{:15.6f}".format('Inertia Ixx, Iyy, Izz', *inertia) )

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def force ( box, r, d ):
    """Takes in box, and r & d arrays, and calculates forces, torques and potentials etc."""

    import numpy as np
    
    # It is assumed that positions are in units where box = 1
    # Forces are calculated in units where sigma = 1 and epsilon = 1
    # Note that this is the force-shifted LJ potential with a linear smoothing term
    # S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)

    n, ndim = r.shape
    assert ndim==3, 'Dimension error for r'
    nn, nna, ndim = d.shape
    assert nna==na and ndim==3, 'Dimension error for d'
    assert n==nn, 'Dimension mismatch for r and d'
    
    sr2_ovr       = 1.77                       # Overlap threshold (pot > 100)
    rm_cut_box    = ( r_cut + diameter ) / box # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box**2              # squared
    assert rm_cut_box<0.5, 'rm_cut/box too large'
    r_cut_sq = r_cut ** 2

    # Initialize
    f     = np.zeros_like(r)
    tau   = np.zeros_like(r)
    total = PotentialType ( pot=0.0, vir=0.0, ovr=False )

    if fast:
        for i in range(n-1):
            rij       = r[i,:]-r[i+1:,:]   # Separation vectors for j>i
            rij       = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
            rij       = rij * box          # Now in sigma=1 units
            for a in range(na):
                for b in range(na):
                    rab      = rij + d[i,a,:] - d[i+1:,b,:] # All atom-atom vectors for given a and b
                    rab_sq   = np.sum(rab**2,axis=1)        # Squared separations
                    in_range = rab_sq < r_cut_sq            # Set flags for within cutoff
                    sr2      = 1.0 / rab_sq                 # (sigma/rab)**2
                    ovr      = sr2 > sr2_ovr                # Set flags for any overlaps
                    rmag     = np.sqrt(rab_sq)
                    sr6      = sr2 ** 3
                    sr12     = sr6 ** 2
                    pot      = np.where ( in_range,
                        4.0*(sr12-sr6) + lambda1 + lambda2*rmag, 0.0 ) # force-shifted pair potentials
                    virab    = np.where ( in_range,
                        24.0*(2.0*sr12-sr6) - lambda2*rmag, 0.0 ) # pair virials
                    fab      = virab * sr2
                    fab      = rab * fab[:,np.newaxis] # atom-atom pair forces

                    total = total + PotentialType ( pot=np.sum(pot), vir=np.sum(rij*fab), ovr=np.any(ovr) )
                    fia         = np.sum(fab,axis=0)
                    f[i,:]      = f[i,:]    + fia
                    f[i+1:,:]   = f[i+1:,:] - fab
                    tau[i,:]    = tau[i,:]    + np.cross ( d[i,a,:],    fia )
                    tau[i+1:,:] = tau[i+1:,:] - np.cross ( d[i+1:,b,:], fab )

    else:
        for i in range(n-1): # Outer loop
            for j in range(i+1,n): # Inner loop
                rij = r[i,:]-r[j,:]      # Separation vector
                rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
                rij_sq = np.sum(rij**2)  # Squared separation

                if rij_sq < rm_cut_box_sq: # Check within cutoff
                    rij    = rij * box       # Now in sigma=1 units
                    for a in range(na):
                        for b in range(na):
                            rab      = rij + d[i,a,:] - d[j,b,:] # Atom-atom vector for given a and b
                            rab_sq   = np.sum(rab**2)            # Squared separation
                            if rab_sq < r_cut_sq:                # Test within potential cutoff
                                sr2      = 1.0 / rab_sq          # (sigma/rab)**2
                                ovr      = sr2 > sr2_ovr         # Set flag for overlap
                                rmag     = np.sqrt(rab_sq)
                                sr6      = sr2 ** 3
                                sr12     = sr6 ** 2
                                pot      = 4.0*(sr12-sr6) + lambda1 + lambda2*rmag # force-shifted pair potential
                                virab    = 24.0*(2.0*sr12-sr6) - lambda2*rmag      # pair virial
                                fab      = virab * sr2
                                fab      = rab * fab # atom-atom pair force

                                total    = total + PotentialType ( pot=pot, vir=np.sum(rij*fab), ovr=ovr )
                                f[i,:]   = f[i,:] + fab
                                f[j,:]   = f[j,:] - fab
                                tau[i,:] = tau[i,:] + np.cross ( d[i,a,:], fab )
                                tau[j,:] = tau[j,:] - np.cross ( d[j,b,:], fab )

    # Multiply results by numerical factors
    total.vir = total.vir / 3.0 # Divide virial by 3
    
    return total, f, tau
