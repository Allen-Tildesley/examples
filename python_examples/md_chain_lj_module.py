#!/usr/bin/env python3
# md_chain_lj_module.py

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

"""Force & constraint routines for MD, LJ chain."""

fast = True # Change this to replace NumPy force evaluation with slower Python

class PotentialType:
    """A composite variable for interactions."""

    def __init__(self, pot, ovr):
        self.pot = pot # the potential energy
        self.ovr = ovr # a flag indicating overlap (i.e. pot too high to use)

    def __add__(self, other):
        pot = self.pot +  other.pot
        ovr = self.ovr or other.ovr
        return PotentialType(pot,ovr)

def introduction():
    """Prints out introductory statements at start of run."""

    print('LJ chain, no cutoff, no shift, no periodic box')
    print('Diameter, sigma = 1')
    print('Well depth, epsilon = 1')
    print('All atomic masses the same m = 1')
    print('All bond lengths the same')
    if fast:
        print('Fast NumPy force routine')
    else:
        print('Slow Python force routine')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def zero_cm ( r, v ):
    """Routine to set centre-of-mass at the origin and zero the total momentum."""

    import numpy as np

    n, d = r.shape
    assert d==3, 'r dimension error in zero_cm'
    cr = np.sum ( r, axis=0 ) / n
    n, d = v.shape
    assert d==3, 'v dimension error in zero_cm'
    cv = np.sum ( v, axis=0 ) / n
    return r-cr, v-cv

def worst_bond ( bond, r ):
    """Returns max amount by which constraint is violated."""

    import numpy as np
    
    rij      = r[:-1,:]-r[1:,:]               # All nearest neighbour separation vectors
    rij_mag  = np.sqrt(np.sum(rij**2,axis=1)) # Current bond lengths
    diff     = np.fabs(rij_mag-bond)          # Absolute amounts by which constraint is violated
    return np.max(diff)                       # Find maximum

def force ( r ):
    """Takes in coordinate array, returns forces, potential etc."""

    import numpy as np

    # The Lennard-Jones energy and sigma parameters are taken to be epsilon = 1, sigma = 1
    # Positions are assumed to be in these units
    # Forces are calculated in the same units and stored in the array f
    # NO box, NO periodic boundaries

    n, d = r.shape
    assert d==3, 'Dimension error in force'

    sr2_ovr = 1.77 # Overlap threshold (pot > 100)

    # Initialize
    f     = np.zeros_like(r)
    total = PotentialType (  pot=0.0, ovr=False )

    if fast:
        for i in range(n-2):
            rij = r[i,:]-r[i+2:,:]         # Separation vectors for j>i+1
            rij_sq = np.sum(rij**2,axis=1) # Squared separations for j>1
            sr2    = 1.0 / rij_sq          # (sigma/rij)**2
            ovr    = sr2 > sr2_ovr         # Overlap if too close
            sr6    = sr2 ** 3
            sr12   = sr6 ** 2
            pot    = sr12 - sr6              # LJ pair potential
            vir    = pot + sr12              # LJ pair virial
            fij    = vir * sr2               # LJ scalar part of forces
            fij    = rij * fij[:,np.newaxis] # LJ pair forces
            total     = total + PotentialType ( pot=np.sum(pot), ovr=np.any(ovr) )
            f[i,:]    = f[i,:] + np.sum(fij,axis=0)
            f[i+2:,:] = f[i+2:,:] - fij
    else:
        for i in range(n-2): # Outer loop
            for j in range(i+2,n): # Inner loop, skipping nearest neighbour
                rij    = r[i,:]-r[j,:]  # Separation vector
                rij_sq = np.sum(rij**2) # Squared separation
                sr2    = 1.0 / rij_sq   # (sigma/rij)**2
                ovr    = sr2 > sr2_ovr  # Overlap if too close
                sr6    = sr2 ** 3
                sr12   = sr6 ** 2
                pot    = sr12 - sr6      # LJ pair potential
                vir    = pot + sr12      # LJ pair virial
                fij    = rij * vir * sr2 # LJ Pair forces

                total  = total  + PotentialType ( pot=pot, ovr=ovr )
                f[i,:] = f[i,:] + fij
                f[j,:] = f[j,:] - fij

    # Multiply results by numerical factors
    f         = f         * 24.0       # 24*epsilon
    total.pot = total.pot * 4.0        # 4*epsilon
    
    return total, f

def spring ( k_spring, bond, r ):
    """Calculates bond spring potential and forces for atomic chain."""

    import numpy as np
    
    # NO box, NO periodic boundaries

    n, d = r.shape
    assert d==3, 'Dimension error in spring'

    # Initialize
    g         = np.zeros_like(r)
    total_spr = 0.0

    if fast:
        rij       = r[:-1,:]-r[1:,:]             # All nearest neighbour separation vectors
        rij_sq    = np.sum(rij**2,axis=1)        # Squared separations
        rij_mag   = np.sqrt(rij_sq)              # Separations
        pair_pot  = (rij_mag-bond)**2            # Spring pair potentials without numerical factor
        gij       = ( bond - rij_mag ) / rij_mag # Factor determining magnitude of forces
        gij       = rij * gij[:,np.newaxis]      # Spring pair forces without numerical factor
        total_spr = total_spr + np.sum(pair_pot)
        g[:-1,:]  = g[:-1,:] + gij
        g[1:,:]   = g[1:,:]  - gij
    else:
        for i in range(n-1):                               # Loop over atoms
            j = i+1                                        # Nearest neighbour
            rij       = r[i,:] - r[j,:]                    # Separation vector
            rij_sq    = np.sum(rij**2)                     # Squared separation
            rij_mag   = np.sqrt(rij_sq)                    # Separation
            pair_pot  = (rij_mag-bond)**2                  # Spring pair potential without numerical factor
            gij       = rij * ( bond - rij_mag ) / rij_mag # Spring pair force without numerical factor
            total_spr = total_spr + pair_pot
            g[i,:]    = g[i,:] + gij
            g[j,:]    = g[j,:] - gij

    # Multiply results by numerical factors
    total_spr = total_spr * 0.5 * k_spring
    g         = g * k_spring

    return total_spr, g

def rattle_a ( dt, bond, r_old, r, v ):
    """First part of velocity Verlet algorithm with constraints."""

    import numpy as np
    
    # This subroutine iteratively adjusts the positions stored in the array r
    # and the velocities stored in the array v, to satisfy the bond constraints

    # On entry to this routine we assume:
    # r_old stores the positions at the start of the step
    # r stores the positions following the unconstrained drift and
    # v stores the velocities following the first unconstrained half-kick
    # The returned arrays r and v will hold the constrained values

    n, d = r.shape
    assert d==3, 'Dimension error in rattle_a'

    tol = 1.0e-9
    tol2 = 2.0 * tol
    dot_tol = 1.0e-9
    iter_max = 500

    iter  = 0
    done  = False
    moved = np.full(n,True,dtype=np.bool_) # Ensures that we look at each bond at least once
    move  = np.empty_like(moved)

    while True: # Iterative loop until done

        if done:
            break # done is equivalent to not np.any ( moved )

        done = True
        move[:] = False

        for i in range(n-1): # Loop over each constraint in turn
            j = i + 1 # Partner atom in this constraint

            if moved[i] or moved[j]: # Test whether need to re-examine ij
                rij = r[i,:] - r[j,:]             # Current bond vector
                diffsq = bond**2 - np.sum(rij**2) # Amount by which constraint is violated

                if abs(diffsq) > tol2*bond**2: # Test whether constraint not already satisfied

                    rij_old = r_old[i,:] - r_old[j,:] # Old vector determines direction of constraint force
                    dot = np.dot(rij_old,rij) # This should be of the order of bond**2

                    assert dot > dot_tol*bond**2, "{}{:15.6f}".format('Constraint failure',dot)

                    # In the following formulae, inverse masses are all unity
                    g       = diffsq / ( 4.0 * dot ) 
                    dr      = rij_old * g    # Position adjustment
                    r[i,:]  = r[i,:] + dr    # Adjust i position
                    r[j,:]  = r[j,:] - dr    # Adjust j position
                    v[i,:]  = v[i,:] + dr/dt # Adjust i velocity
                    v[j,:]  = v[j,:] - dr/dt # Adjust j velocity
                    move[i] = True           # Flag that we moved i
                    move[j] = True           # Flag that we moved j
                    done    = False          # Flag that we moved something

        # Prepare for next iteration
        moved = move.copy()
        iter  = iter + 1
        assert iter <= iter_max, "{}{:15d}{:15d}".format('Too many iterations', iter, iter_max)

    return r, v

def rattle_b ( dt, bond, r, v ):
    """Second part of velocity Verlet with constraints."""

    import numpy as np
    
    # This subroutine iteratively adjusts the velocities stored in the array v
    # to satisfy the time derivatives of the bond constraints
    # Also returns constraint contribution to virial

    # On entry to this routine we assume:
    # r stores the positions at the end of the step with constraints applied
    # v stores the velocities following the second unconstrained half-kick
    # The returned array v will hold the constrained values

    n, d = r.shape
    assert d==3, 'r dimension error in rattle_b'

    tol = 1.0e-9
    tol2 = 2.0 * tol
    dot_tol = 1.0e-9
    iter_max = 500

    iter  = 0
    done  = False
    moved = np.full(n,True,dtype=np.bool_) # Ensures that we look at each bond at least once
    move  = np.empty_like(moved)
    wc    = 0.0

    while True: # Iterative loop until done

        if done:
            break

        done = True
        move[:] = False

        for i in range(n-1): # Loop over each constraint in turn
            j = i + 1 # Partner atom in this constraint

            if moved[i] or moved[j]: # Test whether need to re-examine ij
                vij = v[i,:] - v[j,:]
                rij = r[i,:] - r[j,:]
                dot = np.dot ( rij, vij )

                # In the following formulae, inverse masses are all unity
                g  = -dot / ( 2.0 * bond**2 )

                if abs(g) > tol: # Test whether constraint already satisfied

                    wc      = wc + g * bond**2 # Contribution to virial
                    dv      = rij * g          # Velocity adjustment
                    v[i,:]  = v[i,:] + dv      # Adjust velocity i
                    v[j,:]  = v[j,:] - dv      # Adjust velocity j
                    move[i] = True             # Flag that we moved i
                    move[j] = True             # Flag that we moved j
                    done    = False            # Flag that we moved something

        # Prepare for next iteration
        moved = move.copy()
        iter  = iter + 1
        assert iter <= iter_max, "{}{:15d}{:15d}".format('Too many iterations', iter, iter_max)

    wc = wc / (0.5*dt) / 3.0 # Scale factors for virial

    return wc, v

def milcshake_a ( dt, bond, r_old, r, v ):
    """First part of velocity Verlet algorithm with constraints."""

    # This subroutine iteratively adjusts the positions stored in the array r
    # and the velocities stored in the array v, to satisfy the bond constraints
    # using a tri-diagonal solver
    # See AG Bailey, CP Lowe, and AP Sutton, J Comput Phys, 227, 8949 (2008)
    # and AG Bailey, CP Lowe, and AP Sutton, Comput Phys Commun, 180, 594 (2009)

    # On entry to this routine we assume:
    # r_old stores the positions at the start of the step
    # r stores the positions following the unconstrained drift and
    # v stores the velocities following the first unconstrained half-kick
    # The returned arrays r and v will hold the constrained values

    import numpy as np
    import scipy.linalg as la

    n, d = r.shape
    assert d==3, 'r dimension error in milcshake_a'

    k = n-1 # Number of constraints
    
    tol = 1.0e-9
    iter_max = 500

    r_new = r.copy() # Saves unconstrained positions

    # Old and new (non-constrained) bond vectors
    rij_old = r_old[:-1,:] - r_old[1:,:]
    rij_new = r_new[:-1,:] - r_new[1:,:]

    # Elements of tridiagonal matrix (dot products of old and new bond vectors)
    # In this example, all masses are equal to unity. Let k=n-1 be number of constraints
    tridiag = np.zeros((3,k), dtype=np.float_)
    tridiag[0,1:]  = -2.0*np.sum ( rij_old[1:,:] *rij_new[:-1,:], axis=1 )       # leading zero to pad, then k-1 elements of upper-diagonal
    tridiag[1,:]   =  2.0*np.sum ( rij_old[:,:]  *rij_new[:,:],   axis=1 ) / 0.5 # k elements of diagonal
    tridiag[2,:-1] = -2.0*np.sum ( rij_old[:-1,:]*rij_new[1:,:],  axis=1 )       # k-1 elements of lower-diagonal, then trailing zero to pad

    # Set up rhs of constraint equation
    rijsq  = np.sum(rij_new**2,axis=1)
    rhs    = bond**2 - rijsq
    rhsold = rhs.copy()

    iter = 0

    while True: # Iterative loop until done

        # Test for done
        max_error = np.max(np.fabs(rijsq-bond**2))/(2.0*bond**2)
        if max_error <= tol:
            break

        # Reset tridiagonal elements (may have been over-written by solver)
        tridiag_tmp = tridiag.copy()
        lam = la.solve_banded((1,1),tridiag_tmp,rhs) 

        # Constraint effects on position from lambda multipliers
        r = r_new.copy()
        r[:-1,:] = r[:-1,:] + lam[:,np.newaxis]*rij_old
        r[1:,:]  = r[1:,:]  - lam[:,np.newaxis]*rij_old

        # New bond vectors
        rij = r[:-1,:] - r[1:,:]

        # Prepare for next iteration
        rijsq  = np.sum(rij**2,axis=1)
        rhs    = bond**2 - rijsq + rhsold
        rhsold = rhs.copy()    

        iter = iter + 1
        assert iter <= iter_max, "{}{:15d}{:15d}".format('Too many iterations', iter, iter_max)

    # Effect of constraints on velocities
    v[:-1,:] = v[:-1,:] + lam[:,np.newaxis]*rij_old/dt
    v[1:,:]  = v[1:,:]  - lam[:,np.newaxis]*rij_old/dt

    return r, v

def milcshake_b ( dt, bond, r, v ):
    """Second part of velocity Verlet algorithm with constraints."""

    # This subroutine adjusts the velocities stored in the array v
    # to satisfy the time derivatives of the bond constraints
    # using a tri-diagonal solver: here we use dgtsv from LAPACK.
    # See AG Bailey, CP Lowe, and AP Sutton, J Comput Phys, 227, 8949 (2008)
    # and AG Bailey, CP Lowe, and AP Sutton, Comput Phys Commun, 180, 594 (2009)
    # Also returns constraint contribution to virial

    # On entry to this routine we assume:
    # r stores the positions at the end of the step with constraints applied
    # v stores the velocities following the second unconstrained half-kick
    # The returned array v will hold the constrained values

    import numpy as np
    import scipy.linalg as la

    n, d = r.shape
    assert d==3, 'r dimension error in milcshake_b'

    k = n-1 # Number of constraints

    # Relative velocities and bond vectors
    vij = v[:-1,:] - v[1:,:]
    rij = r[:-1,:] - r[1:,:]
    rhs = -np.sum(vij*rij,axis=1)

    # Elements of tridiagonal matrix (dot products of bond vectors)
    # In this example, all masses are equal to unity. Let k=n-1 be number of constraints
    tridiag = np.zeros((3,k), dtype=np.float_)
    tridiag[0,1:]  = -np.sum ( rij[1:,:] *rij[:-1,:], axis=1 )       # leading zero to pad, then k-1 elements of upper-diagonal
    tridiag[1,:]   =  np.sum ( rij[:,:]  *rij[:,:],   axis=1 ) / 0.5 # k elements of diagonal
    tridiag[2,:-1] = -np.sum ( rij[:-1,:]*rij[1:,:],  axis=1 )       # k-1 elements of lower-diagonal, then trailing zero to pad

    lam = la.solve_banded((1,1),tridiag,rhs) 

    # Effect of constraints on velocities
    v[:-1,:] = v[:-1,:] + lam[:,np.newaxis]*rij
    v[1:,:]  = v[1:,:]  - lam[:,np.newaxis]*rij

    wc = np.sum(lam) * bond**2
    wc = wc / (0.5*dt) / 3.0 # scale factors for virial

    return wc, v
