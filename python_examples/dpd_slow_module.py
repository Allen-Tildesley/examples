#!/usr/bin/env python3
"""Dissipative particle dynamics module (slow version)."""

class PotentialType:
    """A composite variable for interactions."""

    def __init__(self, pot, vir, lap):
        self.pot = pot # the potential energy cut-and-shifted at r_cut
        self.vir = vir # the virial
        self.lap = lap # the Laplacian

    def __add__(self, other):
        pot = self.pot +  other.pot
        vir = self.vir +  other.vir
        lap = self.lap +  other.lap

        return PotentialType(pot,vir,lap)

def introduction():
    """Prints out introductory statements at start of run."""

    print('DPD soft potential')
    print('Diameter, r_cut = 1')
    print('Slow version built around python loops')

def conclusion():
    """Prints out concluding statements at end of run."""

    print('Program ends')

def make_ij ( box, r ):
    """Compiles randomized list of pairs within range, and stores in the array ij.

    Uses a simple Python loop, which will be slow.
    """

    import numpy as np
    from itertools import combinations

    n,d = r.shape
    assert d==3, 'Dimension error in make_ij'

    ij=[] # Initialize list of pairs
    for i,j in combinations(range(n),2): # Double loop over pairs of atoms
        rij = r[i,:]-r[j,:]      # Separation vector
        rij = rij - np.rint(rij) # Periodic boundary conditions
        rij = rij * box          # Now in r_cut=1 units
        rij_sq = np.sum(rij**2)  # Squared distance
        if rij_sq < 1.0:         # If in range ...
            ij.append([i,j])     # add to list

    ij = np.array(ij)              # Convert to array
    ij = np.random.permutation(ij) # Permute (on first index only)
    return ij

def force ( box, a, r ):
    """Takes in box, strength parameter, and coordinate array, and calculates forces and potentials etc.

    Also returns list of pairs in range for use by the thermalization algorithm.
    Uses a simple Python loop, which will be slow.
    """

    import numpy as np
    from itertools import combinations

    # It is assumed that positions are in units where box = 1

    n,d = r.shape
    assert d==3, 'Dimension error in force'

    # Initialize
    f = np.zeros_like(r)
    total = PotentialType ( pot=0.0, vir=0.0, lap=0.0 )
    pairs = [] # Initialize list of pair information

    for i,j in combinations(range(n),2): # Double loop over pairs of atoms
        rij = r[i,:]-r[j,:]      # Separation vector
        rij = rij - np.rint(rij) # Periodic boundary conditions
        rij = rij * box          # Now in r_cut=1 units
        rij_sq = np.sum(rij**2)  # Squared distance
        if rij_sq < 1.0:         # If in range ...
            rij_mag = np.sqrt(rij_sq) # Distance
            wij     = 1.0 - rij_mag   # Weight function
            rij_hat = rij / rij_mag   # Unit separation vector

            pot = 0.5 * wij**2    # Pair potential
            vir = wij * rij_mag   # Pair virial
            lap = 3.0-2.0/rij_mag # Pair Laplacian
            fij = wij * rij_hat   # Pair force

            total  = total + PotentialType ( pot=pot, vir=vir, lap=lap )
            f[i,:] = f[i,:] + fij
            f[j,:] = f[j,:] - fij
            pairs.append((i,j,rij_hat,rij_mag)) # add to list of pair information

    # Multiply results by numerical factors
    total.pot = total.pot * a
    total.vir = total.vir * a / 3.0
    total.lap = total.lap * a * 2.0
    f         = f * a

    return total, f, pairs

def lowe ( box, temperature, gamma_step, v, pairs ):
    """Updates velocities using pairwise Lowe-Andersen thermostat.

    Uses a simple Python loop, which will be slow.
    """

    # It is assumed that the pairs list contains all pairs within range

    import numpy as np

    v_std = np.sqrt(2*temperature) # Standard deviation for relative velocity distribution

    for p in np.random.permutation(len(pairs)):
        zeta = np.random.rand()
        if zeta<gamma_step:
            (i,j,rij_mag,rij_hat) = pairs[p]
            vij    = v[i,:] - v[j,:]             # Relative velocity vector
            v_old  = np.dot(vij,rij_hat)         # Projection of vij along separation
            v_new  = v_std*np.random.randn()     # New projection of vij along separation
            vij    = ( v_new - v_old ) * rij_hat # Change in relative velocity
            v[i,:] = v[i,:] + 0.5 * vij          # New i-velocity
            v[j,:] = v[j,:] - 0.5 * vij          # New j-velocity

    return v

def shardlow ( box, temperature, gamma_step, v, pairs ):
    """Updates velocities using Shardlow integration algorithm.

    Uses a simple Python loop, which will be slow.
    """

    # It is assumed that the pairs list contains all pairs within range

    import numpy as np

    sqrt_gamma_step = np.sqrt(gamma_step)
    v_std = np.sqrt(2*temperature) # Standard deviation for relative velocity distribution

    for p in np.random.permutation(len(pairs)):
        (i,j,rij_mag,rij_hat) = pairs[p]

        wij       = 1 - rij_mag           # Weight function
        sqrt_prob = sqrt_gamma_step * wij # sqrt of p-factor
        prob      = sqrt_prob**2          # p-factor
        v_new     = v_std*np.random.randn()

        # First half step
        vij    = v[i,:] - v[j,:]                            # Relative velocity vector
        v_old  = np.dot(vij,rij_hat)                        # Projection of vij along separation
        vij    = ( sqrt_prob*v_new - prob*v_old ) * rij_hat # Change in relative velocity
        v[i,:] = v[i,:] + 0.5 * vij                         # New i-velocity
        v[j,:] = v[j,:] - 0.5 * vij                         # New j-velocity

        # Second half step
        vij    = v[i,:] - v[j,:]                                     # Relative velocity vector
        v_old  = np.dot(vij,rij_hat)                                 # Projection of vij along separation
        vij    = ( sqrt_prob*v_new - prob*v_old ) * rij_hat/(1+prob) # Change in relative velocity
        v[i,:] = v[i,:] + 0.5 * vij                                  # New i-velocity
        v[j,:] = v[j,:] - 0.5 * vij                                  # New j-velocity

    return v

def p_approx ( a, rho, temperature ):
    """Returns approximate pressure."""

    # This expression is given by Groot and Warren, J Chem Phys 107, 4423 (1997), alpha = 0.101
    # This is the revised formula due to Liyana-Arachchi, Jamadagni, Elke, Koenig, Siepmann,
    # J Chem Phys 142, 044902 (2015)

    import numpy as np

    c1, c2 = 0.0802, 0.7787
    b = [1.705e-8,-2.585e-6,1.556e-4,-4.912e-3,9.755e-2]
    b2 = np.polyval (b,a)                                           # This is B2/a, eqn (10) of above paper
    alpha = b2 / ( 1 + rho**3 ) + ( c1*rho**2 ) / ( 1 + c2*rho**2 ) # This is eqn (14) of above paper
    return rho * temperature + alpha * a * rho**2
