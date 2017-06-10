#!/usr/bin/env python3
# initialize.py

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

def fcc_positions ( n, box, length, soft, quaternions ):
    """Sets up the fcc lattice: four molecules per unit cell."""

    import numpy as np
    from itertools import product

    # Arguments are the number of particles, box length, linear molecule length,
    # a flag for soft interactions (no overlap check)
    # and a flag to indicate quaternion parameters for nonlinear molecules

    # For atoms, for which length=0.0, the e-coordinates are set, but will be ignored
    # For linear molecules, the orientations comply with the alpha-fcc pattern
    # For nonlinear molecules, the 0-element is set to zero

    print('Close-packed fcc lattice positions')

    nc = np.rint ( (n/4)**(1.0/3.0) ).astype(np.int)
    assert n==4*nc**3, "{}{:d}{:d}".format('n, nc mismatch ',n,4*nc**3)
    cell = box / nc  # Unit cell
    box2 = box / 2.0 # Half box length
    r = np.empty((n,3),dtype=np.float_)
    e = np.empty((n,3),dtype=np.float_)

    r_fcc = np.array ( [ [0.25,0.25,0.25],[0.25,0.75,0.75],[0.75,0.75,0.25],[0.75,0.25,0.75] ], dtype=np.float_ )
    e_fcc = np.array ( [ [1.0,1.0,1.0],[1.0,-1.0,-1.0],[-1.0,1.0,-1.0],[-1.0,-1.0,1.0] ],dtype=np.float_)*np.sqrt(1.0/3.0)

    i = 0
    
    for ix, iy, iz in product(range(nc),repeat=3): # triple loop over unit cells
        for a in range(4): # loop over atoms in unit cell
            r[i,:] = r_fcc[a,:] + np.array ( [ix,iy,iz] ).astype(np.float_) # in range 0..nc
            r[i,:] = r[i,:] * cell                                          # in range 0..box
            r[i,:] = r[i,:] - box2                                          # in range -box2..box2
            e[i,:] = e_fcc[a]
            if not soft:
                assert not overlap ( r[i,:], e[i,:], r[:i,:], e[:i,:], box, length ), "Density too high"
            i = i + 1

    if quaternions:
        e=np.insert(e,0,0.0,axis=1) # insert column 0, full of zeros

    return r, e

def ran_positions ( n, box, length, soft, quaternions ):
    """Places atoms at random positions."""

    import numpy as np
    from maths_module import random_quaternion, random_vector
    
    # Unlikely to be useful, unless the interaction potential is soft
    # or the density rather low
    # For atoms, for which length=0.0, the e-coordinates will be ignored

    iter_max = 10000 # Max random placement iterations

    print('Random positions')

    r = np.empty((n,3),dtype=np.float_)
    if quaternions:
        e = np.empty((n,4),dtype=np.float_)
    else:
        e = np.empty((n,3),dtype=np.float_)

    for i in range(n):
        
        iter = 0
        while True: # Loop until non-overlapping position found
            r[i,:] = ( np.random.rand(3) - 0.5 ) * box # In range -box/2..box/2
            if quaternions:
                e[i,:] = random_quaternion()
            else:
                e[i,:] = random_vector()
            if soft:
                break
            if not overlap ( r[i,:], e[i,:], r[:i,:], e[:i,:], box, length ):
                break

            iter = iter + 1
            assert iter <= iter_max, "Too many iterations"

    return r, e

def ran_velocities ( nn, e, temperature, inertia, quaternions ):
    """Chooses translational velocities from Maxwell-Boltzmann (Gaussian) distribution."""

    import numpy as np
    from maths_module import random_perpendicular_vector

    # We set the total momentum to zero
    # We assume unit molecular mass

    # For linear molecules we choose the direction of the angular velocity
    # randomly but perpendicular to the molecular axis.
    # The square of the magnitude of the angular velocity
    # is chosen from an exponential distribution
    # For nonlinear molecules we choose all three components of angular velocity
    # from a Gaussian distribution, assuming equal moments of inertia
    # There is no attempt to set the total angular momentum to zero
    # For atoms, the w array is set here, but ignored later

    print("{:40}{:15.6f}{:15.6f}".format('Velocities at temperature, inertia', temperature, inertia) )

    n, d = e.shape
    assert n==nn, "e shape mismatch {:d}{:d}".format(n,nn)
    if quaternions:
        assert d==4,  "e shape mismatch {:d}".format(d)
    else:
        assert d==3,  "e shape mismatch {:d}".format(d)

    # Linear velocities
    
    v      = np.random.randn ( n, 3 )                  # Unit normal random numbers
    v_cm   = np.sum ( v, axis=0 ) / n                  # Compute centre of mass velocity
    v      = v - v_cm                                  # Set net momentum to zero
    factor = np.sqrt((3*n-3)*temperature/np.sum(v**2)) # sqrt of ratio of kinetic energies
    v      = factor * v

    # Angular velocities

    if quaternions: # Nonlinear molecule, treat as spherical top
        w_std_dev = np.sqrt(temperature/inertia)
        w = np.random.randn ( n, 3 ) * w_std_dev
    else:
        w_sq_mean = 2.0 * temperature / inertia
        w = np.empty ( (n,3), dtype=np.float_ )
        for i in range(n):
            w[i,:] = random_perpendicular_vector ( e[i,:] ) # Set direction of angular velocity
            w[i,:] = w[i,:] * np.sqrt(np.random.exponential(w_sq_mean))

    return v, w

def chain_positions ( n, bond, soft ):
    """Chooses chain positions randomly, at desired bond length, avoiding overlap."""
    
    import numpy as np
    from maths_module import random_vector

    tol = 1.0e-9
    iter_max = 500

    print("{:40}{:15.6f}".format('Chain, randomly oriented bonds = ',bond) )

    r = np.empty ( (n,3), dtype=np.float_ )

    r[0,:] = [0.0,0.0,0.0] # First atom at origin
    r[1,:] = bond*random_vector() # Second atom at random position (bond length away)

    for i in range(2,n): # Loop over atom indices

        iter = 0
        while True: # Loop until non-overlapping position found
            r[i,:] = r[i-1,:] + bond*random_vector() # Subsequent atoms randomly placed (bond length away)
            if soft: # No overlap test
                break
            # Overlap test on all so far except bonded neighbour
            if not chain_overlap ( r[i,:], r[:i-1,:] ):
                break
            iter = iter + 1
            assert iter <= iter_max, 'Too many iterations'
            
    r_cm = np.sum ( r, axis=0 ) / n # Compute centre of mass
    r    = r - r_cm                 # Shift centre of mass to origin

    for i in range(n-1):
        diff_sq = np.sum ( (r[i,:]-r[i+1,:])**2 ) - bond**2
        if np.fabs(diff_sq)> tol:
            print( "{}{:5d}{:5d}{:15.8f}".format('Bond length warning',i,i+1,diff_sq) )

    return r

def chain_velocities ( nn, temperature, constraints, r ):
    """Chooses velocities from Maxwell-Boltzmann (Gaussian) distribution."""

    import numpy as np
    import scipy.linalg as la

    # For simplicity, we just pick each atom velocity randomly and
    # apply bond constraints afterwards
    # In between, we take steps to remove linear and angular momentum
    # since the configuration will be used in MD simulations without periodic boundaries
    # in which case both these quantities are conserved
    # NB there is at present no check for a singular inertia tensor in the angular momentum fix!
    # We assume centre of mass is already at the origin
    # We assume unit molecular mass and employ Lennard-Jones units
    # property                  units
    # energy                    epsilon ( = 1 )
    # molecular mass            m ( = 1 )
    # velocity v                sqrt(epsilon/m)

    tol = 1.e-6

    print( "{:40}{:15.6}".format('Chain velocities at temperature',temperature) )

    n, d = r.shape
    assert n==nn, "r shape mismatch {:d}{:d}".format(n,nn)
    assert d==3,  "r shape mismatch {:d}".format(d)
   
    # Confirm centre of mass is at origin
    r_cm = np.sum ( r, axis=0 ) / n
    assert np.all(r_cm<tol), "{}{:15.8f}{:15.8f}{:15.8f}".format('Centre of mass error',*r_cm)

    v = np.random.randn( n,3 )*np.sqrt(temperature) # Choose 3N random velocities

    # Compute and remove total momentum
    v_cm = np.sum ( v, axis=0 ) / n # Compute centre of mass velocity
    v = v - v_cm                    # Set net momentum to zero

    # Compute total angular momentum and moment of inertia tensor
    ang_mom = np.sum ( np.cross ( r, v ), axis=0 )
    inertia = np.zeros ( (3,3), dtype=np.float_ )
    for i in range(n):
        inertia = inertia - np.outer ( r[i,:], r[i,:] )
        for xyz in range(3):
            inertia[xyz,xyz] = inertia[xyz,xyz] + np.dot ( r[i,:], r[i,:] )

    # Solve linear system to get angular velocity
    ang_vel = la.solve(inertia,ang_mom)

    # Remove angular momentum
    v = v - np.cross ( ang_vel, r )

    if constraints:
        
        #  Apply bond constraints (which should not introduce linear or angular momentum)
        print('Applying velocity constraints relative to bonds')
        v = rattle_b ( r, v )

        # Scale velocities to get correct temperature
        # Number of degrees of freedom is 3*n - (n-1) bonds - 6 for angular and linear momentum
        temp = np.sum(v**2) / ( 3*n - (n-1) - 6 )
        v = v * np.sqrt ( temperature / temp )

    else:
        
        # Scale velocities to get correct temperature
        # Number of degrees of freedom is 3*n - 6 for angular and linear momentum
        temp = np.sum(v**2) / ( 3*n - 6 )
        v = v * np.sqrt ( temperature / temp )

    # Final check on angular and linear momenta
    v_cm = np.sum ( v, axis=0 )
    ang_mom = np.sum ( np.cross(r,v), axis=0 )
    assert not np.any(v_cm>tol), "{}{:15.8f}{:15.8f}{:15.8f}".format('Linear momentum error', *v_cm)
    assert not np.any(ang_mom>tol), "{}{:15.8f}{:15.8f}{:15.8f}".format('Angular momentum error', *ang_mom)

    return v

def rattle_b ( r, v ):
    """A version of velocity Verlet constraint algorithm."""

    import numpy as np
    
    # This subroutine iteratively adjusts the velocities stored in the array v
    # to satisfy the time derivatives of the bond constraints
    
    n, d = r.shape
    assert d==3, 'r dimension error in rattle_b'

    tol = 1.0e-9
    iter_max = 500

    iter  = 0
    done  = False
    moved = np.full(n,True,dtype=np.bool_) # Ensures that we look at each bond at least once
    move  = np.empty_like(moved)

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

                # In the following formulae, inverse masses are all unity
                g  = -0.5*np.dot ( rij, vij ) / np.dot ( rij, rij )

                if abs(g) > tol: # Test whether constraint already satisfied

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

    return v

def overlap ( ri, ei, r, e, box, ell ):
    """This routine checks for overlaps of atoms (ell=0) or spherocylinders (ell>0)."""

    import numpy as np

    tol = 1.0e-6
    
    nj, d = r.shape
    assert d==3, 'r dimension error in overlap'
    assert nj==e.shape[0], 'e dimension error in overlap'

    if nj<1:
        return False

    if ell<tol: # Handle spherical case separately (atoms or nonlinear)
        rij = ri - r
        rij = rij - np.rint(rij/box)*box
        rij_sq = np.sum ( rij**2, axis=1 )
        return np.any ( rij_sq < 1.0 )

    # Otherwise handle the nonspherical case
    ell2 = ell/2.0
        
    for j,rj in enumerate(r):
        rij = ri - rj
        rij = rij - np.rint(rij/box)*box
        rij_sq = sum(rij**2)
        rei = np.dot(rij,ei)
        rej = np.dot(rij,e[j,:])
        eij = np.dot(ei,e[j,:])

        sin_sq = 1.0 - eij**2 # Squared sine of angle between line segments

        if sin_sq < tol: # Guard against nearly-parallel lines
            ci = -rei
            cj =  rej
        else:
            ci = ( - rei + eij * rej ) / sin_sq
            cj = (   rej - eij * rei ) / sin_sq

        ai = np.fabs ( ci )
        aj = np.fabs ( cj )
        if ai > ell2:
            ci = ell2*np.sign(ci)
        if aj > ell2:
            cj = ell2*np.sign(cj)

        if ai > aj:
            cj =  rej + ci * eij
        else:
            ci = -rei + cj * eij

        ai = np.fabs ( ci )
        aj = np.fabs ( cj )
        if ai > ell2:
            ci = ell2*np.sign(ci)
        if aj > ell2:
            cj = ell2*np.sign(cj)

        di =  2.0 * rei + ci - cj * eij
        dj = -2.0 * rej + cj - ci * eij

        sij_sq = rij_sq + ci * di + cj * dj # Squared distance between line segments

        if sij_sq < 1.0:
            return True

    return False

def chain_overlap ( ri, r ):
    """This routine checks for overlaps of atoms."""

    import numpy as np

    # NO box, NO periodic boundary conditions
    
    nj, d = r.shape
    assert d==3, 'r dimension error in chain_overlap'

    if nj<1:
        return False

    rij = ri - r
    rij_sq = np.sum ( rij**2, axis=1 )
    return np.any ( rij_sq < 1.0 )

"""Sets up initial configuration for MD or MC."""

import json
import sys
import numpy as np
import math
from config_io_module import write_cnf_atoms, write_cnf_mols

filename = 'cnf.inp'
atom, linear, nonlinear, chain = 0, 1, 2, 3 # User options
tol = 1.e-6

print('initialize')
print('Sets up initial configuration file for various simulations')
print('Options for molecules are "atom", "linear", "nonlinear", "chain"')
print('Particle mass m=1 throughout')
print('Periodic boundaries')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"n":0, "nc":4, "temperature":1.0, "inertia":1.0, "density":0.75, "length":0.0, "constraints":True,
                "bond":1.122462, "velocities":False, "molecules":"atom", "lattice":True, "soft":False}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
n           = nml["n"]           if "n"           in nml else defaults["n"]
nc          = nml["nc"]          if "nc"          in nml else defaults["nc"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
inertia     = nml["inertia"]     if "inertia"     in nml else defaults["inertia"]
density     = nml["density"]     if "density"     in nml else defaults["density"]
length      = nml["length"]      if "length"      in nml else defaults["length"]
bond        = nml["bond"]        if "bond"        in nml else defaults["bond"]
velocities  = nml["velocities"]  if "velocities"  in nml else defaults["velocities"]
molecules   = nml["molecules"]   if "molecules"   in nml else defaults["molecules"]
lattice     = nml["lattice"]     if "lattice"     in nml else defaults["lattice"]
soft        = nml["soft"]        if "soft"        in nml else defaults["soft"]
constraints = nml["constraints"] if "constraints" in nml else defaults["constraints"]

np.random.seed()

molecules = molecules.lower()
assert  ( "atom" in molecules or "linear" in molecules or
          "nonlin" in molecules or "chain" in molecules    ), 'Unrecognized molecules option'

if "nonlin" in molecules:
    molecule_option = nonlinear
    print('Nonlinear molecules')
elif "linear" in molecules:
    molecule_option = linear
    print('Linear molecules')
elif "atom" in molecules:
    molecule_option = atom
    print('Atoms')
else:
    molecule_option = chain
    print('Atoms in a chain')

if n<= 0: # This is the default
    assert nc>0, "{}{:d}".format('nc must be positive',nc)
    print( "{:40}{:15d}".format('nc',nc)   )
    n = 4*nc**3 # Deduce n from nc
    print( "{:40}{:15d}".format('n',n)   )
else: # n has been specified directly
    print( "{:40}{:15d}".format('n',n)   )


if velocities:
    print('Velocities option selected')

    # Inertia should be positive, even for atoms
    if inertia < tol:
        print("{}{:15.6f}".format('Warning, inertia = ', inertia))
        print('Resetting to 1 ')
        inertia = 1.0
else:
    print('No velocities option selected')

if molecule_option == nonlinear:
    quaternions=True
    print('Periodic boundary conditions')
elif molecule_option == linear:
    quaternions=False
    print('Periodic boundary conditions')
    if length<tol:
        print("{}{:15.6f}".format('Warning, length ',length))
elif molecule_option == atom:
    quaternions=False
    print('Periodic boundary conditions')
    if length>tol:
        print("{}{:15.6f}{}".format('Warning, length ',length,' resetting to zero'))
        length = 0.0
else:
    quaternions=False
    print('NO periodic boundary conditions')
    if length>tol:
        print("{}{:15.6f}{}".format('Warning, length ',length,' resetting to zero'))
        length = 0.0
    if velocities:
        if constraints:
            print('Velocities constrained relative to bonds')
        else:
           print('Velocities not constrained relative to bonds')
    
if soft:
    print('Soft option selected - no overlap checking')

if molecule_option == chain:
    
    print( "{:40}{:15.6f}".format('Bond length',bond    ) )
    r = chain_positions ( n, bond, soft )

    if velocities:
        print( "{:40}{:15.6f}".format('Temperature',temperature    ) )
        v = chain_velocities ( n, temperature, constraints, r )

else:

    # Periodic boundaries apply
    # Box length is deduced from density
    box = ( n / density ) ** ( 1.0/3.0 )
    print( "{:40}{:15.6f}".format('Density',   density) )
    print( "{:40}{:15.6f}".format('Box length',box    ) )

    if lattice:
        r, e = fcc_positions ( n, box, length, soft, quaternions )
    else:
        r, e = ran_positions ( n, box, length, soft, quaternions )

    if velocities:
        print( "{:40}{:15.6f}".format('Temperature',temperature    ) )
        if molecule_option != atom:
            print( "{:40}{:15.6f}".format('Inertia',inertia    ) )
        v, w = ran_velocities ( n, e, temperature, inertia, quaternions )

print("{}{}".format('Writing configuration to filename ',filename))
if molecule_option == atom:
    if velocities:
        write_cnf_atoms ( filename, n, box, r, v )
    else:
        write_cnf_atoms ( filename, n, box, r )
elif molecule_option == chain:
    if velocities:
        write_cnf_atoms ( filename, n, bond, r, v )
    else:
        write_cnf_atoms ( filename, n, bond, r )
else:
    if velocities:
        write_cnf_mols ( filename, n, box, r, e, v, w )
    else:
        write_cnf_mols ( filename, n, box, r, e )
