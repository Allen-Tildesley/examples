#!/usr/bin/env python3
# initialize_chain.py

#------------------------------------------------------------------------------------------------#
# This software was written in 2016/17                                                           #
# by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        #
# and Dominic J. Tildesley <dominic.tildesley@epfl.ch> ("the authors"),                          #
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

def initialize_random ( n, bond, soft ):
    """Chooses chain positions randomly, at desired bond length, avoiding overlap."""
    
    import numpy as np
    from maths_module import random_vector

    tol = 1.0e-9
    iter_max = 500

    print('Chain, bonds randomly oriented, avoiding overlaps')

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
            if not overlap ( r[i,:], r[:i-1,:] ):
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

def initialize_velocities ( nn, temperature, r ):
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
            inertia[xyz,xyz] = inertia[xyz,xyz] - np.dot ( r[i,:], r[i,:] )

    # Solve linear system to get angular velocity
    ang_vel = la.solve(inertia,ang_mom)
    
    # Remove angular momentum
    v = v - np.cross ( ang_vel, r )

    #  Apply bond constraints (which should not introduce linear or angular momentum)
    v = rattle_b ( r, v )

    # Scale velocities to get correct temperature
    # Number of degrees of freedom is 3*n - (n-1) bonds - 6 for angular and linear momentum
    temp = np.sum(v**2) / ( 3*n - (n-1) - 6 )
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

def overlap ( ri, r ):
    """This routine checks for overlaps of atoms."""

    import numpy as np
    
    nj, d = r.shape
    assert d==3, 'r dimension error in overlap'

    if nj<1:
        return False

    rij = ri - r
    rij_sq = np.sum ( rij**2, axis=1 )
    return np.any ( rij_sq < 1.0 )

"""Sets up initial configuration for MD or MC."""

import json
import sys
import numpy as np
from config_io_module import write_cnf_atoms

filename = 'cnf.inp'

print('initialize_chain')
print('Sets up initial configuration file for atomic chain')
print('Particle mass m=1 throughout')
print('NO periodic boundaries')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"n":13, "temperature":1.0, "bond":1.0, "velocities":False, "soft":False}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
n           = nml["n"]           if "n"           in nml else defaults["n"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
bond        = nml["bond"]        if "bond"        in nml else defaults["bond"]
velocities  = nml["velocities"]  if "velocities"  in nml else defaults["velocities"]
soft        = nml["soft"]        if "soft"        in nml else defaults["soft"]

np.random.seed()

print( "{:40}{:15d}".format('n',n)   )

if soft:
    print('Soft option selected - no overlap checking')

r = initialize_random ( n, bond, soft )

if velocities:
    print( "{:40}{:15.6f}".format('Temperature',temperature    ) )
    v = initialize_velocities ( n, temperature, r )

print("{}{}".format('Writing configuration to filename ',filename))
if velocities:
    write_cnf_atoms ( filename, n, bond, r, v )
else:
    write_cnf_atoms ( filename, n, bond, r )
