#!/usr/bin/env python3
# md_nve_hs.py

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

"""Molecular dynamics, NVE ensemble, hard spheres."""

def calc_variables ( ):
    """Calculates all variables of interest.
    
    They are collected and returned as a list, for use in the main program."""
    
    import numpy as np
    from averages_module import VariableType
    
    # Preliminary calculations
    vol = box**3  # Volume
    rho = n / vol # Density

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # We average over the time step
    coll_rate  = VariableType ( nam = 'Collision rate', val = 2.0*col_sum/dt/n, instant = False )
    # ideal + collisional virial / volume averaged over the time step
    p_coll     = VariableType ( nam = 'P', val = rho*temp_kinet + vir_sum/dt/vol, instant = False )

    # Collect together into a list for averaging
    return [ coll_rate, p_coll ]

def advance ( t, box, t_now, coltime, r, v ):
    """Advances positions and reduces collision times."""

    import numpy as np

    # Guard against going back in time (should never happen)
    assert t>0.0, 'Negative time step'
    
    t_now   = t_now + t         # Advance current time by t
    coltime = coltime - t       # Reduce times to next collision by t
    r       = r + t * v / box   # Advance all positions by t (box=1 units)
    r       = r - np.rint ( r ) # Apply periodic boundaries
    return t_now, coltime, r

# Takes in a hard-sphere configuration (positions and velocities)
# Checks for overlaps    
# Conducts molecular dynamics simulation
# Uses no special neighbour lists
# ... so is restricted to small number of atoms
# Assumes that collisions can be predicted by looking at 
# nearest neighbour particles in periodic boundaries
# ... so is unsuitable for low densities

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are stored divided by the box length
# However, input configuration, output configuration, most calculations, and all results 
# are given in units sigma = 1, mass = 1

import json
import sys
import numpy as np
import math
from config_io_module import read_cnf_atoms, write_cnf_atoms
from averages_module  import run_begin, run_end, blk_begin, blk_end, blk_add
from md_nve_hs_module import introduction, conclusion, update, dndate, collide, overlap

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('md_nve_hs')
print('Molecular dynamics, constant-NVE, hard spheres')
print('Particle mass=1 throughout')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":2000, "dt":0.05}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock = nml["nblock"] if "nblock" in nml else defaults["nblock"]
nstep  = nml["nstep"]  if "nstep"  in nml else defaults["nstep"]
dt     = nml["dt"]     if "dt"     in nml else defaults["dt"]

introduction()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock) )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)  )
print( "{:40}{:15.6f}".format('Time step',                 dt)     )

# Read in initial configuration
n, box, r, v = read_cnf_atoms ( cnf_prefix+inp_tag, with_v=True)
print( "{:40}{:15d}  ".format('Number of particles',          n)         )
print( "{:40}{:15.6f}".format('Box length (sigma units)',     box)       )
print( "{:40}{:15.6f}".format('Density',                      n/box**3)  )
r = r / box                    # Convert positions to box units
r = r - np.rint ( r )          # Periodic boundaries
vcm = np.sum ( v, axis=0 ) / n # Centre-of mass velocity
v = v - vcm                    # Set COM velocity to zero
kin        = 0.5 * np.sum ( v**2 )
temp_kinet = 2.0 * kin / ( 3*(n-1) )
v          = v / np.sqrt ( temp_kinet ) # We fix the temperature to be 1.0
kin        = 0.5 * np.sum ( v**2 )
temp_kinet = 2.0 * kin / ( 3*(n-1) )
print( "{:40}{:15.6f}".format('Temperature', temp_kinet)  )
  
# Initial overlap check
assert not overlap ( box, r ), 'Particle overlap in initial configuration'

# Initial search for collision partners >i
coltime = np.full ( n, 1.0e9, dtype=np.float_ )
partner = np.full ( n, n-1, dtype=np.int_ )
for i in range(n-1):
  coltime[i], partner[i] = update ( i, box, r[i:,:], v[i:,:] )

# Initialize arrays for averaging and write column headings
col_sum = 0
vir_sum = 0.0
run_begin ( calc_variables() )

ncoll = 0

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        vir_sum = 0.0 # Zero collisional virial accumulator for this step
        col_sum = 0   # Zero collision counter for this step
        t_now   = 0.0 # Keep track of time within this step

        while True: # Loop over collisions within this step

            i   = np.argmin ( coltime ) # Locate minimum collision time
            j   = partner[i]            # Collision partner
            tij = coltime[i]            # Time to collision

            if t_now + tij > dt:
                t_now, coltime, r = advance ( dt-t_now, box, t_now, coltime, r, v ) # Advance to end of time step
                break                                                               # Exit loop over collisions

            t_now, coltime, r = advance ( tij, box, t_now, coltime, r, v ) # Advance to time of next collision

            # Compute collision dynamics
            v[i,:], v[j,:], vir = collide ( r[i,:], v[i,:], r[j,:], v[j,:], box ) 
            col_sum = col_sum + 1   # Increment collision counter
            vir_sum = vir_sum + vir # Increment collisional virial accumulator

            # Update collision lists
            for k in range(n-1):
                if k==i or k==j or partner[k]==i or partner[k]==j:
                    coltime[k], partner[k] = update ( k, box, r[k:,:], v[k:,:] )
            # Search for partners <i
            if i>0:
                coltime[:i], partner[:i] = dndate ( i, box, r[:i+1,:], v[:i+1,:], coltime[:i], partner[:i] )
            if j>0:
                coltime[:j], partner[:j] = dndate ( j, box, r[:j+1,:], v[:j+1,:], coltime[:j], partner[:j] ) 

        ncoll = ncoll + col_sum

        # Calculate and accumulate variables for this step
        blk_add ( calc_variables() )

    blk_end(blk)                                             # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'       # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r*box, v ) # Save configuration

run_end ( calc_variables() )

print( "{:40}{:15d}  ".format('Total collisions', ncoll)         )
assert not overlap ( box, r ), 'Particle overlap in final configuration'

write_cnf_atoms ( cnf_prefix+out_tag, n, box, r*box, v ) # Save configuration
conclusion()
