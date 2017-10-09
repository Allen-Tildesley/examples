#!/usr/bin/env python3
# md_nvt_poly_lj.py

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

"""Molecular dynamics, NVT ensemble (NVE option), polyatomic molecules."""

def kick_propagator ( t ):
    """Advances velocities and space-fixed angular momenta."""

    # t is time over which to propagate, typically dt/2
    # v and ell must be declared global to be updated in main program
    
    global v, ell
    import numpy as np
    
    v   = v   + t * f   # Linear momenta are equivalent to velocities
    ell = ell + t * tau # Space-fixed angular momenta ell and torques tau

def drift_translate ( t ):
    """Advances positions."""

    # t is time over which to propagate, typically dt
    # r must be declared global to be updated in main program

    global r
    import numpy as np

    r = r + t * v / box   # Drift step (positions in box=1 units)
    r = r - np.rint ( r ) # Periodic boundaries

def drift_rotate ( xyz, t ):
    """Advances quaternion orientations about a specified axis."""
    
    # t is time over which to propagate, typically dt or dt/2
    # e must be declared global to be updated in main program

    global e
    import numpy as np
    from maths_module import rotate_quaternion

    for i, ei in enumerate(e):
        ai     = q_to_a ( ei )                             # Rotation matrix for i
        w_hat  = ai[xyz,:]                                 # Space-fixed axis about which to rotate
        w_mag  = np.dot ( ell[i,:], w_hat ) / inertia[xyz] # Angular velocity about this axis
        e[i,:] = rotate_quaternion ( w_mag*t, w_hat, ei )  # Rotate by specified angle

def ran_velocities ( temperature, inertia ):
    """Returns random velocities and space-fixed angular momenta."""

    import numpy as np
    from maths_module import q_to_a
    
    v      = np.random.randn ( n, 3 )         # Random velocities
    v_cm   = np.average ( v, axis=0 )         # Compute centre of mass velocity
    v      = v - v_cm                         # Set net momentum to zero
    factor = np.sum ( v**2 ) / (3*n-3)        # Estimate of kinetic temperature
    factor = np.sqrt ( temperature / factor ) # Necessary correction factor
    v      = v * factor                       # Make correction

    ell     = np.random.randn ( n, 3 )                  # Random body-fixed angular momenta
    factors = np.sum ( ell**2, axis=0 ) / ( n*inertia ) # Estimate of kinetic temperatures
    factors = np.sqrt ( temperature / factors )         # Necessary correction factors
    ell     = factors * ell                             # Make corrections

    # Convert to space-fixed angular momenta
    for i, ei in enumerate(e):
        ai       = q_to_a ( ei )           # Rotation matrix for i
        ell[i,:] = np.dot ( ell[i,:], ai ) # NB: equivalent to ell_s = ai_T*ell_b, ai_T=transpose of ai

    return v, ell

def calc_variables ( ):
    """Calculates all variables of interest.
    
    They are collected and returned as a list, for use in the main program."""
    
    import numpy as np
    import math
    from averages_module import msd, VariableType
    from maths_module    import q_to_a
    
    # Preliminary calculations (n,r,v,total are taken from the calling program)
    vol   = box**3           # Volume
    rho   = n / vol          # Density
    kin_t = 0.5*np.sum(v**2) # Translational kinetic energy
    kin_r = 0.0
    for i, ei in enumerate(e):
        ai    = q_to_a ( ei )                      # Rotation matrix for i
        ell_i = np.dot ( ai, ell[i,:] )            # Get body-fixed angular momentum
        kin_r = kin_r + np.sum((ell_i**2)/inertia) # Increment kinetic energy
    kin_r = kin_r / 2                              # Rotational kinetic energy

    tmp_t = 2.0 * kin_t / (3*n-3)     # Remove three degrees of freedom for momentum conservation
    tmp_r = 2.0 * kin_r / (3*n)       # 3N degrees of rotational freedom
    eng   = kin_t + kin_r + total.pot # Total energy for simulated system

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Internal energy (shifted-force potential) per atom
    # Total translational and rotational KE plus total PE divided by N
    if nvt:
        e_sf = VariableType ( nam = 'E/N shifted force', val = 3.0*temperature + total.pot/n )
    else:
        e_sf = VariableType ( nam = 'E/N shifted force', val = eng/n )

    # Pressure (shifted-force potential)
    # Ideal gas contribution plus total virial divided by V
    if nvt:
        p_sf = VariableType ( nam = 'P shifted force', val = rho*temperature + total.vir/vol )
    else:
        p_sf = VariableType ( nam = 'P shifted force', val = rho*tmp_t + total.vir/vol )

    # Kinetic translational temperature
    t_t = VariableType ( nam = 'T translational', val = tmp_t )

    # Kinetic rotational temperature
    t_r = VariableType ( nam = 'T rotational', val = tmp_r )

    # Mean-squared deviation of conserved energy per atom
    conserved_msd = VariableType ( nam = 'Conserved MSD', val = eng/n,
                                   method = msd, e_format = True, instant = False )

    # Collect together into a list for averaging
    return [ e_sf, p_sf, t_t, t_r, conserved_msd ]

# Takes in a configuration of atoms (positions, quaternions, velocities and angular momenta)
# Cubic periodic boundary conditions
# Conducts molecular dynamics with optional velocity thermalization
# Uses no special neighbour lists
# The rotational algorithm is described in the text, section 3.3. See A Dullweber, B Leimkuhler, 
# R McLachlan, J Chem Phys 107, 5840 (1997) and TF Miller, M Eleftheriou, P Pattnaik, 
# A Ndirango, D Newns, GJ Martyna, J Chem Phys 116, 8649 (2002).

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in and we assume mass=1 throughout
# However, input configuration, output configuration, most calculations, and all results 
# are given in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1, epsilon = 1

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model is defined in md_poly_lj_module

import json
import sys
import numpy as np
import math
from config_io_module  import read_cnf_mols, write_cnf_mols
from averages_module   import run_begin, run_end, blk_begin, blk_end, blk_add
from maths_module      import q_to_a
from md_poly_lj_module import introduction, conclusion, force, na, db, inertia, PotentialType

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'

print('md_nvt_poly_lj')
print('Molecular dynamics, constant-NVT/NVE ensemble')
print('Molecular mass=1 throughout')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":5000, "dt":0.003, "temperature":1.0, "t_interval":0 }
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
dt          = nml["dt"]          if "dt"          in nml else defaults["dt"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
t_interval  = nml["t_interval"]  if "t_interval"  in nml else defaults["t_interval"]

introduction()
np.random.seed()

# Write out parameters
nvt = t_interval > 0 and t_interval < nstep
print( "{:40}{:15d}  ".format('Number of blocks',          nblock) )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)  )
print( "{:40}{:15.6f}".format('Time step',                 dt)     )
if nvt:
    print( "{:40}        ".format('NVT ensemble'                         )  )
    print( "{:40}{:15d}  ".format('Thermalization interval', t_interval  )  )
    print( "{:40}{:15.6f}".format('Temperature',             temperature )  )
else:
    print( "{:40}        ".format('NVE ensemble'                         )  )
    t_interval = nstep+1

# Read in initial configuration
if nvt:
    n, box, r, e = read_cnf_mols ( cnf_prefix+inp_tag, quaternions=True )
    v, ell = ran_velocities ( temperature, inertia )
else:
    n, box, r, e, v, ell = read_cnf_mols ( cnf_prefix+inp_tag, with_v=True, quaternions=True )
    vcm = np.sum ( v, axis=0 ) / n # Centre-of mass velocity
    v   = v - vcm                  # Set COM velocity to zero

print( "{:40}{:15d}  ".format('Number of particles',          n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
print( "{:40}{:15.6f}".format('Density', n/box**3)  )
r = r / box           # Convert positions to box units
r = r - np.rint ( r ) # Periodic boundaries

# Calculate all bond vectors
d    = np.empty ( (n,na,3), dtype=np.float_ )
norm = np.sqrt ( np.sum(e**2,axis=1) ) # Quaternion norms
e    = e / norm[:,np.newaxis]          # Ensure normalized quaternions
for i, ei in enumerate(e):
    ai       = q_to_a ( ei )     # Rotation matrix for i
    d[i,:,:] = np.dot ( db, ai ) # NB: equivalent to ai_T*db, ai_T=transpose of ai

# Initial forces, potential, etc plus overlap check
total, f, tau = force ( box, r, d )
assert not total.ovr, 'Overlap in initial configuration'

# Initialize arrays for averaging and write column headings
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        if nvt and stp%t_interval == 0:
            v, ell = ran_velocities ( temperature, inertia )

        kick_propagator ( 0.5 * dt ) # Half-kick step

        drift_translate ( dt ) # Drift step for positions

        # Succession of drift steps for rotation about body-fixed axes
        # Depending on the values of the moments of inertia, a different nested
        # sequence of axes may produce better or worse energy conservation
        drift_rotate ( 0, 0.5*dt )
        drift_rotate ( 1, 0.5*dt )
        drift_rotate ( 2,     dt )
        drift_rotate ( 1, 0.5*dt )
        drift_rotate ( 0, 0.5*dt )

        # Calculate all bond vectors
        norm = np.sqrt ( np.sum(e**2,axis=1) ) # Quaternion norms
        e    = e / norm[:,np.newaxis]          # Ensure normalized quaternions
        for i, ei in enumerate(e):
            ai       = q_to_a ( ei )     # Rotation matrix for i
            d[i,:,:] = np.dot ( db, ai ) # NB: equivalent to ai_T*db, ai_T=transpose of ai

        total, f, tau = force ( box, r, d ) # Force evaluation
        assert not total.ovr, 'Overlap in configuration'

        kick_propagator ( 0.5 * dt ) # Half-kick step

        blk_add ( calc_variables() )

    blk_end(blk)                                                    # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'              # Number configuration by block
    write_cnf_mols ( cnf_prefix+sav_tag, n, box, r*box, e, v, ell ) # Save configuration

run_end ( calc_variables() )

total, f, tau = force ( box, r, d ) # Force evaluation
assert not total.ovr, 'Overlap in final configuration'

write_cnf_mols ( cnf_prefix+out_tag, n, box, r*box, e, v, ell ) # Save configuration
conclusion()
