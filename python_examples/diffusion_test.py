#!/usr/bin/env python3
# diffusion_test.py

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

"""Generate test data for diffusion.py."""

def a_propagator ( t, r, v ):
    """A propagator (drift).

    t is the time over which to propagate (typically dt/2)
    r and v are the current positions and velocities
    The function returns the new positions.
    """
    return r + v * t

def o_propagator ( t, v ):
    """O propagator (friction and random contributions).

    t is the time over which to propagate (typically dt)
    v is the current velocity
    gamma and temperature are accessed from the calling program
    The function returns the new velocities.
    """
    import math
    import numpy as np
    
    x = gamma*t
    # Taylor expansion for small x
    c = 1.0-math.exp(-2.0*x) if x>0.0001 else np.polyval([-2.0/3.0,4.0/3.0,-2.0,2.0,0.0],x)
    c = math.sqrt(c)
    return np.exp(-x) * v + c * np.sqrt(temperature) * np.random.randn(n,3)

# diffusion_test program

import json
import sys
import numpy as np
import math
from config_io_module import write_cnf_atoms

print('diffusion_test')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"n":250, "nblock":999, "nstep":25, "dt":0.002, "gamma":1.0, "temperature":1.0, "box":1.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))
    
# Set parameters to input values or defaults
n           = nml["n"]           if "n"           in nml else defaults["n"]
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
dt          = nml["dt"]          if "dt"          in nml else defaults["dt"]
gamma       = nml["gamma"]       if "gamma"       in nml else defaults["gamma"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
box         = nml["box"]         if "box"         in nml else defaults["box"]

# Write out parameters
print( "{:40}{:15d}  ".format('Number of atoms',           n)                 )
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)            )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)             )
print( "{:40}{:15.6f}".format('Time step',                 dt)                )
print( "{:40}{:15.6f}".format('Friction coefficient',      gamma)             )
print( "{:40}{:15.6f}".format('Temperature',               temperature)       )
print( "{:40}{:15.6f}".format('Ideal diffusion coefft',    temperature/gamma) )

np.random.seed()

r = np.random.rand(n,3).astype(np.float_)  # Random positions
r = r - 0.5                                # Now in range (-1/2,1/2)
r = r * box                                # Now in range (-box/2,box/2)
v = np.random.randn(n,3).astype(np.float_) # Random velocities
v = v * np.sqrt(temperature)               # At desired temperature

cnf_prefix = 'cnf.'
sav_tag = '000'
write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r, v )

for blk in range(nblock): # Loop over blocks
    for stp in range(nstep): # Loop over steps
        r = a_propagator(dt/2.0,r,v) # A drift half-step
        v = o_propagator(dt,v)       # O random velocities and friction step
        r = a_propagator(dt/2.0,r,v) # A drift half-step
    r = r - np.rint(r/box)*box                           # Periodic boundaries
    sav_tag = str(blk+1).zfill(3) if blk<999 else 'sav'  # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r, v ) # Save configuration

print('Exact results output to diffusion_exact.out')
with open("diffusion_exact.out","w") as f:
    for blk in range(nblock//2): # Loop up to half the run length
        t    = (blk*nstep)*dt                                                         # Time advances block by block
        vacf = 3.0*temperature * math.exp(-gamma*t)                                   # Velocity autocorrelation function
        rvcf = 3.0*temperature * ( 1.0 - math.exp(-gamma*t) ) / gamma                 # Velocity-displacement correlation
        msd  = 6.0*temperature * ( t - ( 1.0 - math.exp(-gamma*t) ) / gamma ) / gamma # Mean-square displacement
        print("{:15.6f}{:15.8f}{:15.8f}{:15.8f}".format(t, vacf, rvcf, msd), file=f)
