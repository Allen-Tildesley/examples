#!/usr/bin/env python3
# md_npt_lj.py

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

"""Molecular dynamics, NPT ensemble."""

def calc_variables ( ):
    """Calculates all variables of interest.
    
    They are collected and returned as a list, for use in the main program.
    """

    from averages_module import msd, VariableType
    from lrc_module import potential_lrc, pressure_lrc
    import numpy as np
    import math

    # Preliminary calculations (n,r,v,f,total are taken from the calling program)
    vol = box**3                  # Volume
    rho = n / vol                 # Density
    kin = 0.5*np.sum(v**2)        # Kinetic energy
    fsq = np.sum ( f**2 )         # Total squared force
    ext = np.sum(0.5*p_eta**2/q) + np.sum(0.5*p_eta_baro**2/q_baro) + 0.5*p_eps**2/w_eps + temperature * (
          g*eta[0] + np.sum(eta[1:])+ np.sum(eta_baro) )  # Extra terms for conserved variable
    eng = kin + total.pot         # Total energy (cut-and-shifted)

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Internal energy (cut-and-shifted) per atom
    # Total KE plus total cut-and-shifted PE divided by N
    e_s = VariableType ( nam = 'E/N cut&shifted', val = eng/n )

    # Internal energy (full, including LRC) per atom
    # LRC plus total KE plus total cut (but not shifted) PE divided by N
    e_f = VariableType ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + (kin+total.cut)/n )

    # Pressure (cut-and-shifted)
    # Ideal gas contribution plus total virial divided by V
    p_s = VariableType ( nam = 'P cut&shifted', val = rho*temperature + total.vir/vol )

    # Pressure (full, including LRC)
    # LRC plus ideal gas contribution plus total virial divided by V
    p_f = VariableType ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*temperature + total.vir/vol )

    # Kinetic temperature
    t_k = VariableType ( nam = 'T kinetic', val = 2.0*kin/g )

    # Configurational temperature
    # Total squared force divided by total Laplacian
    t_c = VariableType ( nam = 'T config', val = fsq/total.lap )

    # Density
    density = VariableType ( nam = 'Density', val = rho )

    # Conserved energy-like quantity per atom
    # Energy plus PV term plus extra terms defined above
    enp = eng + pressure*vol
    conserved = VariableType ( nam = 'Conserved/N', val = (enp+ext)/n )

    # Volume MSD
    vol_msd = VariableType ( nam = 'Volume MSD', val = vol, method = msd, instant = False )

    # Mean-squared deviation of conserved energy-like quantity
    # Energy plus PV term plus extra terms defined above
    enp = eng + pressure*vol
    conserved_msd = VariableType ( nam = 'Conserved MSD', val = (enp+ext)/n,
                                   method = msd, e_format = True, instant = False )

    # Collect together into a list for averaging
    return [ e_s, p_s, e_f, p_f, t_k, t_c, density, conserved, vol_msd, conserved_msd ]

def u1_propagator ( t ):
    """U1 and U1' combined: position and strain drift step propagator.

    t is the time over which to propagate (typically dt).
    r, v, eps, p_eps, w_eps, and box are accessed from the calling program.
    """

    global r, eps, box
    import numpy as np

    # U1 part
    # The propagator for r looks different from the formula given on p142 of the text.
    # However it is easily derived from that formula, bearing in mind that coordinates in
    # this program are divided by the box length, which is itself updated in this routine.

    x = t * p_eps / w_eps # Time step * time derivative of strain
    c = (1.0-np.exp(-x))/x if x>0.001 else np.polyval([-1/24,1/6,-1/2,1.0],x) # Guard against small values

    r = r + c * t * v / box   # Positions in box=1 units
    r = r - np.rint ( r ) # Periodic boundaries

    # U1' part
    # Because we divide by box above, it is important to update eps afterwards
    # If we did not use box-scaled coordinates this would not matter

    eps = eps + x            # Update strain
    box = box0 * np.exp(eps) # Update box length

def u2_propagator ( t ):
    """U2: velocity kick step propagator.

    t is the time over which to propagate (typically dt/2).
    v, p_eps, w_eps, and g are accessed from the calling program.
    """

    global v
    import numpy as np

    alpha = 1.0 + 3 / g
    x = t * alpha * p_eps / w_eps
    c = (1.0-np.exp(-x))/x if x>0.001 else np.polyval([-1/24,1/6,-1/2,1.0],x) # Guard against small values

    v = v*np.exp(-x) + c * t * f

def u2p_propagator ( t ):
    """U2': strain momentum propagator.

    t is the time over which to propagate (typically dt/2).
    v, p_eps, total, pressure, box, and g are accessed from the calling program.
    """

    global p_eps
    import numpy as np

    alpha = 1.0 + 3 / g
    pv    = alpha * np.sum(v**2) / 3 + total.vir # total.vir is the total virial
    p_eps = p_eps + 3.0 * ( pv - pressure*box**3 ) * t

def u3_propagator ( t ):
    """U3 and U3' combined: thermostat propagator.

    t is the time over which to propagate (typically dt/2).
    v, eta, p_eta, and q are accessed from the calling program.
    """

    global v, eta, eta_baro, p_eps
    import numpy as np

    # U3 part
    
    v = v * np.exp ( -t * p_eta[0] / q[0] )
    eta = eta + t * p_eta / q

    # U3' part

    eta_baro = eta_baro + t * p_eta_baro / q_baro
    p_eps    = p_eps * np.exp ( -t * p_eta_baro[0] / q_baro[0] )

def u4_propagator ( t, j_list ):
    """U4 and U4' combined: thermostat propagator.

    t is the time over which to propagate (typically dt/4).
    j_list determines the order in which to tackle variables
    v, p_eta, g, temperature, p_eps, w_eps, q etc are accessed from the calling program.
    """

    global p_eta, p_eta_baro
    import numpy as np

    # U4 part
    
    for j in j_list:
        if j==0:
            gj = np.sum(v**2) - g*temperature # The driver Gj for p_eta[0] is different
        else:
            gj = ( p_eta[j-1]**2 / q[j-1] ) - temperature
        if j==m-1:
            p_eta[j]  = p_eta[j] + t * gj # The equation for p_eta[M-1] is different
        else:
            x = t * p_eta[j+1]/q[j+1]
            c = (1.0-np.exp(-x))/x if x>0.001 else np.polyval([-1/24,1/6,-1/2,1.0],x) # Guard against small values
            p_eta[j] = p_eta[j]*np.exp(-x) + t * gj * c

    # U4' part
    
    for j in j_list:
        if j==0:
            gj = p_eps**2/w_eps - temperature # The driver Gj for p_eta_baro[0] is different
        else:
            gj = ( p_eta_baro[j-1]**2 / q_baro[j-1] ) - temperature
        if j==m-1:
            p_eta_baro[j]  = p_eta_baro[j] + t * gj # The equation for p_eta_baro[M-1] is different
        else:
            x = t * p_eta_baro[j+1]/q_baro[j+1]
            c = (1.0-np.exp(-x))/x if x>0.001 else np.polyval([-1/24,1/6,-1/2,1.0],x) # Guard against small values
            p_eta_baro[j] = p_eta_baro[j]*np.exp(-x) + t * gj * c
  
# Takes in a configuration of atoms (positions, velocities)
# Cubic periodic boundary conditions
# Conducts molecular dynamics using a measure-preserving algorithm for NPT
# Nose-Hoover chains are used, following Martyna et al, Molec Phys, 87, 1117 (1996)
# and Tuckerman et al J Phys A, 39, 5629 (2006)
# To keep this example reasonably simple, we do not subdivide the timesteps with a
# Suzuki-Yoshida decomposition, as described in those papers
# Uses no special neighbour lists

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in and we assume mass=1 throughout
# However, input configuration, output configuration, most calculations, and all results 
# are given in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1, epsilon = 1

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model is defined in md_lj_module

import json
import sys
import numpy as np
import math
from config_io_module import read_cnf_atoms, write_cnf_atoms
from averages_module  import run_begin, run_end, blk_begin, blk_end, blk_add
from md_lj_module     import introduction, conclusion, force, PotentialType

cnf_prefix = 'cnf.'
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'
m = 3 # Number of Nose-Hoover chain variables

print('md_npt_lj')
print('Molecular dynamics, constant-NPT ensemble')
print('Particle mass=1 throughout')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nblock":10, "nstep":1000, "r_cut":2.5, "dt":0.005, "temperature":1.0,
            "pressure":0.99, "tau":2.0, "tau_baro":2.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
nblock      = nml["nblock"]      if "nblock"      in nml else defaults["nblock"]
nstep       = nml["nstep"]       if "nstep"       in nml else defaults["nstep"]
r_cut       = nml["r_cut"]       if "r_cut"       in nml else defaults["r_cut"]
dt          = nml["dt"]          if "dt"          in nml else defaults["dt"]
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
pressure    = nml["pressure"]    if "pressure"    in nml else defaults["pressure"]
tau         = nml["tau"]         if "tau"         in nml else defaults["tau"]
tau_baro    = nml["tau_baro"]    if "tau_baro"    in nml else defaults["tau_baro"]

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)       )
print( "{:40}{:15.6f}".format('Potential cutoff distance', r_cut)       )
print( "{:40}{:15.6f}".format('Time step',                 dt)          )
print( "{:40}{:15.6f}".format('Specified temperature',     temperature) )
print( "{:40}{:15.6f}".format('Specified pressure',        pressure)    )
print( "{:40}{:15.6f}".format('Thermostat timescale',      tau)         )
print( "{:40}{:15.6f}".format('Barostat timescale',        tau_baro)    )
print( "{:40}{:15d}  ".format('Nose-Hoover chain length',  m)           )

# Read in initial configuration
n, box, r, v = read_cnf_atoms ( cnf_prefix+inp_tag, with_v=True)
print( "{:40}{:15d}  ".format('Number of particles',          n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
print( "{:40}{:15.6f}".format('Density', n/box**3)  )
r = r / box                    # Convert positions to box units
r = r - np.rint ( r )          # Periodic boundaries
vcm = np.sum ( v, axis=0 ) / n # Centre-of mass velocity
v = v - vcm                    # Set COM velocity to zero

# Initial values of thermal variables
g    = 3*(n-1)
q    = np.full(m,temperature * tau**2,dtype=np.float_) 
q[0] = g * temperature * tau**2
fmt_string = "{:15.6f}"*m
fmt_string = "{:40}"+fmt_string
print(fmt_string.format('Thermal inertias Q',*q))
eta   = np.zeros(m,dtype=np.float_)
p_eta = np.random.randn(m)*np.sqrt(temperature)
p_eta = p_eta * np.sqrt(q)
q_baro = np.full(m,temperature * tau_baro**2,dtype=np.float_)
print(fmt_string.format("Barostat thermal inertias Q'",*q_baro))
eta_baro = np.zeros(m,dtype=np.float_)
p_eta_baro = np.random.randn(m)*np.sqrt(temperature)
p_eta_baro = p_eta_baro * np.sqrt(q_baro)
w_eps = g * temperature * tau_baro**2
print( "{:40}{:15.6f}".format('Barostat inertia W', w_eps)  )
box0  = box # Reference box length for strain
eps   = 0.0 # Initial strain; generally eps = log ( box/box0 )
p_eps = np.random.randn() * np.sqrt(temperature*w_eps) # strain momentum

# Initial forces, potential, etc plus overlap check
total, f = force ( box, r_cut, r )
assert not total.ovr, 'Overlap in initial configuration'

# Initialize arrays for averaging and write column headings
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        u4_propagator ( dt/4, list(reversed(range(m))) ) # Inwards order
        u3_propagator ( dt/2 )
        u4_propagator ( dt/4, list(range(m)) ) # Outwards order

        u2p_propagator ( dt/2 )
        u2_propagator  ( dt/2 )

        u1_propagator ( dt )

        total, f = force ( box, r_cut, r ) # Force evaluation
        assert not total.ovr, 'Overlap in configuration'

        u2_propagator  ( dt/2 )
        u2p_propagator ( dt/2 )

        u4_propagator ( dt/4, list(reversed(range(m))) ) # Inwards order
        u3_propagator ( dt/2 )
        u4_propagator ( dt/4, list(range(m)) ) # Outwards order

        blk_add ( calc_variables() )

    blk_end(blk)                                             # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'       # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r*box, v ) # Save configuration

run_end ( calc_variables() )

total, f = force ( box, r_cut, r ) # Force evaluation
assert not total.ovr, 'Overlap in final configuration'

write_cnf_atoms ( cnf_prefix+out_tag, n, box, r*box, v ) # Save configuration
conclusion()
