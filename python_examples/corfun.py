#!/usr/bin/env python3
# corfun.py

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

"""Time correlation function, directly and by FFT"""

def c_exact(t):
    """Returns exact time correlation function for given random process.

    t is the time (floating point).
    """
    # See AD Baczewski and SD Bond J Chem Phys 139 044107 (2013)
    # In general the exact correlation function may be obtained from the inverse Laplace transform
    # C(s) = (kT/m) * 1 / ( s + M(s) ) where both C(t) and M(t) are sums of exponentials in t
    # M(s) = sum_p m_p*kappa_p / (s+kappa_p) p = 1..Np (note amplitude mp*kappa_p)
    # C(s) = sum_p c_p / (s+k_p) p = 1..Np+1
    # The coefficients c_p and k_p may be determined in terms of the m_p and kappa_p
    # by solving an equation of order Np+1 and using partial fractions
    # Here we just do this for the simplest case Np=1
    # Agreement with the simulated function is only expected to within statistical error

    import math
    
    if kappa > 4.0 * m: # Real roots
        root = SQRT ( 0.25*kappa**2 - m*kappa )
        kp  = 0.5*kappa + root
        km  = 0.5*kappa - root
        cp  = -km/(kp-km)
        cm  =  kp/(kp-km)
        c   = cp*math.exp(-kp*t) + cm*math.exp(-km*t)

    else: # Complex roots
        omega = math.sqrt ( m*kappa - 0.25*kappa**2 )
        c     = math.exp(-0.5*kappa*t) * ( math.cos(omega*t) + (0.5*kappa/omega)*math.sin(omega*t) )

    return c

import json
import sys
import numpy as np
import math
import time

print('corfun')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
# nstep can be made larger, but the default direct method of analysis is very slow (see below)
defaults = {"nt":1000, "origin_interval":1, "nstep":2**15, "nequil":10000, "delta":0.01, "temperature":1.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))
    
# Set parameters to input values or defaults
nt              = nml["nt"]              if "nt"              in nml else defaults["nt"]              # Max time for correlation function
origin_interval = nml["origin_interval"] if "origin_interval" in nml else defaults["origin_interval"] # Interval between origins
nstep           = nml["nstep"]           if "nstep"           in nml else defaults["nstep"]           # Number of steps
nequil          = nml["nequil"]          if "nequil"          in nml else defaults["nequil"]          # Number of equilibration timesteps
delta           = nml["delta"]           if "delta"           in nml else defaults["delta"]           # Timestep for simulation
temperature     = nml["temperature"]     if "temperature"     in nml else defaults["temperature"]     # Temperature for simulation

# Write out parameters
n0 = nt//origin_interval + 1
print( "{:40}{:15d}  ".format('Number of steps in run',    nstep)           )
print( "{:40}{:15d}  ".format('Equilibration steps',       nequil)          )
print( "{:40}{:15d}  ".format('Max correlation time nt',   nt)              )
print( "{:40}{:15d}  ".format('Origin interval',           origin_interval) )
print( "{:40}{:15d}  ".format('Number of time origins n0', n0)              )
print( "{:40}{:15.6f}".format('Time step delta',           delta)           )
print( "{:40}{:15.6f}".format('Temperature',               temperature)     )

# The memory function model is defined here
# Values used by Baczewski and Bond in their example are
# (m,kappa)
# (1.0,1.0)  (underdamped)
# (0.5,2.0)  (critically damped)
# (0.25,4.0) (overdamped)
m     = 1.0
kappa = 1.0
print( "{:40}{:15.6f}".format('m',     m)     )
print( "{:40}{:15.6f}".format('kappa', kappa) )

# Coefficients used in algorithm
x = delta*kappa
e = math.exp(-x) # theta in B&B paper
if x > 0.0001:
   b = 1.0 - math.exp(-2.0*x)
   d = 1.0 - math.exp(-x)
else: # Taylor expansions for low x (coefficients of decreasing powers)
   b = np.polyval([-2.0/3.0,4.0/3.0,-2.0,2.0,0.0],x)
   d = np.polyval([-1.0/24.0,1.0/6.0,-1.0/2.0,1.0,0.0],x)

b      = math.sqrt ( b )
b      = b * math.sqrt ( kappa/2.0 ) # alpha in B&B paper
stddev = math.sqrt(2.0*temperature)

# Data generation
cpu_1 = time.process_time()

# Initial values
vt = 0.0
s  = 0.0
v  = np.empty(nstep,dtype=np.float_)
np.random.seed()

for t in range(-nequil, nstep): # Loop over steps including an equilibration period
    # Velocity Verlet type algorithm for vt and auxiliary variable s
    vt = vt + 0.5 * delta * s
    s  = e * s - d * m * vt + b * math.sqrt(m) * stddev * np.random.randn()
    vt = vt + 0.5 * delta * s

    if t >= 0:
        v[t] = vt # Store velocities, after equilibration

cpu_2 = time.process_time()
print( "{:40}{:15.6f}".format('CPU time to generate data', cpu_2-cpu_1 ) )

# Data analysis (direct method)
# This is very slow.
# A much faster direct method uses the NumPy correlate function (see below)

c  = np.zeros(nt+1,dtype=np.float_)
n  = np.zeros(nt+1,dtype=np.float_)
t0 = np.empty(n0,dtype=np.int_)
v0 = np.empty(n0,dtype=np.float_)
mk   = -1 # Storage location of time origin
full = False

# Main loop correlating data
for t in range(nstep):
    if t%origin_interval == 0:
        mk = mk + 1
        if mk >= n0:
            full = True
            mk = mk - n0 # Overwrite older values
        t0[mk] = t    # Store time origins
        v0[mk] = v[t] # Store velocity at time origins

    # Correlate with all time origins stored so far
    nk = n0 if full else mk+1
    for k in range(nk): # Loop over time origins
        dt = t - t0[k]
        if 0 <= dt <= nt : # Check that dt is in range
            c[dt] = c[dt] + v[t] * v0[k] # Increment correlation function
            n[dt] = n[dt] + 1.0          # Increment normalizing factor

assert np.all(n>0.5), 'Normalization array error' # Should never happen
c = c / n # Normalise by number of increments

cpu_3 = time.process_time()
print( "{:40}{:15.6f}".format('CPU time for direct method', cpu_3-cpu_2 ) )

# Data analysis (FFT method)
fft_len = 2*nstep # Actual length of FFT data

# Prepare data for FFT
fft_inp = np.zeros(fft_len,dtype=np.complex_) # Fill input array with zeros
fft_inp[0:nstep] = v                          # Put data into first part (real only)

fft_out = np.fft.fft(fft_inp) # Forward FFT

fft_out = fft_out * np.conj ( fft_out ) # Square modulus

fft_inp = np.fft.ifft(fft_out) # Backward FFT (the factor of 1/fft_len is built in)

# Normalization factors associated with number of time origins
n = np.linspace(nstep,nstep-nt,nt+1,dtype=np.float_)
assert np.all(n>0.5), 'Normalization array error' # Should never happen
c_fft = fft_inp[0:nt+1].real / n

cpu_4 = time.process_time()
print( "{:40}{:15.6f}".format('CPU time for FFT method', cpu_4-cpu_3 ) )

# Data analysis (much faster direct method using NumPy library routine)
c_full = np.correlate(v,v,mode='full')
mid = c_full.size//2
c_lib = c_full[mid:mid+nt+1]/n
cpu_5 = time.process_time()
print( "{:40}{:15.6f}".format('CPU time for library correlation routine', cpu_5-cpu_4 ) )

with open("corfun.out","w") as f:
    for t in range(nt+1):
        print("{:10d}{:15.8f}{:15.8f}{:15.8f}{:15.8f}".format(t,c[t],c_fft[t],c_lib[t],c_exact(t*delta)), file=f)
