#!/usr/bin/env python3
# qmc_walk_sho.py

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

"""Quantum ground state wavefunction by random walk."""

# Program to calculate the ground state wavefunction
# for a particle in a harmonic potential, V=(x**2)/2,
# by solving the corresponding diffusion equation in imaginary time.
#
# In atomic units, mass=1, hbar=1, diffusion coefficient D = hbar**2/(2*m) = 1/2,
# so root-mean-squared displacement in time step ds is sqrt(2*D*ds) = sqrt(ds)
# Walkers are created and destroyed, depending on the potential energy relative to the trial energy et
# The program adjusts et to attempt to converge on a target number of walkers
# hence obtaining the exact ground state value et = 0.5
# This type of simulation is very sensitive to the initial guess for et,
# and the resulting time evolution is quite noisy:
# results are output as averages over the production period.
# The simulated ground state wave function may be compared with the exact result for this simple problem.
#
# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults
# The default parameters run the simulation with "et":0.5, the exact groundstate energy for this potential.
# You can then try, say, "et":0.6 and "et":0.4, observing the behaviour of the number of walkers in each case.
#
# Many of the parameters, and the updating scheme for et, are taken from the following paper:
# I Kostin, B Faber and K Schulten, Amer J Phys, 64, 633 (1996).

import json
import sys
import importlib
import numpy as np

print('qmc_walk_sho')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"n_max":2000, "n_target":500, "production_steps":20000, "equilibration_steps":20000,
             "output_interval":200, "et":0.5, "ds":0.1, "x_max":10.0, "n_bin":200}

for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ', list(defaults.keys()))
    
# Set parameters to input values or defaults
n_max               = nml["n_max"]               if "n_max"               in nml else defaults["n_max"]
n_target            = nml["n_target"]            if "n_target"            in nml else defaults["n_target"]
production_steps    = nml["production_steps"]    if "production_steps"    in nml else defaults["production_steps"]
equilibration_steps = nml["equilibration_steps"] if "equilibration_steps" in nml else defaults["equilibration_steps"]
output_interval     = nml["output_interval"]     if "output_interval"     in nml else defaults["output_interval"]
et                  = nml["et"]                  if "et"                  in nml else defaults["et"]
ds                  = nml["ds"]                  if "ds"                  in nml else defaults["ds"]
x_max               = nml["x_max"]               if "x_max"               in nml else defaults["x_max"]
n_bin               = nml["n_bin"]               if "n_bin"               in nml else defaults["n_bin"]

# Write out parameters
print( "{:40}{:15d}  ".format('Max number of walkers',                   n_max)               )
print( "{:40}{:15d}  ".format('Target number of walkers',                n_target)            )
print( "{:40}{:15d}  ".format('Number of steps for production',          production_steps)    )
print( "{:40}{:15d}  ".format('Number of steps for equilibration',       equilibration_steps) )
print( "{:40}{:15d}  ".format('Output interval',                         output_interval)     )
print( "{:40}{:15.6f}".format('Initial estimate of ground state energy', et)                  )
print( "{:40}{:15.6f}".format('Time step',                               ds)                  )
print( "{:40}{:15.6f}".format('Max x for wavefunction',                  x_max)               )
print( "{:40}{:15d}  ".format('Number of bins for wavefunction',         n_bin)               )

assert n_target <= n_max, 'n_target exceeds n_max at start'

x       = np.empty(n_max,dtype=np.float_) # position of each walker
v       = np.empty(n_max,dtype=np.float_) # potential energy of each walker
k       = np.empty(n_max,dtype=np.float_) # floating point number of replicas to make
replica = np.empty(n_max,dtype=np.int_)   # integer number of replicas to make
alive   = np.empty(n_max,dtype=np.bool_)  # flag those walkers still alive
bin     = np.zeros(n_bin,dtype=np.int_)   # histogram bins for wavefunction

np.random.seed()

# Set up an initial delta distribution of walkers at origin
n = n_target
x[0:n] = 0.0

# Set up energy averages
step_count = 0
et_avg     = 0.0
pot_avg    = 0.0

print("{:>15}{:>15}{:>15}{:>15}".format('step','nwalkers','e0','<v>'))

for step in range(-equilibration_steps,production_steps):
    x[0:n] = x[0:n] + np.sqrt(ds)*np.random.randn(n)          # Move each walker forward by ds in imaginary time (Brownian dynamics)
    v[0:n] = 0.5 * x[0:n]**2                                  # Calculate potential energy of each walker
    k[0:n] = np.exp(-ds*(v[0:n]-et)) + np.random.rand(n)      # New number of replicas
    replica[0:n] = np.minimum(np.floor(k[0:n]).astype(int),3) # Integerize and guard against too rapid growth

    # Remove walkers with replica = 0
    alive[0:n]       = replica[0:n] > 0
    nlive            = np.count_nonzero( alive[0:n] )
    x[0:nlive]       = np.extract ( alive[0:n], x[0:n])
    v[0:nlive]       = np.extract ( alive[0:n], v[0:n] )
    replica[0:nlive] = np.extract ( alive[0:n], replica[0:n] )

    # Replicate walkers with replica > 1
    n  = nlive
    for i in range(nlive):
        n_add = replica[i] - 1
        if n_add > 0:
           if n + n_add > n_max:
               print("{}{:5d}{:5d}{:5d}".format('n+n_add exceeds n_max', n, n_add, n_max) )
               sys.exit()
           x[n:n+n_add] = x[i]
           v[n:n+n_add] = v[i]
           n = n + n_add

    pot = np.sum(v[0:n])/n # potential averaged over walkers

    # Stop if n is too large or too small
    if n < 3:
        print('n is too small, increase the value of et')
    elif n > n_max-3:
        print('n is too large, decrease the value of et')

    # Periodically write the current number of walkers, trial energy, average potential
    if step%output_interval== 0:
        print("{:15d}{:15d}{:15.6f}{:15.6f}".format(step, n, et, pot) )

    # Calculate the distribution of walkers on the line and energy averages
    if step >= 0:
        hist, edges = np.histogram ( x[0:n], bins=n_bin, range=(-x_max,x_max) )
        bin = bin + hist
        et_avg     = et_avg     + et
        pot_avg    = pot_avg    + pot
        step_count = step_count + 1

    # Reset trial energy following Kozstin et al (1996)
    et = pot + (1.0-(n/n_target))/ds

# Normalize the wavefunction so that the integral over psi**2 = 1.0
dx     = edges[1]-edges[0]
psi    = bin.astype(np.float_)   # un-normalized wave function
factor = np.sum ( psi**2 ) * dx  # integral, assuming that -x_max .. +x_max catches everything
psi    = psi / np.sqrt( factor ) # normalizing factor

# Print averages
et_avg  = et_avg / step_count
pot_avg = pot_avg / step_count
print( "{:40}{:15.6f}".format('Average trial energy', et_avg ) )
print( "{:40}{:15.6f}".format('Average potential energy', pot_avg ) )

# Print ground state wave function
print("{:>15}{:>15}{:>15}".format('x','psi(x)','psi_exact(x)'))
const = (1.0/np.pi) ** (1.0/4.0)
for i in range(n_bin):
    x_bin  = 0.5*(edges[i]+edges[i+1])
    pexact = const * np.exp(- 0.5 * x_bin * x_bin)
    print("{:15.6f}{:15.6f}{:15.6f}".format(x_bin, psi[i], pexact) )
