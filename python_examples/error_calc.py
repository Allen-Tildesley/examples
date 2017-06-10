#!/usr/bin/env python3
# error_calc.py

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

"""Estimated error in correlated data."""

# In this program we analyse synthetic data, a time series whose properties are exactly known
# Define underlying process by generalized Langevin equation (GLE)
# with memory function expressed as a decaying exponential.
# See G Ciccotti and JP Ryckaert Mol Phys 40 141 (1980)
# and AD Baczewski and SD Bond J Chem Phys 139 044107 (2013).

import json
import sys
import numpy as np
import math

print('error_calc')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"nstep":2**16, "nequil":10000, "nrepeat":50, "delta":0.01, "variance":1.0, "average":1.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))
    
# Set parameters to input values or defaults
nstep    = nml["nstep"]    if "nstep"    in nml else defaults["nstep"]    # Number of steps
nequil   = nml["nequil"]   if "nequil"   in nml else defaults["nequil"]   # Number of equilibration timesteps
nrepeat  = nml["nrepeat"]  if "nrepeat"  in nml else defaults["nrepeat"]  # Number of simulation repeats for brute force empirical calculation
delta    = nml["delta"]    if "delta"    in nml else defaults["delta"]    # Timestep for simulation
variance = nml["variance"] if "variance" in nml else defaults["variance"] # Desired variance of data (equivalent to temperature)
average  = nml["average"]  if "average"  in nml else defaults["average"]  # Desired average value of data

# Write out parameters
print( "{:40}{:15d}  ".format('Number of steps in run', nstep)    )
print( "{:40}{:15d}  ".format('Equilibration steps',    nequil)   )
print( "{:40}{:15d}  ".format('Number of repeats',      nrepeat)  )
print( "{:40}{:15.6f}".format('Time step delta',        delta)    )
print( "{:40}{:15.6f}".format('Desired average value',  average)  )
print( "{:40}{:15.6f}".format('Desired variance',       variance) )

# The memory function model is defined here
# Values used by Baczewski and Bond in their example are
# (m,kappa)
# (1.0,1.0)  (underdamped)
# (0.5,2.0)  (critically damped)
# (0.25,4.0) (overdamped)
m     = 0.25
kappa = 4.0
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
stddev = math.sqrt(2.0*variance)     # NB stddev of random forces, not data

# For this process, the results of interest can be calculated exactly
# The time correlation function is known, and hence the statistical inefficiency (SI)
# From this, the error in the mean for a run of any length can be calculated
tcor = 1.0 / m      # Correlation time is determined by memory function
tcor = tcor / delta # Express this in timesteps
si   = 2*tcor       # Statistical inefficiency (SI)
err  = np.sqrt(si*variance/nstep)
print( "{:40}{:15.6f}".format('Exact correlation time in steps', tcor      ) )
print( "{:40}{:15.6f}".format('Run length / correlation time',   nstep/tcor) )
print( "{:40}{:15.6f}".format('Exact value of SI',               si        ) )
print( "{:40}{:15.6f}".format('Exact error estimate',            err       ) )

# Data generation
np.random.seed()
a = np.empty(nstep,dtype=np.float_)
a_avg = 0.0 # Zero average accumulator
a_var = 0.0 # Zero mean-squared accumulator

for irepeat in range(nrepeat): # Loop over repeats of simulation
    
    # Initial values, hopefully not critical
    at = 0.0
    s  = 0.0

    for t in range(-nequil, nstep): # Loop over steps including an equilibration period
        # Velocity Verlet type algorithm for at and auxiliary variable s
        at = at + 0.5 * delta * s
        s  = e * s - d * m * at + b * math.sqrt(m) * stddev * np.random.randn()
        at = at + 0.5 * delta * s

        if t >= 0:
            a[t] = average + at # Store values (after equilibration, adding average)

    a_run = np.sum(a) / nstep # The run average
    a_avg = a_avg + a_run     # Average over runs
    a_var = a_var + a_run**2  # Mean squared value over runs

a_avg = a_avg / nrepeat             # Mean value
a_var = a_var / nrepeat             # Mean-squared value
a_var = a_var - a_avg**2            # Mean-squared deviation
a_var = a_var * nrepeat/(nrepeat-1) # Bias correction
err   = np.sqrt(a_var)              # Empirical standard deviation in run average
print( "{:40}{:15.6f}".format('Empirical error estimate', err ) )
print( 'This should be in reasonable agreement with exact error estimate' )

# Now analyse the last run, as if it were the only one that had been carried out
# This is what usually happens; we rarely have the luxury of many independent runs

# Simple calculation of average and variance
a_avg   = np.sum(a)/nstep        # Sample average
a       = a - a_avg              # Centre the data
a_var_1 = np.sum(a**2)/(nstep-1) # Bias-corrected sample variance
a_err   = np.sqrt(a_var_1/nstep) # Error estimate neglecting any correlations
print( "{:40}{:15.6f}".format('Sample average value', a_avg ) )
print( "{:40}{:15.6f}".format('Deviation from exact average', a_avg-average ) )
print( 'Deviation should (typically) lie within +/- exact error estimate' )
print( "{:40}{:15.6f}".format('Sample variance', a_var_1 ) )
print( "{:40}{:15.6f}".format('Error estimate neglecting correlations', a_err ) )
print( 'This should be very over-optimistic!' )

# We must take account of the correlations between successive values in time
# The two common methods which follow are actually very similar
# They differ in the choice of block lengths

# Traditional block analysis
# The rationale here is that 20 (say) independent block averages should be enough to get a reasonable
# estimate of the desired quantities, and that there is not a lot to gain by looking at more (shorter) blocks.
# Some workers just assume that the run may be divided into 10 or 20 blocks, which they hope will be independent.
# This is exactly what we do in our example programs, just for simplicity.
# We cannot recommend this in general, unless there is good reason to support the assumption of independence.
# If the 20 blocks are not independent, then attention should be focused on fewer (longer) blocks,
# rather than more (shorter) ones, and a plot of squared error estimate, or statistical inefficiency,
# vs 1/tblock carried out to extrapolate to tblock=nstep. The loop below provides the data for that plot.

print( "{:>15}{:>15}{:>15}{:>15}".format('tblock', 'nblock', 'error estimate', 'estimate of SI') )
for nblock in range(20,3,-1):
    tblock = nstep//nblock              # Block length in steps (rounded down)
    trun   = nblock*tblock              # Run length in steps, accounting for rounding
    a_run  = np.sum(a[0:trun+1]) / trun # Average of data
    a_var  = 0.0                        # Zero mean-square block average accumulator
    for stp in range(0,trun,tblock): # Loop over block starting steps
        a_blk = np.sum(a[stp:stp+tblock]-a_run)/tblock # Block average of deviation
        a_var = a_var + a_blk**2                       # Mean-square block average
    a_var = a_var / (nblock-1)       # Bias-corrected variance of block averages
    a_err = np.sqrt(a_var/nblock)    # Estimate of error from block-average variance
    si    = tblock * a_var / a_var_1 # Statistical inefficiency
    print( "{:15d}{:15d}{:15.6f}{:15.6f}".format(tblock, nblock, a_err, si) )

print('Plateau at large tblock (small nblock)')
print('should agree quite well with exact error estimate')
print('Can plot SI or error**2 against 1/tblock')

# Method of Flyvbjerg and Petersen, J Chem Phys, 91, 461 (1989)
# Note that, in this implementation, we over-write the data array a
# Advantages of this method are the very neat method for reducing the data
# and the formal scaling analysis that accompanies the blocking transformation
# The basic calculation is the same as for the traditional block-averaging,
# but the block length, and hence tblock, change by a factor 2 each time.
# Advocates of Flyvbjerg and Petersen might argue that the additional data points
# at low nblock (high tblock) which are calculated in the traditional method
# do not actually provide much new (independent) information.
# One may attempt an extrapolation of si or a_err**2 as a function of 1/tblock
# but F&P suggest estimating a plateau in a plot vs number of blocking transformations
# which is essentially log2(tblock)

nblock = nstep
tblock = 1
print( "{:>15}{:>15}{:>15}{:>15}".format('tblock', 'nblock', 'error estimate', 'estimate of SI') )

while True: # Loop over number, and hence length, of blocks
    nblock = nblock // 2 # Halve the number of blocks, rounding down if nblock is odd
    tblock = tblock*2    # Double the block length
    if nblock < 3:
        break
    a[0:nblock] = ( a[0:2*nblock-1:2] + a[1:2*nblock:2] ) / 2.0 # Blocking transformation, halving the data set
    a_avg       = np.sum ( a[0:nblock] ) / nblock               # Re-compute sample average
    a[0:nblock] = a[0:nblock] - a_avg                           # Re-centre in case of dropped points (odd nblock)
    a_var       = np.sum ( a[0:nblock]**2 ) / (nblock-1)        # Bias-corrected variance of block averages
    a_err       = np.sqrt ( a_var / nblock )                    # Estimate of error from block average variance
    si          = tblock * a_var / a_var_1                      # Statistical inefficiency
    print( "{:15d}{:15d}{:15.6f}{:15.6f}".format(tblock, nblock, a_err, si) )

print('Plateau at large tblock (small nblock)')
print('should agree quite well with exact error estimate')
print('Can plot SI or error**2 against 1/tblock or log2(tblock)')

