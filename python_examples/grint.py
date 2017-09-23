#!/usr/bin/env python3
# grint.py

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

def calculate ( r ):
    """Carry out the histogramming for rho2."""

    import numpy as np

    zcr      = np.empty_like(r)    # Create array to be histogrammed
    zcr[:,0] = r[:,2]              # All the zi values (unchanged throughout)
    r2       = np.zeros_like(rho2) # Create array to hold rho2 snapshot
    bins     = ( 2*nz+1, 2*nc+1, nr )
    ranges   = [[-z_max,z_max],[-1.0,1.0],[0.0,nr*dr]]

    for shift in range(1,n): # Loop over shifts, covering all ij and ji, omitting shift=n
        rij      = r - np.roll(r,shift,axis=0) # Difference between N pairs of coordinates
        rij      = rij - np.rint(rij/box)*box  # Apply periodic boundaries
        rij_sq   = np.sum(rij**2,axis=1)       # Squared separations
        rij_mag  = np.sqrt(rij_sq)             # Separation distances
        weights  = np.reciprocal(rij_sq)       # Incorporate 1/rij**2 factor here
        zcr[:,1] = rij[:,2] / rij_mag          # Separation vector cosines
        zcr[:,2] = rij_mag                     # Separation distances

        h, edges = np.histogramdd ( zcr, bins=bins, range=ranges, weights=weights )
        r2       = r2 + h

    return r2

def f1tanh ( z, z_gl, width, rho_g, rho_l ):
    """1-tanh fit function."""
    import numpy as np
    t = np.tanh ( ( z - z_gl ) / width )
    return 0.5 * ( rho_l + rho_g ) + 0.5 * ( rho_l - rho_g ) * t

def f2tanh ( z, z_gl, z_lg, width, rho_g, rho_l ):
    """2-tanh fit function."""
    import numpy as np
    t1 = np.tanh ( ( z - z_gl ) / width )
    t2 = np.tanh ( ( z - z_lg ) / width )
    return rho_g + 0.5 * ( rho_l - rho_g ) * ( t1 - t2 )
    
"""g(z,c,r) in a planar interface."""

import json
import sys
import numpy as np
from scipy.optimize   import curve_fit
from config_io_module import read_cnf_atoms
import os.path

# Reads a trajectory from a sequence of configuration files
# Calculates pair distribution function for a planar interface in the xy plane,
# including dependence on z and symmetry breaking with respect to z direction

# Single-particle density profile in box coordinates is written to a file 'den.out'
# The combined profile for both interfaces, relative to interface position, is written to rho.out
# Slices through the pair distribution function g2(z,c,r) where z=z1 is the z-coordinate of atom 1
# and c=cos(theta) is the angle of the r12 vector, are written out to files as a function of r.
# Assuming that the liquid phase is more-or-less central in the box, the interfaces are combined
# in the analysis and oriented so that z<0 is towards the gas and z>0 is towards the liquid.
# The cosine is defined so that c<0 corresponds to z1<z2 and c>0 to z1>z2.

# For illustration and simplicity, we adopt a scheme of formatted files of the same kind
# as those that are saved at the end of each block of our MD simulation examples
# We assume that the initial configuration of a run has been copied to cnf.000
# and subsequent configurations are called cnf.001 cnf.002 etc., up to (at most) cnf.999
# Obviously, in a practical application, a binary trajectory file would fulfil this role.

# Cubic periodic boundary conditions are assumed
# r and box are assumed to be in the same units (e.g. LJ sigma)
# box is assumed to be unchanged throughout

# Values of basic parameters are read from standard input using JSON format

print('grint')
# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"dz":0.2, "nz":15, "nc":6, "dr":0.02, "z_mid":0.0, "iz_max":10, "zskip":5, "cskip":3, "nr":200 }
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ', list(defaults.keys()))
    
# Set parameters to input values or defaults
dz     = nml["dz"]     if "dz"     in nml else defaults["dz"]
nz     = nml["nz"]     if "nz"     in nml else defaults["nz"]
nc     = nml["nc"]     if "nc"     in nml else defaults["nc"]
dr     = nml["dr"]     if "dr"     in nml else defaults["dr"]
z_mid  = nml["z_mid"]  if "z_mid"  in nml else defaults["z_mid"]
iz_max = nml["iz_max"] if "iz_max" in nml else defaults["iz_max"]
zskip  = nml["zskip"]  if "zskip"  in nml else defaults["zskip"]
cskip  = nml["cskip"]  if "cskip"  in nml else defaults["cskip"]
nr     = nml["nr"]     if "nr"     in nml else defaults["nr"]

dc = 2.0 / (2*nc+1) # Cosine spacing to cover the range (-1,1)

# Write out parameters
print ( "{:40}{:15.6f}".format('Spacing in z-direction',       dz)     )
print ( "{:40}{:15d}  ".format('Number of z points',           nz)     )
print ( "{:40}{:15.6f}".format('+/- zmax',                     nz*dz)  )
print ( "{:40}{:15.6f}".format('Spacing in cos(theta)',        dc)     )
print ( "{:40}{:15d}  ".format('Number of cos(theta) points',  nc)     )
print ( "{:40}{:15.6f}".format('Spacing in r',                 dr)     )
print ( "{:40}{:15d}  ".format('Number of r points',           nr)     )
print ( "{:40}{:15.6f}".format('rmax',                         nr*dr)  )
print ( "{:40}{:15d}  ".format('Output z skip',                zskip)  )
print ( "{:40}{:15d}  ".format('Output z max',                 iz_max) )
print ( "{:40}{:15d}  ".format('Output c skip',                cskip)  )
print ( "{:40}{:15.6f}".format('Liquid slab midpoint (guess)', z_mid)  )

# Read in configuration
cnf_prefix = 'cnf.'
if not os.path.isfile(cnf_prefix+'000'):
    print(cnf_prefix+'000 does not exist')
    sys.exit()
n, box, r = read_cnf_atoms(cnf_prefix+'000')

# We must remember that artefacts are expected whenever z approaches the "other" interface
# This depends on widths of the two phases, on nz*dz, and on nr*dr
# It's up to you if you ignore this warning
if ( nz*dz + nr*dr ) > 0.25*box:
    print ( "{:40}{:15.6f}{:15.6f}".format('Warning: max z > box/4 = ', (nz*dz+nr*dr), 0.25*box)  )

# We define dz_box to fit the box exactly
nk     = np.rint(box/dz).astype(np.int)
dz_box = box / nk
area   = box**2
z_box  = np.linspace ( (-box+dz_box)/2, (box-dz_box)/2, nk )

z_max = ( nz+0.5 ) * dz # +/- binning range of z around interface
z     = np.linspace ( -nz*dz, nz*dz, 2*nz+1 ) # Coordinates around interface

cos_vals = np.linspace ( -1.0+dc/2, 1.0-dc/2, 2*nc+1 )
r_vals   = np.linspace ( dr/2, (nr-0.5)*dr, nr )

# Zero the accumulator arrays
dens = np.zeros ( nk, dtype=np.float_ )
rho1 = np.zeros ( 2*nz+1, dtype=np.float_ )
rho2 = np.zeros ( ( 2*nz+1, 2*nc+1, nr ), dtype=np.float_ )

norm = 0
t    = 0

# Initial guesses at slab fit parameters
# These will be passed on at each step, assuming that changes are small
# 5 parameters are z_gl, z_lg, width, rho_g, rho_l
c2tanh = [ -0.25*box, 0.25*box, 1.0, 0.0, 0.8 ]

while True: # Loop until configurations or naming scheme exhausted
    if t >= 1000:
        break
    sav_tag   = str(t).zfill(3)
    file_name = cnf_prefix+sav_tag
    if not os.path.isfile(file_name):
        break
    print('Processing file '+file_name)
    n, box, r = read_cnf_atoms(file_name)
    r[:,2] = r[:,2] - z_mid         # Place liquid slab approximately in centre of box
    r      = r - np.rint(r/box)*box # Apply periodic boundary conditions

    d, edges      = np.histogram ( r[:,2],bins=nk,range=(-box/2.0,box/2.0) ) # Histogram density
    dens          = dens + d
    d             = d / ( area*dz_box )
    c2tanh, covar = curve_fit ( f2tanh, z_box, d, c2tanh ) # Fit profile to estimate z_gl, z_lg

    rz = np.copy(r[:,2]) # Save z-coordinates of all atoms

    # Process gas-liquid interface
    r[:,2] = rz - c2tanh[0]                   # Shift gas-liquid interface to origin
    r[:,2] = r[:,2] - np.rint(r[:,2]/box)*box # Apply PBC
    d1, edges = np.histogram ( r[:,2],bins=2*nz+1,range=(-z_max,z_max) ) # Histogram density
    rho1      = rho1 + d1
    rho2      = rho2 + calculate(r)

    # Process liquid-gas interface
    r[:,2]    = c2tanh[1] - rz                   # Shift liquid-gas interface to origin and reflect
    r[:,2]    = r[:,2] - np.rint(r[:,2]/box)*box # Apply PBC
    d1, edges = np.histogram ( r[:,2],bins=2*nz+1,range=(-z_max,z_max) ) # Histogram density
    rho1      = rho1 + d1
    rho2      = rho2 + calculate(r)

    norm  = norm + 1
    t     = t + 1 # Ready for next file
    z_mid = z_mid + 0.5*(c2tanh[0]+c2tanh[1]) # Refine estimate of slab midpoint for next time

# Normalize (including factor for 2 interfaces)
dens = dens / ( norm * area * dz_box )
rho1 = rho1 / ( 2.0 * norm * area * dz )
rho2 = rho2 / ( 2.0 * 2.0*np.pi * norm * area * dr * dc * dz )

# Fit the averaged density profile
c2tanh, covar = curve_fit ( f2tanh, z_box, dens, c2tanh )

# Fit the single particle density
# 4 parameters are z_gl, width, rho_g, rho_l
c1tanh = [ 0.0, c2tanh[2], c2tanh[3], c2tanh[4] ]
c1tanh, covar = curve_fit ( f1tanh, z, rho1, c1tanh )

# Convert rho2 to g2, normalizing by fitted single-particle densities at z1 and z2
g2 = np.empty_like(rho2)

for iz, z1 in enumerate(z): # Loop over z coordinates around interface
    rho1_z1 = f1tanh ( z1, *c1tanh ) # Use fitted single-particle density (an approximation)

    for ic, c in enumerate(cos_vals): # Loop over cos(theta)

        for ir, rij_mag in enumerate(r_vals): # Loop over radial distance
            z2      = z1 - c * rij_mag
            rho1_z2 = f1tanh ( z2, *c1tanh ) # Use fitted single-particle density (an approximation)

            g2[iz,ic,ir] = rho2[iz,ic,ir] / ( rho1_z1 * rho1_z2 ) 

print('Box average density profile output to den.out')
fit = np.array([f2tanh(zval,*c2tanh) for zval in z_box])
np.savetxt('den.out',np.column_stack((z_box,dens,fit)),fmt="%15.8f")

print('Single-particle density profile rho1 output to rho.out')
fit = np.array([f1tanh(zval,*c1tanh) for zval in z])
np.savetxt('rho.out',np.column_stack((z,rho1,fit)),fmt="%15.8f")

print('Pair distribution function g2 output in selected slices')
print('Each slice has fixed z=z1 and c=cos(theta)')
print('Filenames have the form g2_ZZZ_CCC.out')

iz_max = min(iz_max,nz)
iz_max = iz_max - ( iz_max % zskip ) # Round to multiple of zskip
assert iz_max < 100, 'The output filename format will only cope with iz_max<100'
print('ZZZ             z1')
for iz in range(-iz_max,iz_max+1,zskip):
    print("{:+03d}{:+15.5f}".format(iz,z[nz+iz]))
print('-ve sign means z1 on gas side, +ve sign means z1 on liquid side')

ic_max = nc - ( nc % cskip ) # Round to multiple of cskip
assert ic_max < 100, 'The output filename format will only cope with ic_max<100'
print('CCC     cos(theta)')
for ic in range(-ic_max,ic_max+1, cskip):
    print("{:+03d}{:+15.5f}".format(ic,cos_vals[nc+ic]))
print('-ve sign means z1<z2, +ve sign means z1>z2')

for iz in range(-iz_max,iz_max+1,zskip):
    ztag = "{:+03d}".format(iz)
    for ic in range(-ic_max,ic_max+1, cskip):
        ctag = "{:+03d}".format(ic)
        filename = 'g2_'+ztag+'_'+ctag+'.out'
        np.savetxt(filename,np.column_stack((r_vals,g2[nz+iz,nc+ic,:])),fmt="%15.8f")
