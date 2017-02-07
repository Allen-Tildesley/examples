#!/usr/bin/env python3
# fft3dwrap program

import json
import sys
import numpy as np
import math

print('fft3dwrap')

# Read parameters in JSON format
allowed_nml_keys = ["sc2","box"]

try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

for key in nml:
    if key not in allowed_nml_keys:
        print('Warning', key, 'not in allowed_nml_keys',allowed_nml_keys)
    
# Set parameters to input values or defaults
sc2 = nml["sc2"] if "sc2" in nml else 2**4 # Not essential to be a power of 2, but usually more efficient
box = nml["box"] if "box" in nml else 6.0  # Large enough to accommodate the chosen 3D Gaussian, for good comparison with analytical result

# Write out parameters
sc = sc2 * 2
dr = box / sc
dk = (2.0 * np.pi) / dr / sc
kbox = (2.0 * np.pi) / dr
print ( "{:40}{:15d}  ".format('Grid points in each dimension (sc)',    sc)  )
print ( "{:40}{:15.6f}".format('Periodic repeat length in r (box)',    box)  )
print ( "{:40}{:15.6f}".format('Grid spacing in real space (dr)',       dr)  )
print ( "{:40}{:15.6f}".format('Periodic repeat length in k (kbox)',   kbox) )
print ( "{:40}{:15.6f}".format('Grid spacing in reciprocal space (dk)', dk)  )

# Set limits for printing sample values
isq_max = 15
imax = min ( sc2, math.ceil(math.sqrt(isq_max)) )

# Set up indices in wraparound convention
imesh = np.arange(sc)
imesh[sc2:] = imesh[sc2:] - sc

# Set up mesh of points in r-space with subset for printing sample values
rmesh  = (dr*imesh).tolist()
rprint = rmesh[:imax]

# Set up mesh of points in k-space with subset for printing sample values
kmesh  = (dk*imesh).tolist()
kprint = kmesh[:imax]

# Set up initial arrays
# For a function of the form f(x)*f(y)*f(z) such as this, there are faster ways of initializing the array
# However, we retain the generic triple loop to make it easy to choose a different function

fft_inp = np.empty((sc,sc,sc),dtype=np.complex_)

print()
print("Initial real-space Gaussian")
print("{:>15}{:>15}{:>15}{:>15}{:>15}".format('   ix   iy   iz', '|r|', 'Gaussian(r)', 'FFTinp (real)', 'FFTinp (imag)'))

# Triple loop over all xyz grid points in r-space
for ix, rx in enumerate(rmesh):
    for iy, ry in enumerate(rmesh):
        for iz, rz in enumerate(rmesh):
            r_sq = rx**2+ry**2+rz**2    # Squared magnitude of r vector
            g = math.exp(-math.pi*r_sq) # Setup 3D Gaussian (decay parameter chosen to be pi)
            fft_inp[ix,iy,iz] = g       # Store input function
            if ix**2+iy**2+iz**2 <= isq_max:
                f = fft_inp[ix,iy,iz]
                r = math.sqrt(r_sq)
                print("{:5d}{:5d}{:5d}{:15.6f}{:15.6f}{:15.6f}{:15.6f}".format(ix,iy,iz,r,g,f.real,f.imag))

# Forward FFT
fft_out = np.fft.fftn(fft_inp)

print()
print("Reciprocal-space transform")
print("{:>15}{:>15}{:>15}{:>15}{:>15}".format('   ix   iy   iz', '|k|', 'Gaussian(k)', 'FFTout (real)', 'FFTout (imag)'))

# Triple loop over subset of xyz grid points in k-space
for ix, kx in enumerate(kprint):
    for iy, ky in enumerate(kprint):
        for iz, kz in enumerate(kprint):
            if ix**2+iy**2+iz**2 <= isq_max:
                k_sq = kx**2+ky**2+kz**2             # Squared magnitude of k vector
                g    = math.exp(-k_sq/(4.0*math.pi)) # Analytical transform of Gaussian
                f    = fft_out[ix,iy,iz]*dr**3       # FFT out, scaled by dr**3
                k    = math.sqrt(k_sq)               # Magnitude of k vector
                print("{:5d}{:5d}{:5d}{:15.6f}{:15.6f}{:15.6f}{:15.6f}".format(ix,iy,iz,k,g,f.real,f.imag))

# Backward FFT
fft_inp = np.fft.ifftn(fft_out)

# Write some elements of data in real space after the back transform
# For NumPy FFT, the normalising factor 1/sc**3 for the inverse transform is built in
# Compare with the (real) input data

print()
print("Back transform to real space")
print("{:>15}{:>15}{:>15}{:>15}{:>15}".format('   ix   iy   iz', '|r|', 'Gaussian(r)', 'FFTbak (real)', 'FFTbak (imag)'))

# Triple loop over subset of xyz grid points in r-space
for ix, rx in enumerate(rprint):
    for iy, ry in enumerate(rprint):
        for iz, rz in enumerate(rprint):
            if ix**2+iy**2+iz**2 <= isq_max:
                r_sq = rx**2+ry**2+rz**2 # Squared magnitude of r vector
                g    = math.exp(-math.pi*r_sq) # Original 3D Gaussian
                f    = fft_inp[ix,iy,iz]
                r    = math.sqrt(r_sq)
                print("{:5d}{:5d}{:5d}{:15.6f}{:15.6f}{:15.6f}{:15.6f}".format(ix,iy,iz,r,g,f.real,f.imag))
