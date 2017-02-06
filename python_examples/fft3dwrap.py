#!/usr/bin/env python3
# fft3dwrap program

import json
import sys
import numpy as np

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
dk = (2.0 * np.pi) / dr / sc # interval in reciprocal space
print ( "{:40}{:15d}  ".format('Grid points in each dimension (sc)',    sc)  )
print ( "{:40}{:15.6f}".format('Periodic repeat length (box)',          box) )
print ( "{:40}{:15.6f}".format('Grid spacing in real space (dr)',       dr)  )
print ( "{:40}{:15.6f}".format('Grid spacing in reciprocal space (dk)', dk)  )

rmesh = np.linspace ( 0.0, box, num=sc, endpoint=False, dtype=np.float_ )
rmesh[sc2:] = rmesh[sc2:] - box # wraparound positions
kbox = (2.0 * np.pi) / dr
kmesh = np.linspace ( 0.0, kbox, num=sc, endpoint=False, dtype=np.float_ )
kmesh[sc2:] = kmesh[sc2:] - kbox # wraparound k values

# Set up initial arrays
r=np.zeros(3,dtype=np.float_)
k=np.zeros(3,dtype=np.float_)
fft_inp = np.zeros((sc,sc,sc),dtype=np.complex_)
out_max=15

print("Initial real-space Gaussian")
print("{:>15}{:>15}{:>15}{:>15}{:>15}".format('   ix   iy   iz', '|r|', 'Gaussian(r)', 'FFT (real)', 'FFT (imag)'))

# Triple loop over xyz grid points
for ix in range(sc):
    r[0]=rmesh[ix]
    for iy in range(sc):
        r[1]=rmesh[iy]
        for iz in range(sc):
            r[2]=rmesh[iz]
            r_sq = np.sum(r**2) # Squared distance from origin
            g = np.exp(-np.pi*r_sq) # Setup 3D Gaussian (decay parameter chosen to be pi)
            fft_inp[ix,iy,iz] = g
            if ix**2+iy**2+iz**2 <= out_max:
                f=fft_inp[ix,iy,iz]
                r_mag = np.sqrt(r_sq)
                print("{:5d}{:5d}{:5d}{:15.6f}{:15.6f}{:15.6f}{:15.6f}".format(ix,iy,iz,r_mag,g,f.real,f.imag))

# Forward FFT
fft_out = np.fft.fftn(fft_inp)

print("Reciprocal-space transform")
print("{:>15}{:>15}{:>15}{:>15}{:>15}".format('   ix   iy   iz', '|k|', 'Gaussian(k)', 'FFT (real)', 'FFT (imag)'))

# Triple loop over xyz grid points
for ix in range(sc):
    k[0]=kmesh[ix]
    for iy in range(sc):
        k[1]=kmesh[iy]
        for iz in range(sc):
            k[2]=kmesh[iz]
            if ix**2+iy**2+iz**2 <= out_max:
                k_sq = np.sum(k**2) # Squared magnitude of wave vector
                g = np.exp(-k_sq/(4.0*np.pi)) # Analytical transform of Gaussian
                f=fft_out[ix,iy,iz]*dr**3
                k_mag = np.sqrt(k_sq)
                print("{:5d}{:5d}{:5d}{:15.6f}{:15.6f}{:15.6f}{:15.6f}".format(ix,iy,iz,k_mag,g,f.real,f.imag))

# Backward FFT
fft_inp = np.fft.ifftn(fft_out)

# Write some elements of data in real space after the back transform
# For NumPy FFT, the normalising factor 1/sc**3 for the inverse transform is built in
# Compare with the (real) input data

print("Back transform to real space")
print("{:>15}{:>15}{:>15}{:>15}{:>15}".format('   ix   iy   iz', '|r|', 'Gaussian(r)', 'FFT (real)', 'FFT (imag)'))

# Triple loop over xyz grid points
for ix in range(sc):
    r[0]=rmesh[ix]
    for iy in range(sc):
        r[1]=rmesh[iy]
        for iz in range(sc):
            r[2]=rmesh[iz]
            if ix**2+iy**2+iz**2 <= out_max:
                r_sq = np.sum(r**2) # Squared distance from origin
                g = np.exp(-np.pi*r_sq) # Original 3D Gaussian
                f=fft_inp[ix,iy,iz]
                r_mag = np.sqrt(r_sq)
                print("{:5d}{:5d}{:5d}{:15.6f}{:15.6f}{:15.6f}{:15.6f}".format(ix,iy,iz,r_mag,g,f.real,f.imag))
