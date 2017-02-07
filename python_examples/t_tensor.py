#!/usr/bin/env python3
# t_tensor program
"""Electrostatic interactions: T-tensors compared with angles.

The dipole moment of molecule 1 is aligned along the axial vector e1
The quadrupole tensor, quad1, is diagonal and traceless with
quad1_xx = -0.5*quad1_mag, quad1_yy = -0.5*quad1_mag, and quad1_zz = quad1_mag in the molecule-fixed system.
Similarly for molecule 2
The vector r12 = r1-r2 points from 2 to 1.
"""

def t2_tensor ( r, r3 ):
    """Returns second-rank 3x3 interaction tensor.

    Supplied arguments should be the unit vector from 2 to 1 and
    the third power of the modulus of that vector.
    """

    import numpy as np
    t2 = 3.0 * np.outer(r,r) # Starting point: outer product
    t2 = t2 - np.identity(3) # Make traceless

    t2 = t2 / r3 # Scale by third power of distance

    return t2

def t3_tensor ( r, r4 ):
    """Returns third-rank 3x3x3 interaction tensor (note positive sign).

    Supplied arguments should be the unit vector from 2 to 1 and
    the fourth power of the modulus of that vector.
    """

    import numpy as np

    t3 = 15.0 * np.einsum('i,j,k->ijk',r,r,r) # Starting point: outer product of three vectors

    for i in range(3):
        t3[i,i,i] = t3[i,i,i] - 9.0 * r[i] # Correction for all indices the same
        for j in range(3):
            if j == i:
                continue
            t3[i,i,j] = t3[i,i,j] - 3.0 * r[j] # Correction for two indices the same
            t3[i,j,i] = t3[i,j,i] - 3.0 * r[j] # Correction for two indices the same
            t3[j,i,i] = t3[j,i,i] - 3.0 * r[j] # Correction for two indices the same

    t3 = t3 / r4 # Scale by fourth power of distance

    return t3

import json
import sys
import numpy as np
from maths_module import random_vector

print('t_tensor')
print('Calculation of electrostatic interactions between linear molecules')
print('using T-tensors and Euler angles')

# Read parameters in JSON format
allowed_nml_keys = ["d_min","d_max","mu1_mag","mu2_mag","quad1_mag","quad2_mag"]

try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

for key in nml:
    if key not in allowed_nml_keys:
        print('Warning', key, 'not in allowed_nml_keys',allowed_nml_keys)
    
# Set parameters to input values or defaults
d_min     = nml["d_min"]     if "d_min"     in nml else 0.5 # Minimum separation
d_max     = nml["d_max"]     if "d_max"     in nml else 1.5 # Maximum separation
mu1_mag   = nml["mu1_mag"]   if "mu1_mag"   in nml else 1.0 # Dipole moment of molecule 1
mu2_mag   = nml["mu2_mag"]   if "mu2_mag"   in nml else 1.0 # Dipole moment of molecule 2
quad1_mag = nml["quad1_mag"] if "quad1_mag" in nml else 1.0 # Quadrupole moment of molecule 1
quad2_mag = nml["quad2_mag"] if "quad2_mag" in nml else 1.0 # Quadrupole moment of molecule 2

# Write out parameters
print ( "{:40}{:15.6f}".format('Min separation d_min',            d_min)      )
print ( "{:40}{:15.6f}".format('Max separation d_max',            d_max)      )
print ( "{:40}{:15.6f}".format('Dipole moment of molecule 1',     mu1_mag)    )
print ( "{:40}{:15.6f}".format('Dipole moment of molecule 2',     mu2_mag)    )
print ( "{:40}{:15.6f}".format('Quadrupole moment of molecule 1', quad1_mag)  )
print ( "{:40}{:15.6f}".format('Quadrupole moment of molecule 2', quad2_mag)  )

np.random.seed()

# Choose orientations at random
e1 = random_vector()
e2 = random_vector()

# Place atom 2 at origin and atom 1 in a random direction within desired distance range
r12_hat = random_vector()
r12_mag = np.random.rand()
r12_mag = d_min + (d_max-d_min)*r12_mag # Magnitude of r12
r12     = r12_hat * r12_mag             # Within desired range of origin

r12_sq  = r12_mag**2       # Squared distance
r12_cu  = r12_sq * r12_mag # Cubed distance
r12_fo  = r12_sq * r12_sq  # Fourth power of distance

c1  = np.dot( e1, r12_hat ) # Cosine of angle between e1 and r12
c2  = np.dot( e2, r12_hat ) # Cosine of angle between e2 and r12
c12 = np.dot( e1, e2   )    # Cosine of angle between e1 and e2

print ( "{:40}{:10.6f}{:10.6f}{:10.6f}".format('Displacement r12', *r12)  )
print ( "{:40}{:10.6f}{:10.6f}{:10.6f}".format('Orientation e1',   *e1)   )
print ( "{:40}{:10.6f}{:10.6f}{:10.6f}".format('Orientation e2',   *e2)   )

# Dipole vectors in space-fixed frame
mu1 = mu1_mag * e1                             
mu2 = mu2_mag * e2

# Quadrupole tensors in space-fixed frame (traceless)
quad1 = 1.5 * np.outer ( e1, e1 )
quad1 = quad1 -0.5*np.identity(3)
quad1 = quad1_mag * quad1
quad2 = 1.5 * np.outer ( e2, e2 )
quad2 = quad2 -0.5*np.identity(3)
quad2 = quad2_mag * quad2

tt2 = t2_tensor ( r12_hat, r12_cu ) # Calculate the second rank tensor, T2
tt3 = t3_tensor ( r12_hat, r12_fo ) # Calculate the third rank tensor, T3

# Headings
print()
print("{:>60}{:>40}".format('.....Result from T tensor','.....Result from Euler angles') )

# Calculate the dipole-dipole energy
vddt = -np.einsum('i,ij,j',mu1,tt2,mu2 ) # Contract both dipoles with T2
vdde = mu1_mag * mu2_mag * ( c12 - 3.0 * c1 * c2 ) / r12_cu
print("{:30}{:30.6f}{:40.6f}".format('Dipole-dipole energy',vddt,vdde) )

# Calculate the dipole-quadrupole energy
vdqt = -(1.0/3.0) * np.einsum ('i,ijk,jk',mu1,tt3,quad2) # Contract dipole on 1 with T3 and quadrupole on 2
vdqe = 1.5 * mu1_mag * quad2_mag * ( c1 * (1.0 - 5.0 * c2 * c2) + 2.0 * c2 * c12 ) / r12_fo
print("{:30}{:30.6f}{:40.6f}".format('Dipole-quadrupole energy',vdqt,vdqe) )

# Calculate the quadrupole-dipole energy
vqdt = (1.0/3.0) * np.einsum ('i,ijk,jk',mu2,tt3,quad1) # Contract dipole on 2 with T3 and quadrupole on 1
vqde = - 1.5 * quad1_mag * mu2_mag * ( c2 * (1.0 - 5.0 * c1 * c1 ) + 2.0 * c1 * c12 ) / r12_fo
print("{:30}{:30.6f}{:40.6f}".format('Quadrupole-dipole energy',vqdt,vqde) )

# Calculate the dipole-dipole force
f12t = -np.einsum ( 'ijk,j,k->i',tt3,mu1,mu2) # Contract T3 with both dipoles
f12e = (3.0 * mu1_mag * mu2_mag / r12_fo ) * ((c12 - 5.0 *c1 * c2) * r12_hat + c2 * e1 + c1 * e2 )
print("{:30}{:10.6f}{:10.6f}{:10.6f}{:20.6f}{:10.6f}{:10.6f}".format('Dipole-dipole force',*f12t,*f12e) )
