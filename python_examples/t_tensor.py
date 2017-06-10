#!/usr/bin/env python3
# t_tensor.py

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

"""Electrostatic interactions: T-tensors compared with angles."""

# The dipole moment of molecule 1 is aligned along the axial vector e1
# The quadrupole tensor, quad1, is diagonal and traceless with
# quad1_xx = -0.5*quad1_mag, quad1_yy = -0.5*quad1_mag, and quad1_zz = quad1_mag in the molecule-fixed system.
# Similarly for molecule 2
# The vector r12 = r1-r2 points from 2 to 1.

# Forces are calculated by differentiating the T-tensor, giving the next higher rank T-tensor
# Torques are calculated from the angular dependence of dipole, quadrupole etc.
# potential V = mu_i g_i => torque tau_i = -epsilon_ijk mu_j g_k (=-cross_product)
# potential V = Q_ij G_ij => torque tau_l = -2 epsilon_lij Q_ik G_jk
# where ijkl are Cartesian indices and epsilon is the Levi-Civita symbol
# It is just necessary to identify the constants g_i, G_ij, in terms of the T tensor and the
# multipole on the other molecule.

# NB in the text, eqn (1.15), the signs of the odd-rank terms in the energy are wrong.
# See https://github.com/Allen-Tildesley/corrections. The correct formulae are used here.

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

def t4_tensor ( r, r5 ):
    """Returns fourth-rank 3x3x3x3 interaction tensor

    Supplied arguments should be the unit vector from 2 to 1 and
    the fifth power of the modulus of that vector.
    """

    import numpy as np
    from itertools import product

    # Define 3x3 unit matrix or Kronecker delta
    u = np.zeros((3,3),dtype=np.float_)
    u[0,0]=u[1,1]=u[2,2]=1.0

    t4 = 105.0 * np.einsum('i,j,k,l->ijkl',r,r,r,r) # Starting point: outer product of four vectors

    for i,j,k,l in product ( range(3), repeat=4 ):
        t4[i,j,k,l] = t4[i,j,k,l] - 15.0 * (
            r[i] * r[j] * u[k,l] + r[i] * r[k] * u[j,l] 
            + r[i] * r[l] * u[j,k] + r[j] * r[k] * u[i,l] 
            + r[j] * r[l] * u[i,k] + r[k] * r[l] * u[i,j]
            ) + 3.0 * ( u[i,j] * u[k,l] + u[i,k] * u[j,l] + u[i,l] * u[j,k] )

    t4 = t4 / r5 # Scale by fifth power of distance

    return t4

def t5_tensor ( r, r6 ):
    """Returns fifth-rank 3x3x3x3x3 interaction tensor

    Supplied arguments should be the unit vector from 2 to 1 and
    the sixth power of the modulus of that vector.
    """

    import numpy as np
    from itertools import product

    # Define 3x3 unit matrix or Kronecker delta
    u = np.zeros((3,3),dtype=np.float_)
    u[0,0]=u[1,1]=u[2,2]=1.0

    t5 = 945.0 * np.einsum('i,j,k,l,m->ijklm',r,r,r,r,r) # Starting point: outer product of five vectors

    for i,j,k,l,m in product ( range(3), repeat=5 ):
        t5[i,j,k,l,m] = t5[i,j,k,l,m] - 105.0 * (
            r[i] * r[j] * r[k] * u[l,m] + r[i] * r[j] * r[l] * u[k,m]     
            + r[i] * r[j] * r[m] * u[k,l] + r[i] * r[k] * r[l] * u[j,m]     
            + r[i] * r[k] * r[m] * u[j,l] + r[i] * r[l] * r[m] * u[j,k]     
            + r[j] * r[k] * r[l] * u[i,m] + r[j] * r[k] * r[m] * u[i,l]     
            + r[j] * r[l] * r[m] * u[i,k] + r[k] * r[l] * r[m] * u[i,j]
            ) + 15.0 * ( 
            r[i] * ( u[j,k] * u[l,m] + u[j,l] * u[k,m] + u[j,m] * u[k,l] ) 
            + r[j] * ( u[i,k] * u[l,m] + u[i,l] * u[k,m] + u[i,m] * u[k,l] ) 
            + r[k] * ( u[i,j] * u[l,m] + u[i,l] * u[j,m] + u[i,m] * u[j,l] ) 
            + r[l] * ( u[i,j] * u[k,m] + u[i,k] * u[j,m] + u[i,m] * u[j,k] ) 
            + r[m] * ( u[i,j] * u[k,l] + u[i,k] * u[j,l] + u[i,l] * u[j,k] )
            )

    t5 = t5 / r6 # Scale by sixth power of distance

    return t5

def skew ( a ):
    """Returns contraction of supplied 3x3 matrix with Levi-Civita tensor."""

    b = np.empty(3,dtype=np.float_)
    b[0] = a[1,2] - a[2,1]
    b[1] = a[2,0] - a[0,2]
    b[2] = a[0,1] - a[1,0]
    return b

import json
import sys
import numpy as np
from maths_module import random_vector

print('t_tensor')
print('Calculation of electrostatic interactions between linear molecules')
print('using T-tensors and Euler angles')

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"d_min":0.5, "d_max":1.5, "mu1_mag":1.0, "mu2_mag":1.0, "quad1_mag":1.0, "quad2_mag":1.0}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ', list(defaults.keys()))
    
# Set parameters to input values or defaults
d_min     = nml["d_min"]     if "d_min"     in nml else defaults["d_min"]     # Minimum separation
d_max     = nml["d_max"]     if "d_max"     in nml else defaults["d_max"]     # Maximum separation
mu1_mag   = nml["mu1_mag"]   if "mu1_mag"   in nml else defaults["mu1_mag"]   # Dipole moment of molecule 1
mu2_mag   = nml["mu2_mag"]   if "mu2_mag"   in nml else defaults["mu2_mag"]   # Dipole moment of molecule 2
quad1_mag = nml["quad1_mag"] if "quad1_mag" in nml else defaults["quad1_mag"] # Quadrupole moment of molecule 1
quad2_mag = nml["quad2_mag"] if "quad2_mag" in nml else defaults["quad2_mag"] # Quadrupole moment of molecule 2

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

# The T tensors of each rank: T2, T3, T4, T5
tt2 = t2_tensor ( r12_hat, r12_mag**3 )
tt3 = t3_tensor ( r12_hat, r12_mag**4 ) 
tt4 = t4_tensor ( r12_hat, r12_mag**5 ) 
tt5 = t5_tensor ( r12_hat, r12_mag**6 ) 

# Headings
print()
print("{:>66}{:>40}{:>40}".format('.....Result from T tensor','.....Result from Euler angles','.........Difference') )

print('\nDipole-dipole')
e_fmt = "{:30}{:36.6f}{:40.6f}{:40.2e}"
f_fmt = "{:30}{:12.6f}{:12.6f}{:12.6f}{:16.6f}{:12.6f}{:12.6f}{:16.2e}{:12.2e}{:12.2e}"

# Calculate the dipole-dipole energy
v12t = -np.einsum('i,ij,j',mu1,tt2,mu2 ) # Contract both dipoles with T2
v12e = (mu1_mag*mu2_mag/r12_mag**3) * ( c12-3.0*c1*c2 )
print(e_fmt.format('Energy',v12t,v12e,v12t-v12e) )

# Calculate the dipole-dipole force
f12t = -np.einsum ( 'ijk,j,k->i',tt3,mu1,mu2) # Contract T3 with both dipoles
f12e = (3.0*mu1_mag*mu2_mag/r12_mag**4) * ( (c12-5.0*c1*c2)*r12_hat + c2*e1 + c1*e2 )
print(f_fmt.format('Force',*np.concatenate((f12t,f12e,f12t-f12e)) ))

# Calculate the dipole-dipole torques
g   = -np.einsum( 'ij,j->i',tt2,mu2 ) # Contract T2 with dipole 2
t1t = -np.cross ( mu1, g )            # Cross-product result with dipole 1
g   = e2 - 3.0*c2*r12_hat             # Compare result from angles
t1e = -(mu1_mag*mu2_mag/r12_mag**3) * np.cross ( e1, g )
print(f_fmt.format('Torque on 1',*np.concatenate((t1t,t1e,t1t-t1e)) ))
g   = -np.einsum( 'ij,j->i',tt2,mu1) # Contract T2 with dipole 1
t2t = -np.cross ( mu2, g )           # Cross-product result with dipole 2
g   = e1 - 3.0*c1 * r12_hat          # Compare result from angles
t2e = -(mu1_mag*mu2_mag/r12_mag**3) * np.cross ( e2, g )
print(f_fmt.format('Torque on 2',*np.concatenate((t2t,t2e,t2t-t2e)) ))

print('\nDipole-quadrupole')

# Calculate the dipole-quadrupole energy
v12t = -(1.0/3.0) * np.einsum ('i,ijk,jk',mu1,tt3,quad2) # Contract dipole 1 with T3 and quadrupole 2
v12e = (1.5*mu1_mag*quad2_mag/r12_mag**4) * ( c1*(1.0-5.0*c2**2) + 2.0*c2*c12 ) 
print(e_fmt.format('Energy',v12t,v12e,v12t-v12e) )

# Calculate the dipole-quadrupole force
f12t = -(1.0/3.0) * np.einsum( 'ijkl,j,kl->i', tt4, mu1, quad2 ) # Contract T4 with dipole 1 and quadrupole 2
f12e = -(1.5*mu1_mag*quad2_mag/r12_mag**5) * (  # Compare result from angles
          ( 35.0*c1*c2**2 - 10.0*c2*c12 - 5.0*c1 ) * r12_hat 
            + ( 1.0 - 5.0*c2**2 ) * e1
            + ( 2.0*c12 - 10.0*c1*c2 ) * e2
          )
print(f_fmt.format('Force',*np.concatenate((f12t,f12e,f12t-f12e)) ))

# Calculate the dipole-quadrupole torques
g   = -(1.0/3.0)*np.einsum('ijk,jk->i', tt3, quad2 ) # Contract T3 with quadrupole 2
t1t = -np.cross ( mu1, g )                           # Cross-product result with dipole 1
g   =  (1.0-5.0*c2**2) * r12_hat + 2.0*c2 * e2       # Compare result from angles
t1e = -(1.5*mu1_mag*quad2_mag/r12_mag**4) * np.cross ( e1, g )
print(f_fmt.format('Torque on 1',*np.concatenate((t1t,t1e,t1t-t1e)) ))
gg  = -(1.0/3.0)*np.einsum('ijk,k->ij', tt3, mu1 ) # Contract T3 with dipole 1
gg  = np.einsum('ik,jk->ij', quad2, gg )           # Contract result with quadrupole 2
t2t = -2.0*skew ( gg )                             # Contract with Levi-Civita symbol
g   =  (c12-5.0*c1*c2) * r12_hat + c2 * e1 # Compare result from angles
t2e = -(3.0*mu1_mag*quad2_mag/r12_mag**4) * np.cross ( e2, g )
print(f_fmt.format('Torque on 2',*np.concatenate((t2t,t2e,t2t-t2e)) ))

print('\nQuadrupole-dipole')

# Calculate the quadrupole-dipole energy
v12t = (1.0/3.0) * np.einsum ('i,ijk,jk',mu2,tt3,quad1) # Contract dipole on 2 with T3 and quadrupole on 1
v12e = -(1.5*quad1_mag*mu2_mag/r12_mag**4) * ( c2*(1.0-5.0*c1**2) + 2.0*c1*c12 ) 
print(e_fmt.format('Energy',v12t,v12e,v12t-v12e) )

# Calculate the quadrupole-dipole force
f12t = (1.0/3.0) *  np.einsum( 'ijkl,j,kl->i', tt4, mu2, quad1 ) # Contract T4 with quadrupole 1 and dipole 2
f12e = (1.5*quad1_mag*mu2_mag/r12_mag**5) * (  # Compare result from angles
        ( 35.0*c2*c1**2 - 10.0*c1*c12 - 5.0*c2 ) * r12_hat 
        + ( 1.0-5.0*c1**2 ) * e2
        + ( 2.0*c12 - 10.0*c1*c2 ) * e1
        )
print(f_fmt.format('Force',*np.concatenate((f12t,f12e,f12t-f12e)) ))

# Calculate the quadrupole-dipole torques
gg  = (1.0/3.0)*np.einsum( 'ijk,k->ij', tt3, mu2 ) # Contract T3 with dipole 2
gg  = np.einsum('ik,jk->ij', quad1, gg )           # Contract result with quadrupole 1
t1t = -2.0*skew ( gg )                             # Contract with Levi-Civita symbol
g   = (c12-5.0*c1*c2) * r12_hat + c1 * e2 # Compare result from angles
t1e = (3.0*quad1_mag*mu2_mag/r12_mag**4) * np.cross ( e1, g )
print(f_fmt.format('Torque on 1',*np.concatenate((t1t,t1e,t1t-t1e)) ))
g   = (1.0/3.0)*np.einsum( 'ijk,jk->i', tt3, quad1 ) # Contract T3 with quadrupole 1
t2t = -np.cross ( mu2, g )                           # Cross product result with dipole 2
g   = (1.0-5.0*c1**2) * r12_hat + 2.0*c1 * e1 # Compare result from angles
t2e = (1.5*quad1_mag*mu2_mag/r12_mag**4) * np.cross ( e2, g )
print(f_fmt.format('Torque on 2',*np.concatenate((t2t,t2e,t2t-t2e)) ))

print('\nQuadrupole-quadrupole')

# Calculate the quadrupole-quadrupole energy
v12t = (1.0/9.0) * np.einsum( 'ijkl,ij,kl', tt4, quad1, quad2 ) # Contract T4 with both quadrupoles
v12e = (0.75*quad1_mag*quad2_mag/r12_mag**5) * (  # Compare result from angles
          1.0 - 5.0*c1**2 - 5.0*c2**2 + 2.0*c12**2 + 35.0*(c1*c2)**2 - 20.0*c1*c2*c12
        )
print(e_fmt.format('Energy',v12t,v12e,v12t-v12e) )

# Calculate the quadrupole-quadrupole force
f12t = (1.0/9.0) * np.einsum( 'ijklm,jk,lm->i', tt5, quad1, quad2 ) # Contract T5 with both quadrupoles
f12e = (0.75*quad1_mag*quad2_mag/r12_mag**6) * (  # Compare result from angles
           ( 5.0 - 35.0*c1**2 - 35.0*c2**2 + 10.0*c12**2 + 315.0*(c1*c2)**2 - 140.0*c1*c2*c12 ) * r12_hat 
         + ( 10.0*c1 - 70.0*c1*c2**2 + 20.0*c2*c12 ) * e1
         + ( 10.0*c2 - 70.0*c2*c1**2 + 20.0*c1*c12 ) * e2
        )          
print(f_fmt.format('Force',*np.concatenate((f12t,f12e,f12t-f12e)) ))

# Calculate the quadrupole-quadrupole torques
gg  = (1.0/9.0)*np.einsum( 'ijkl,kl->ij', tt4, quad2 ) # Contract T4 with quadrupole 2
gg  = np.einsum( 'ik,jk->ij', quad1, gg )              # Contract result with quadrupole 1
t1t = -2.0*skew ( gg )                                 # Contract with Levi-Civita symbol
g   = 2.5*(c1*(7.0*c2**2-1.0)-2.0*c2*c12) * r12_hat - (5.0*c1*c2-c12) * e2 # Compare result from angles
t1e = -(3.0*quad1_mag*quad2_mag/r12_mag**5) * np.cross ( e1, g )
print(f_fmt.format('Torque on 1',*np.concatenate((t1t,t1e,t1t-t1e)) ))
gg  = (1.0/9.0)*np.einsum( 'ijkl,kl->ij', tt4, quad1 ) # Contract T4 with quadrupole 1
gg  = np.einsum( 'ik,jk->ij', quad2, gg )              # Contract result with quadrupole 2
t2t = -2.0*skew ( gg )                                 # Contract with Levi-Civita symbol
g   = 2.5*(c2*(7.0*c1**2-1.0)-2.0*c1*c12) * r12_hat -(5.0*c1*c2-c12) * e1 # Compare result from angles
t2e = -(3.0*quad1_mag*quad2_mag/r12_mag**5) * np.cross ( e2, g )
print(f_fmt.format('Torque on 2',*np.concatenate((t2t,t2e,t2t-t2e)) ))
