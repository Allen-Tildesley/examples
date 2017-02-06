#!/usr/bin/env python3
""" hit_and_miss program to illustrate Monte Carlo evaluation of an integral.

No parameters need be supplied by the user. The exact value of the integral is 5/3.
For details, see Chapter 4 of the text.
"""

import numpy as np

print('hit_and_miss')
np.random.seed()

r_0 = np.array([1.0,2.0,3.0],dtype='f8')
v_0 = np.prod ( r_0 )

tau_hit  = 0
tau_shot = 1000000

for tau in range(tau_shot):
    zeta = np.random.rand(3) # uniform in range (0,1)
    r = zeta * r_0           # uniform in v_0
    if r[1] < ( 2.0 - 2.0*r[0] ) and r[2] < ( 1.0 + r[1] ) :
        tau_hit += 1 # in polyhedron

v = v_0 * tau_hit / tau_shot
print ( "{}{:10.5f}".format('Estimate =', v) )

