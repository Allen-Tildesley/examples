#!/usr/bin/env python3
""" sample_mean program to illustrate Monte Carlo evaluation of an integral.

No parameters need be supplied by the user. The exact value of the integral is 5/3.
For details, see Chapter 4 of the text.
"""

import numpy as np

print('sample_mean')
np.random.seed()

r_0 = np.array([1.0,2.0],dtype=np.float_)
a_0 = np.prod ( r_0 )

f = 0.0
tau_max = 1000000

for tau in range(tau_max):
    zeta = np.random.rand(2) # uniform in range (0,1)
    r = zeta * r_0           # uniform in xy rectangle
    if r[1] < ( 2.0 - 2.0*r[0] ) :
        f += (1.0+r[1]) # value of z in xy triangle

v = a_0 * f / tau_max
print ( "{}{:10.5f}".format('Estimate =', v) )

