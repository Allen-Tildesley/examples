#!/usr/bin/env python3
# sample_mean.py

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

"""sample-mean program to illustrate Monte Carlo evaluation of an integral."""

# No parameters need be supplied by the user. The exact value of the integral is 5/3.
# For details, see Chapter 4 of the text.

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

