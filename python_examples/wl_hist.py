#!/usr/bin/env python3
# wl_hist.py

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

"""Wang-Landau histogram post-processing."""

# Program to post-process histograms produced by mc_chain_wl_sw.f90
# Reads the following from standard input, until end-of-file:
#    histogram filename (20 characters max)
#    number of atoms in chain (solely to make PE and Cv intensive)
#    one temperature per line (as many lines as you like)

import sys
import numpy as np
import math

hist_file = input("Enter histogram file name: ")
print('Reading histograms from '+hist_file)
e,h,g,s = np.loadtxt(hist_file,dtype=np.float_,unpack=True)

n = input("Enter number of atoms: ")
n = int(n)
print("{:10}{:15d}".format('N = ',n))
print('Enter temperatures one by one')

print("{:>15}{:>15}{:>15}{:>15}".format('T', 'Rg', 'PE/N', 'Cv(ex)/N'))

for line in sys.stdin:
    t = float(line)

    # Locate maximum Boltzmann factor (helps avoid overflows)
    qs_max = np.argmax ( s - e / t )
    s_max  = s[qs_max]

    # Compute Boltzmann weights including density of states
    boltz = s - s_max - e / t
    boltz = np.exp ( boltz )

    # Calculate averages
    norm  = np.sum ( boltz )
    g_avg = np.sum ( boltz * g ) / norm
    e_avg = np.sum ( boltz * e ) / norm
    e_msd = np.sum ( boltz * (e - e_avg )**2  ) / norm
    e_avg = e_avg / n            # Energy per atom
    e_msd = e_msd / ( n * t**2 ) # Heat capacity per atom

    print("{:15.6f}{:15.6f}{:15.6f}{:15.6f}".format(t, g_avg, e_avg, e_msd))
