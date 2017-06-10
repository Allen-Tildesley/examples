#!/usr/bin/env python3
# eos_hs.py

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

"""Equation of State for hard sphere potential."""

import json
import sys
import numpy as np

# This program uses the function derived by H Hansen-Goos, J Chem Phys, 16, 164506 (2016)
# which is claimed to be an excellent fit to simulation data over the whole fluid density range
# That paper also gives references to previous approximate equations of state (such as the
# venerable Carnahan-Starling equation).

# The coefficients appear in Table I of Hansen-Goos (2016)

a = 8.0
b = [ -0.096495, 0.248245, 0.041212, -1.265575, -2.635232, 47.0/3.0, -19.0, 9.0 ]

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"density":0.75}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ', list(defaults.keys()))
    
# Set parameters to input values or defaults
density = nml["density"] if "density" in nml else defaults["density"]

# Write out parameters
eta = np.pi * density / 6.0 # Packing fraction
print ( "{:40}{:15.6f}".format('Density rho',          density) )
print ( "{:40}{:15.6f}".format('Packing fraction eta', eta    ) )

# Equation (6) of Hansen-Goos (2016)
z = a * np.log ( 1.0-eta ) / eta
z = z + np.polyval ( b, eta ) / ( 1.0 - eta ) ** 3 # Compressibility factor P/(rho*kT)
p = z * density                                    # Pressure P / kT
print ( "{:40}{:15.6f}".format('Pressure P',                     p ) )
print ( "{:40}{:15.6f}".format('Compressibility factor Z = P/(rho*kT)', z ) )
