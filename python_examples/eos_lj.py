#!/usr/bin/env python3
# eos_lj.py

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

"""Equation of State for Lennard-Jones pair potential."""

import json
import sys
import numpy as np
from eos_lj_module import a_res_full, a_res_cutshift
from lrc_module    import potential_lrc, pressure_lrc, pressure_delta

# The routines in the above module use the fitting function described and parametrized in
# M Thol, G Rutkai, R Span, J Vrabec, R Lustig, Int J Thermophys 36, 25 (2015)
# M Thol, G Rutkai, A Koester, R Lustig, R Span, J Vrabec, J Phys Chem Ref Data 45, 023101 (2016)
# Those authors also supply C++ codes (in the supplementary information of those papers)
# They are NOT responsible for this Fortran code, which was written independently by Michael P Allen
# A similar notation, consistent with the papers, is retained for clarity.

# Formulae for P, E/N etc in terms of the scaled free energy derivatives a_res(0,1) etc
# may be found in the above papers

r_cut = 2.5

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
defaults = {"temperature":1.0, "density":0.75}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ', list(defaults.keys()))
    
# Set parameters to input values or defaults
temperature = nml["temperature"] if "temperature" in nml else defaults["temperature"]
density     = nml["density"]     if "density"     in nml else defaults["density"]

# Write out parameters
print ( "{:40}{:15.6f}".format('Temperature T', temperature    ) )
print ( "{:40}{:15.6f}".format('Density rho',   density) )

# Results for full potential from Thol et al (2016) fitting formula
print('')
print('Full Lennard-Jones potential')
print('')

a_res = a_res_full ( temperature, density )

for (i,j), aij in np.ndenumerate ( a_res ):
    if i+j > 2: # Only interested in some of the results
        continue
    print ( "{:4}{:1d}{:<35d}{:15.6f}".format('Ares', i, j, aij ) )

p  = density * temperature * ( 1.0 + a_res[0,1] )
e  = temperature * ( 1.5 + a_res[1,0] )
cv = 1.5 - a_res[2,0]
cp = 2.5 - a_res[2,0]+(1.0+a_res[0,1]-a_res[1,1])*(1.0+a_res[0,1]-a_res[1,1])/(1.0+2.0*a_res[0,1]+a_res[0,2]) - 1.0
mu = temperature * ( np.log(density) + a_res[0,0] + a_res[0,1] )
z  = density * np.exp ( a_res[0,0] + a_res[0,1] )
print('')
print ( "{:40}{:15.6f}".format('Pressure P',            p  ) )
print ( "{:40}{:15.6f}".format('Energy E/N',            e  ) )
print ( "{:40}{:15.6f}".format('Heat capacity Cv/NkB',  cv ) )
print ( "{:40}{:15.6f}".format('Heat capacity Cp/NkB',  cp ) )
print ( "{:40}{:15.6f}".format('Chemical potential mu', mu ) )
print ( "{:40}{:15.6f}".format('Activity z',            z  ) )

# Estimates for cut (but not shifted) potential by reverse-application of long-range & delta corrections
print('')
print('Lennard-Jones potential cut (but not shifted) at 2.5 sigma')

p  = p - pressure_lrc ( density, r_cut ) + pressure_delta ( density, r_cut )
e  = e - potential_lrc ( density, r_cut )
mu = mu - 2.0 * potential_lrc ( density, r_cut )
z  = z * np.exp ( -2.0* potential_lrc ( density, r_cut ) / temperature )
print('')
print ( "{:40}{:15.6f}".format('Pressure P',            p  ) )
print ( "{:40}{:15.6f}".format('Energy E/N',            e  ) )
print ( "{:40}{:15.6f}".format('Chemical potential mu', mu ) )
print ( "{:40}{:15.6f}".format('Activity z',            z  ) )

# Results for cut-and-shifted potential from Thol et al (2015) fitting formula
print('')
print('Lennard-Jones potential cut-and-shifted at 2.5 sigma')
print('')

a_res = a_res_cutshift ( temperature, density )

for (i,j), aij in np.ndenumerate ( a_res ):
    if i+j > 2: # Only interested in some of the results
        continue
    print ( "{:4}{:1d}{:<35d}{:15.6f}".format('Ares', i, j, aij ) )

p  = density * temperature * ( 1.0 + a_res[0,1] )
e  = temperature * ( 1.5 + a_res[1,0] )
cv = 1.5 - a_res[2,0]
cp = 2.5 - a_res[2,0]+(1.0+a_res[0,1]-a_res[1,1])*(1.0+a_res[0,1]-a_res[1,1])/(1.0+2.0*a_res[0,1]+a_res[0,2]) - 1.0
mu = temperature * ( np.log(density) + a_res[0,0] + a_res[0,1] )
z  = density * np.exp ( a_res[0,0] + a_res[0,1] )
print('')
print ( "{:40}{:15.6f}".format('Pressure P',            p  ) )
print ( "{:40}{:15.6f}".format('Energy E/N',            e  ) )
print ( "{:40}{:15.6f}".format('Heat capacity Cv/NkB',  cv ) )
print ( "{:40}{:15.6f}".format('Heat capacity Cp/NkB',  cp ) )
print ( "{:40}{:15.6f}".format('Chemical potential mu', mu ) )
print ( "{:40}{:15.6f}".format('Activity z',            z  ) )
