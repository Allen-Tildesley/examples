#!/usr/bin/env python3
# lrc_module.py

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

"""Long-range and delta corrections for potential energy and pressure."""

# The purpose of this module is simply to gather in one place the common
# functions for long-range and delta corrections for the Lennard-Jones potential.
# If a different potential is used in the simulation, a different file
# (with the same module name) containing different expressions should be substituted

def potential_lrc ( density, r_cut ):
    """Calculates long-range correction for Lennard-Jones potential per atom."""

    import math
    
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut**3 
    return math.pi * ( (8.0/9.0)  * sr3**3  - (8.0/3.0)  * sr3 ) * density

def pressure_lrc ( density, r_cut ):
    """Calculates long-range correction for Lennard-Jones pressure."""

    import math
    
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut**3 
    return math.pi * ( (32.0/9.0) * sr3**3  - (16.0/3.0) * sr3 ) * density**2

def pressure_delta ( density, r_cut ):
    """Calculates correction for Lennard-Jones pressure due to discontinuity in the potential at r_cut."""

    import math
    
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut**3 
    return math.pi * (8.0/3.0) * ( sr3**3  - sr3 ) * density**2

