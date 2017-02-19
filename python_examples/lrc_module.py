#!/usr/bin/env python3
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

