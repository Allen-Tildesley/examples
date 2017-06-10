#!/usr/bin/env python3
# adjust.py

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

"""Utility program to allow user to change the density or kinetic energy of MC or MD configuration."""

import json
import sys
import numpy as np
from config_io_module import read_cnf_atoms, read_cnf_mols, write_cnf_atoms, write_cnf_mols

# Takes in a configuration of atoms or molecules
# positions, possibly orientations, and optionally, velocities and angular velocities
# Cubic periodic boundary conditions
# Adjusts the density by an amount delta_rho
# and kinetic energy per particle by an amount delta_kin

# Input configuration and output configuration are assumed to be in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1
# There is nothing here specific to Lennard-Jones
# We assume unit mass and adjust only the translational kinetic energy

filename = 'cnf.inp'
atom, linear, nonlinear, chain = 0, 1, 2, 3 # User options
tol = 1.e-9

# Read parameters in JSON format
try:
    nml = json.load(sys.stdin)
except json.JSONDecodeError:
    print('Exiting on Invalid JSON format')
    sys.exit()

# Set default values, check keys and typecheck values
# If neither delta_rho nor delta_kin is changed, program will just write out information
# Options for molecules are 'atom', 'chain', 'linear', 'nonlinear'
defaults = {"delta_rho":0.0, "delta_kin":0.0, "velocities":False, "molecules":"atom"}
for key, val in nml.items():
    if key in defaults:
        assert type(val) == type(defaults[key]), key+" has the wrong type"
    else:
        print('Warning', key, 'not in ',list(defaults.keys()))

# Set parameters to input values or defaults
delta_rho  = nml["delta_rho"]  if "delta_rho"  in nml else defaults["delta_rho"]
delta_kin  = nml["delta_kin"]  if "delta_kin"  in nml else defaults["delta_kin"]
velocities = nml["velocities"] if "velocities" in nml else defaults["velocities"]
molecules  = nml["molecules"]  if "molecules"  in nml else defaults["molecules"]

molecules = molecules.lower()
assert  ( "atom" in molecules or "linear" in molecules or
          "nonlin" in molecules or "chain" in molecules    ), 'Unrecognized molecules option'

if "nonlin" in molecules:
    molecule_option = nonlinear
    print('Nonlinear molecules')
elif "linear" in molecules:
    molecule_option = linear
    print('Linear molecules')
elif "atom" in molecules:
    molecule_option = atom
    print('Atoms')
else:
    molecule_option = chain
    print('Atoms in a chain')

assert velocities or np.fabs(delta_kin) < tol, 'No kinetic energy change possible'
assert molecule_option != chain or np.fabs(delta_rho) < tol, 'No density change possible'

atomic = molecule_option==atom or molecule_option==chain
quaternions = molecule_option==nonlinear

if atomic:
    if velocities:
        n, box, r, v = read_cnf_atoms ( filename, with_v=True )
    else:
        n, box, r = read_cnf_atoms ( filename, with_v=False )
else:
    if velocities:
        n, box, r, e, v, w = read_cnf_mols ( filename, with_v=True, quaternions=quaternions )
    else:
        n, box, r, e = read_cnf_mols ( filename, with_v=False, quaternions=quaternions )

if molecule_option == chain:
    print("{:40}{:15.6f}".format('Molecular bond length',box) ) # box plays role of bond
else:
    rho = n / box**3
    print("{:40}{:15.6f}".format('Simulation box length',box) )
    print("{:40}{:15.6f}".format('Density',rho) )
        
if velocities:
    kin = 0.5*np.sum(v**2) / n
    print("{:40}{:15.6f}".format('Kinetic energy',kin) )

if np.fabs(delta_rho) < tol and np.fabs(delta_kin) < tol:
    print('No changes requested')
else:
    if np.fabs(delta_rho) > tol:
        assert rho+delta_rho > 0.0, 'New requested density would be negative'
        scale = ( rho / (rho+delta_rho) )**(1.0/3.0)
        box   = box * scale
        r     = r * scale
        rho   = n / box**3
        print("{:40}{:15.6f}".format('Density',rho) )
    if np.fabs(delta_kin) > tol:
        assert kin+delta_kin > 0.0, 'New requested kinetic energy would be negative'
        scale  = np.sqrt ( (kin+delta_kin) / kin )
        v      = v * scale
        kin    = 0.5*np.sum(v**2) / n
        print("{:40}{:15.6f}".format('New kinetic energy',kin) )

    if atomic:
        if velocities:
            write_cnf_atoms ( filename, n, box, r, v )
        else:
            write_cnf_atoms ( filename, n, box, r )
    else:
        if velocities:
            write_cnf_mols ( filename, n, box, r, e, v, w )
        else:
            write_cnf_mols ( filename, n, box, r, e )
