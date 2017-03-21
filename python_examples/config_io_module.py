#!/usr/bin/env python3
# config_io_module.py

#------------------------------------------------------------------------------------------------#
# This software was written in 2016/17                                                           #
# by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        #
# and Dominic J. Tildesley <dominic.tildesley@epfl.ch> ("the authors"),                          #
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

"""Routines for atomic/molecular configuration I/O."""

def read_cnf_atoms ( filename, with_v=False ):
    """Read in atomic configuration."""
    
    import numpy as np
    import sys
    with open(filename,"r") as f:
        n=int(f.readline()) # Number of atoms
        box=float(f.readline()) # Simulation box length (assumed cubic)
    rv=np.loadtxt(filename,skiprows=2) # The rest of the file should be uniformly formatted
    rows, cols = rv.shape
    if rows != n:
        print("{}{:5}{}{:5}".format('Number of rows',rows,' not equal to',n))
        sys.exit()
    if cols < 3:
        print("{}{:5}{}".format('Number of cols',cols,' less than 3'))
        sys.exit()
    r = rv[:,0:3].astype(np.float_) # Coordinate array
    if with_v:
        if cols < 6:
            print("{}{:5}{}".format('Number of cols',cols,' less than 6'))
            sys.exit()
        v = rv[:,3:6].astype(np.float_) # Velocity array
        return n, box, r, v
    else:
        return n, box, r

def write_cnf_atoms ( filename, n, box, *args ):
    """Write out atomic configuration."""

    import numpy as np
    import sys

    nargs=len(args)
    assert nargs>0, "No array arguments supplied in write_cnf_atoms"
    assert args[0].shape[0]==n, "r shape mismatch {:5d}{:5d}".format(n,args[0].shape[0])
    assert args[0].shape[1]==3, "r shape mismatch {:5d}{:5d}".format(3,args[0].shape[1])
    if nargs==1:
        rv=args[0] # Just r
    elif nargs==2:
        assert args[1].shape[0]==n, "v shape mismatch {:5d}{:5d}".format(n,args[1].shape[0])
        assert args[1].shape[1]==3, "v shape mismatch {:5d}{:5d}".format(3,args[1].shape[1])
        rv=np.concatenate((args[0],args[1]),axis=1) # Both r and v
    else:
        print('Too many array arguments in write_cnf_atoms')
        sys.exit()

    my_header="{:15d}\n{:15.8f}".format(n,box)
    np.savetxt(filename,rv,header=my_header,fmt='%15.10f',comments='')
    
