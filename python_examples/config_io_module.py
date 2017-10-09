#!/usr/bin/env python3
# config_io_module.py

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

"""Routines for atomic/molecular configuration I/O."""

def read_cnf_atoms ( filename, with_v=False ):
    """Read in atomic configuration."""
    
    import numpy as np
    
    with open(filename,"r") as f:
        n=int(f.readline()) # Number of atoms
        box=float(f.readline()) # Simulation box length (assumed cubic)
    rv=np.loadtxt(filename,skiprows=2) # The rest of the file should be uniformly formatted
    rows, cols = rv.shape
    assert rows == n, "{:d}{}{:d}".format(rows,' rows not equal to ',n)
    assert cols >= 3, "{:d}{}".format(cols,' cols less than 3')
    r = rv[:,0:3].astype(np.float_) # Coordinate array
    if with_v:
        assert cols >= 6, "{:d}{}".format(cols,' cols less than 6')
        v = rv[:,3:6].astype(np.float_) # Velocity array
        return n, box, r, v
    else:
        return n, box, r

def read_cnf_mols ( filename, with_v=False, quaternions=False ):
    """Read in molecular configuration."""
    
    import numpy as np
    
    with open(filename,"r") as f:
        n=int(f.readline()) # Number of atoms
        box=float(f.readline()) # Simulation box length (assumed cubic)
    revw=np.loadtxt(filename,skiprows=2) # The rest of the file should be uniformly formatted
    rows, cols = revw.shape
    assert rows == n, "{:d}{}{:d}".format(rows,' rows not equal to ',n)
    cols_re = 7 if quaternions else 6
    assert cols >= cols_re, "{:d}{}{:d}".format(cols,' cols less than ',cols_re)
    r = revw[:,0:3].astype(np.float_) # Coordinate array
    e = revw[:,3:cols_re].astype(np.float_) # Orientation array
    if with_v:
        assert cols >= cols_re+6, "{:d}{}{:d}".format(cols,' cols less than',cols_re+6)
        v = revw[:,cols_re  :cols_re+3].astype(np.float_) # Velocity array
        w = revw[:,cols_re+3:cols_re+6].astype(np.float_) # Angular velocity/momentum array
        return n, box, r, e, v, w
    else:
        return n, box, r, e

def write_cnf_atoms ( filename, nn, box, *args ):
    """Write out atomic configuration."""

    import numpy as np

    nargs=len(args)
    assert nargs==1 or nargs==2, "{}{:d}".format('Wrong number of arguments ',nargs)

    my_header="{:15d}\n{:15.8f}".format(nn,box)

    n, d = args[0].shape
    assert n==nn, "r shape mismatch {:d}{:d}".format(n,nn)
    assert d==3,  "r shape mismatch {:d}".format(d)

    if nargs==1: # Positions only
        r=args[0] # Just r
        np.savetxt(filename,r,header=my_header,fmt='%15.10f',comments='')

    else: # Positions and velocities
        n, d = args[1].shape
        assert n==nn, "v shape mismatch {:d}{:d}".format(n,nn)
        assert d==3, "v shape mismatch {:d}".format(d)
        rv=np.concatenate((args[0],args[1]),axis=1) # Both r and v
        np.savetxt(filename,rv,header=my_header,fmt='%15.10f',comments='')

def write_cnf_mols ( filename, nn, box, *args ):
    """Write out molecular configuration."""

    import numpy as np

    nargs=len(args)
    assert nargs==2 or nargs==4, "{}{:d}".format('Wrong number of arguments ',nargs)

    my_header="{:15d}\n{:15.8f}".format(nn,box)

    n, d = args[0].shape # Positions
    assert n==nn, "r shape mismatch {:d}{:d}".format(n,nn)
    assert d==3,  "r shape mismatch {:d}".format(d)

    n, d = args[1].shape # Orientations: vectors or quaternions
    assert n==nn, "e shape mismatch {:d}{:d}".format(n,nn)
    assert d==3 or d==4, "e shape mismatch {:d}".format(d)

    if nargs==2: # Positions and orientations
        re=np.concatenate((args[0],args[1]),axis=1) # Just r and e
        np.savetxt(filename,re,header=my_header,fmt='%15.10f',comments='')

    else: # Positions, orientations, velocities and angular velocities/momenta
        n, d = args[2].shape
        assert n==nn, "v shape mismatch {:d}{:d}".format(n,nn)
        assert d==3, "v shape mismatch {:d}".format(d)
        n, d = args[3].shape
        assert n==nn, "w shape mismatch {:d}{:d}".format(n,nn)
        assert d==3, "w shape mismatch {:d}".format(d)
        revw=np.concatenate((args[0],args[1],args[2],args[3]),axis=1) # All of r, e, v and w
        np.savetxt(filename,revw,header=my_header,fmt='%15.10f',comments='')
