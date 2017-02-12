#!/usr/bin/env python3
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

