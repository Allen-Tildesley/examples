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
    
