#!/usr/bin/env python3
# mc_nvt_lj_re.py

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

"""Monte Carlo, NVT ensemble, replica exchange."""

def calc_variables ( ):
    """Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.
    """

    # In this example we simulate using the cut (but not shifted) potential
    # The values of < p_c >, < e_c > and density should be consistent (for this potential)
    # For comparison, long-range corrections are also applied to give
    # estimates of < e_f > and < p_f > for the full (uncut) potential
    # The value of the cut-and-shifted potential is not used, in this example

    import numpy as np
    import math
    from averages_module import msd, VariableType
    from lrc_module      import potential_lrc, pressure_lrc, pressure_delta
    from mc_lj_module    import force_sq
    
    # Preliminary calculations (n,r,total are taken from the calling program)
    vol = box**3                     # Volume
    rho = n / vol                    # Density
    fsq = force_sq ( box, r_cut, r ) # Total squared force

    # Variables of interest, of class VariableType, containing three attributes:
    #   .val: the instantaneous value
    #   .nam: used for headings
    #   .method: indicating averaging method
    # If not set below, .method adopts its default value of avg
    # The .nam and some other attributes need only be defined once, at the start of the program,
    # but for clarity and readability we assign all the values together below

    # Move and exchange acceptance ratios
    m_r = VariableType ( nam = 'Move ratio', val = m_ratio, instant = False )
    x_r = VariableType ( nam = 'Swap ratio', val = x_ratio, instant = False )

    # Internal energy per atom for simulated, cut, potential
    # Ideal gas contribution plus cut (but not shifted) PE divided by N
    e_c = VariableType ( nam = 'E/N cut', val = 1.5*temperature + total.pot/n )

    # Internal energy per atom for full potential with LRC
    # LRC plus ideal gas contribution plus cut (but not shifted) PE divided by N
    e_f = VariableType ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + 1.5*temperature + total.pot/n )

    # Pressure for simulated, cut, potential
    # Delta correction plus ideal gas contribution plus total virial divided by V
    p_c = VariableType ( nam = 'P cut', val = pressure_delta(rho,r_cut) + rho*temperature + total.vir/vol )

    # Pressure for full potential with LRC
    # LRC plus ideal gas contribution plus total virial divided by V
    p_f = VariableType ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*temperature + total.vir/vol )

    # Configurational temperature
    # Total squared force divided by total Laplacian
    t_c = VariableType ( nam = 'T config', val = fsq/total.lap )

    # Heat capacity (full)
    # MSD potential energy divided by temperature and sqrt(N) to make result intensive; LRC does not contribute
    # We add ideal gas contribution, 1.5, afterwards
    c_f = VariableType ( nam = 'Cv/N full', val = total.pot/(temperature*math.sqrt(n)),
                         method = msd, add = 1.5, instant = False )

    # Collect together into a list for averaging
    return [ m_r, x_r, e_c, p_c, e_f, p_f, t_c, c_f ]

# Takes in a configuration of atoms (positions)
# Cubic periodic boundary conditions
# Conducts Monte Carlo at the given temperature
# Uses no special neighbour lists
# Uses MPI to run replica exchange of configurations with neighbouring temperatures
# Assume at most 100 processes numbered 0 to 99

# All processes write to their standard output, but the default in MPI is for all this output
# to be collated (in an undefined order) and written to a single channel. We assume that the program
# will be run with a command-line which includes an option for each process to write to separate files, such as
#          mpiexec -n 8 -output-filename out ./mc_nvt_lj_re.py < mc.inp
# where the standard output files are named out##, the ## part being determined by the process rank.
# If your implementation does not have this option, you should edit the code to explicitly open a file for
# standard output, with a process-rank-dependent name, and associate the output_unit with it.

# Note that configurations are read, saved, and written to files named cnf##.inp etc
# NB a program intended for real-world application would be much more careful about
# closing down all the MPI processes cleanly in the event of an error on any one process.
# However, the documentation for mpi4py currently states that the MPI Init and Finalize are
# automatically called on importing the MPI module and on Python process exit, respectively.
# Therefore, we don't call them explicitly here.

# Reads several variables and options from standard input using JSON format
# Leave input empty "{}" to accept supplied defaults

# Positions r are divided by box length after reading in
# However, input configuration, output configuration, most calculations, and all results
# are given in simulation units defined by the model
# For example, for Lennard-Jones, sigma = 1, epsilon = 1

# Despite the program name, there is nothing here specific to Lennard-Jones
# The model is defined in mc_lj_module

import json
import sys
import os
import numpy as np
import math
from mpi4py           import MPI
from config_io_module import read_cnf_atoms, write_cnf_atoms
from averages_module  import run_begin, run_end, blk_begin, blk_end, blk_add
from maths_module     import random_translate_vector, metropolis
from mc_lj_module     import introduction, conclusion, potential, potential_1, PotentialType

comm  = MPI.COMM_WORLD
m     = comm.Get_rank ( )
nproc = comm.Get_size ( )
assert nproc <= 100, 'Number of processes is too large'
m_tag      = str(m).zfill(2) # Convert rank into character form
inp_tag    = 'inp'
out_tag    = 'out'
sav_tag    = 'sav'
cnf_prefix = 'cnf'+m_tag+'.' # Configuration filename includes process rank and dot
msg1_id, msg2_id, msg3_id, msg4_id = 999, 888, 777, 666 # MPI message identifiers

print('mc_nvt_lj_re')
print('Monte Carlo, constant-NVT, replica exchange')
print('Simulation uses cut (but not shifted) potential')
print( "{:40}{:15d}".format('This is process rank', m)       )
print( "{:40}{:15d}".format('Number of processes is', nproc) )

ioerr = False
if m==0: # Process 0 reads data from standard input
    # Read parameters in JSON format
    try:
        nml = json.load(sys.stdin)
    except json.JSONDecodeError:
        ioerr = True

ioerr = comm.bcast(ioerr,root=0) # Process 0 sends error outcome to all other processes
if ioerr:
    print('Exiting on Invalid JSON format')
    sys.exit()

if m==0: # Process 0 sets default values, check keys and typecheck values
    # Empirical choices for temperature and dr_max give approx 20% swap rate and 35-40% move rate for 256 LJ atoms
    defaults = {"nblock":10, "nstep":1000, "r_cut":2.5,
                "every_temperature":[ 1.00*(1.14)**(i-1) for i in range(nproc) ],
                "every_dr_max":[ 0.15*(1.11)**(i-1) for i in range(nproc) ] }
    for key, val in nml.items():
        if key in defaults:
            if type(val) != type(defaults[key]):
                print(key+" has the wrong type")
                ioerr = True
        else:
            print('Warning', key, 'not in ',list(defaults.keys()))

ioerr = comm.bcast(ioerr,root=0) # Process 0 sends error outcome to all other processes
if ioerr:
    print('Exiting on Invalid JSON format')
    sys.exit()

if m==0: # Process 0 checks temperature and dr_max lists
    if 'every_temperature' in nml:
        if len(nml['every_temperature']) != nproc:
            print('temperature list length must match number of processes')
            ioerr = True
        for temperature in nml['every_temperature']:
            if type(temperature) != float:
                print('temperatures must be floats')
                ioerr = True
    if 'every_dr_max' in nml:
        if len(nml['every_dr_max']) != nproc:
            print('dr_max list length must match number of processes')
            ioerr = True
        for dr_max in nml['every_dr_max']:
            if type(dr_max) != float:
                print('dr_max values must be floats')
                ioerr = True

ioerr = comm.bcast(ioerr,root=0) # Process 0 sends error outcome to all other processes
if ioerr:
    print('Exiting on invalid data')
    sys.exit()

if m==0: # Process 0 sets parameters to input values or defaults
    nblock            = nml["nblock"]            if "nblock"            in nml else defaults["nblock"]
    nstep             = nml["nstep"]             if "nstep"             in nml else defaults["nstep"]
    r_cut             = nml["r_cut"]             if "r_cut"             in nml else defaults["r_cut"]
    every_temperature = nml["every_temperature"] if "every_temperature" in nml else defaults["every_temperature"]
    every_dr_max      = nml["every_dr_max"]      if "every_dr_max"      in nml else defaults["every_dr_max"]
else:
    nblock            = None
    nstep             = None
    r_cut             = None
    every_temperature = None
    every_dr_max      = None

# Process 0 sends run parameters to all other processes
nblock            = comm.bcast(nblock,root=0)
nstep             = comm.bcast(nstep,root=0)
r_cut             = comm.bcast(r_cut,root=0)
every_temperature = comm.bcast(every_temperature,root=0)
every_dr_max      = comm.bcast(every_dr_max,root=0)

every_beta  = [1.0/temperature for temperature in every_temperature] # All the inverse temperatures
temperature = every_temperature[m]    # Temperature for this process
dr_max      = every_dr_max[m]         # Max displacement for this process
beta        = every_beta[m]           # Inverse temperature for this process

introduction()
np.random.seed()

# Write out parameters
print( "{:40}{:15d}  ".format('Number of blocks',          nblock)      )
print( "{:40}{:15d}  ".format('Number of steps per block', nstep)       )
print( "{:40}{:15.6f}".format('Specified temperature',     temperature) )
print( "{:40}{:15.6f}".format('Potential cutoff distance', r_cut)       )
print( "{:40}{:15.6f}".format('Maximum displacement',      dr_max)      )

exists = os.path.isfile(cnf_prefix+inp_tag) # Check that initial configuration file exists
if not exists:
    print(cnf_prefix+inp_tag+' file does not exist')
all_exist = comm.allreduce(exists,op=MPI.LAND)
if not all_exist:
    print('Exiting: one or more configuration files do not exist')
    sys.exit()

# Read in initial configuration
n, box, r = read_cnf_atoms ( cnf_prefix+inp_tag)
print( "{:40}{:15d}  ".format('Number of particles', n) )
print( "{:40}{:15.6f}".format('Box length', box)  )
print( "{:40}{:15.6f}".format('Density', n/box**3)  )
r = r / box           # Convert positions to box units
r = r - np.rint ( r ) # Periodic boundaries

# Initial energy and overlap check
total = potential ( box, r_cut, r )
assert not total.ovr, 'Overlap in initial configuration'

# Initialize arrays for averaging and write column headings
m_ratio = 0.0
x_ratio = 0.0
run_begin ( calc_variables() )

for blk in range(1,nblock+1): # Loop over blocks

    blk_begin()

    for stp in range(nstep): # Loop over steps

        moves = 0

        for i in range(n): # Loop over atoms
            rj = np.delete(r,i,0) # Array of all the other atoms
            partial_old = potential_1 ( r[i,:], box, r_cut, rj ) # Old atom potential, virial etc
            assert not partial_old.ovr, 'Overlap in current configuration'

            ri = random_translate_vector ( dr_max/box, r[i,:] ) # Trial move to new position (in box=1 units)
            ri = ri - np.rint ( ri )                            # Periodic boundary correction
            partial_new = potential_1 ( ri, box, r_cut, rj )    # New atom potential, virial etc

            if not partial_new.ovr: # Test for non-overlapping configuration
                delta = partial_new.pot - partial_old.pot # Use cut (but not shifted) potential
                delta = delta / temperature
                if metropolis ( delta ): # Accept Metropolis test
                    total = total + partial_new - partial_old # Update total values
                    r[i,:] = ri                               # Update position
                    moves = moves + 1                         # Increment move counter

        m_ratio = moves / n

        x_ratio = 0.0
        for updown in range(2): # Loop to look one way then the other
            if m%2 == updown: # Look up, partner is m+1
                if m+1 < nproc: # Ensure partner exists
                    other_beta = every_beta[m+1] # We already know the other beta
                    other_pot = comm.recv(source=m+1,tag=msg1_id) # Receive pot from other process
                    delta = -(beta - other_beta) * ( total.pot - other_pot ) # Delta for Metropolis decision
                    swap  = metropolis ( delta ) # Decision taken on this process
                    comm.send(swap,dest=m+1,tag=msg2_id) # Send decision to other process
                    if swap: # Exchange configurations
                        comm.Sendrecv_replace(r,dest=m+1,sendtag=msg3_id,source=m+1,recvtag=msg4_id)
                        total = potential ( box, r_cut, r ) # Alternatively, we could get this from m+1
                        x_ratio = 1.0
            else: # Look down, partner is m-1
                if m-1 >= 0: # Ensure partner exists
                    comm.send(total.pot,dest=m-1,tag=msg1_id) # Send pot to other process
                    swap=comm.recv(source=m-1,tag=msg2_id) # Receive decision from other process
                    if swap: # Exchange configurations
                        comm.Sendrecv_replace(r,dest=m-1,sendtag=msg4_id,source=m-1,recvtag=msg3_id)
                        total = potential ( box, r_cut, r ) # Alternatively, we could get this from m-1

        blk_add ( calc_variables() )

    blk_end(blk)                                          # Output block averages
    sav_tag = str(blk).zfill(3) if blk<1000 else 'sav'    # Number configuration by block
    write_cnf_atoms ( cnf_prefix+sav_tag, n, box, r*box ) # Save configuration

run_end ( calc_variables() )

total = potential ( box, r_cut, r ) # Double check book-keeping
assert not total.ovr, 'Overlap in final configuration'

write_cnf_atoms ( cnf_prefix+out_tag, n, box, r*box ) # Save configuration
conclusion()
