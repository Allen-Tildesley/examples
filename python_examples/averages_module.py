#!/usr/bin/env python3
# averages_module.py

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

"""Calculation of run averages with output to stdout."""

# Options for averaging methods
avg = 0
msd = 1
cke = 2

# Internal variables
col_width = 15              # Must be large enough to allow sensible format
nam_width = 2*col_width+1   # At most two column widths plus spacer
colf_fmt  = ' {:15.6f}'     # Format for floats; we assume that 6 dp will be sufficient
cole_fmt  = ' {:15.4e}'     # Alternative format for floats, suitable for small numbers
head_fmt  = ' {:>15}'       # Format for heading strings
col1a_fmt = '{:>15}'        # Format for column 1 strings
col1i_fmt = '{:15d}'        # Format for column 1 integers
sngl_fmt  = '{:40}{:15.6f}' # Format for single line output

class VariableType:
    """Class encapsulating the essential information for simulation averages."""

    def __init__(self, nam, val, method=avg, add=0.0, e_format=False, instant=True):
        self.nam      = nam
        self.val      = val
        self.method   = method
        self.add      = add
        self.e_format = e_format
        self.instant  = instant

def time_stamp():
    """Function to print date, time, and cpu time information."""

    import time

    print("{:45}{}".format("Date:",time.strftime("%Y/%m/%d")))
    print("{:47}{}".format("Time:",time.strftime("%H:%M:%S")))
    print("{:40}{:15.6f}".format("CPU time:",time.process_time()))
    
def run_begin ( variables ):
    """Set up averaging variables based on supplied list of names & other attributes."""
    
    import numpy as np
    global n_avg, headings, subheads, line_fmt, line_width, method, add, run_nrm, run_avg, run_err

    need_header=True
    for variable in variables:
        if variable.instant:
            if need_header:
                print('Initial values')
                need_header=False
            print(sngl_fmt.format(variable.nam,variable.val))

    n_avg = len(variables)

    # First column plus a column for each variable; allow one space between columns
    line_width = col_width + n_avg * ( col_width + 1 )

    # Store variable names in module variables
    # Attempt to split name string at first space
    # Build up format string for line of averages
    headings=[]
    subheads=[]
    line_fmt=[]
    for variable in variables:
        parts=variable.nam.strip().split(' ',maxsplit=1)
        headings.append(parts[0])
        if len(parts)>1:
            subheads.append(parts[1])
        else:
            subheads.append(' ')
        if variable.e_format:
            line_fmt.append(cole_fmt)
        else:
            line_fmt.append(colf_fmt)

    line_fmt=''.join(line_fmt) # Change line format list into string
    headings_fmt = "{:>15}"+head_fmt*n_avg

    # Store method options and add-constants in module variables
    method=np.array([variable.method for variable in variables],dtype=np.int_)
    add   =np.array([variable.add    for variable in variables],dtype=np.float_)

    # Zero averages and error accumulators at start of run
    run_nrm = 0.0
    run_avg = np.zeros(n_avg,dtype=np.float_)
    run_err = np.zeros(n_avg,dtype=np.float_)

    # Write headings
    print()
    print('Run begins')
    time_stamp()
    print()
    print('='*line_width)
    print(headings_fmt.format('Block',*headings))
    print(headings_fmt.format('     ',*subheads))
    print('-'*line_width)
    
def blk_begin():
    """Zero average variables at start of each block."""
    
    import numpy as np
    global blk_nrm, blk_avg, blk_msd

    blk_nrm = 0.0
    blk_avg = np.zeros(n_avg,dtype=np.float_)
    blk_msd = np.zeros(n_avg,dtype=np.float_)

def blk_add ( variables ):
    """Increment block-average variables."""
    
    import numpy as np
    global blk_nrm, blk_avg, blk_msd

    assert len(variables)==n_avg, 'Mismatched variable arrays'
    
    values = np.array([variable.val for variable in variables],dtype=np.float_)
    blk_avg = blk_avg + values
    blk_msd = blk_msd + values**2
    blk_nrm = blk_nrm + 1.0

def blk_end ( blk ):
    """Write out averages at end of every block."""
    
    import numpy as np
    global run_avg, run_err, run_nrm, blk_avg, blk_msd, blk_rnm
    
    assert blk_nrm>0.5, 'Block accumulation error'

    blk_avg = blk_avg / blk_nrm # Normalize block averages
    blk_msd = blk_msd / blk_nrm # Normalize block averages of squared variables

    # Replace blk_avg by mean-squared deviations plus optional constant where required
    mask = np.logical_or(method == msd,method == cke)
    blk_avg = np.where(mask,add+blk_msd-blk_avg**2,blk_avg)
    if np.any ( method == cke ):
        cke_calc() # Call special routine for Cv from KE fluctuations

    run_avg = run_avg + blk_avg    # Increment run averages
    run_err = run_err + blk_avg**2 # Increment error accumulators
    run_nrm = run_nrm + 1.0        # Increment run normalizer

    # Write out block averages
    print((col1i_fmt+line_fmt).format(blk,*blk_avg))

def run_end ( variables ):
    """Write out averages and error estimates at end of run."""
    
    import numpy as np
    global run_avg, run_err
    
    assert run_nrm>0.5, 'Run accumulation error'

    # NB, these are the crudest possible error estimates, based on the wholly unjustified
    # assumption that the blocks are statistically independent
    # For a discussion of errors, see Chapter 8 and the error_calc.py example

    run_avg = run_avg / run_nrm    # Normalize run averages
    run_err = run_err / run_nrm    # Normalize error estimates
    run_err = run_err - run_avg**2 # Compute fluctuations of block averages

    # Normalize and get estimated errors guarding against roundoff
    run_err = np.where ( run_err > 0.0, np.sqrt(run_err/run_nrm), 0.0 )

    print('-'*line_width)
    print((col1a_fmt+line_fmt).format('Run averages',*run_avg))
    print((col1a_fmt+line_fmt).format('Run errors',  *run_err))
    print('='*line_width)
    print()
    print('Run ends')
    time_stamp()
    print()

    need_header=True
    for variable in variables:
        if variable.instant:
            if need_header:
                print('Final values')
                need_header=False
            print(sngl_fmt.format(variable.nam,variable.val))

def cke_calc():
    """Special fluctuation formula for microcanonical heat capacity."""
    
    import numpy as np
    global blk_avg
    
    # Locate variable corresponding to kinetic temperature

    found=False
    for head,subhead,temperature in zip(headings,subheads,blk_avg):
        if 'T' in head+subhead and 'kin' in head+subhead:
            found=True
            break

    assert found, 'Could not find T kin variable'

    # Apply special fluctuation formula for microcanonical ensemble heat capacity
    # blk_avg[i] should contain mean-squared total KE, divided by N

    blk_avg = np.where ( method==cke, 9.0 / ( 6.0 - 4.0 * blk_avg / temperature**2 ), blk_avg )
