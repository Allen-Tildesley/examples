#!/usr/bin/env python3
# eos_lj_module.py

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

"""Routines for Lennard-Jones fitted equations of state."""
  
# The routines in this module use the fitting function described and parametrized in
# M Thol, G Rutkai, R Span, J Vrabec, R Lustig, Int J Thermophys 36, 25 (2015)
# M Thol, G Rutkai, A Koester, R Lustig, R Span, J Vrabec, J Phys Chem Ref Data 45, 023101 (2016)
# Those authors also supply C++ codes (in the supplementary information of those papers)
# They are NOT responsible for this Python code, which was written independently by Michael P Allen
# A similar notation, consistent with the papers, is retained for clarity.


def power ( tau, delta, c ):
    """The power basis function."""

    import numpy as np
    
    # f[0,0] is n*(tau**t)*(delta**d)
    # f[i,:] is differentiated i times with respect to tau, and then multiplied by tau**i
    # f[:,j] is differentiated j times with respect to delta, and then multiplied by delta**j

    f = np.full ( (3,3), c['n'] * (tau**c['t']) * (delta**c['d']), dtype=np.float_ )
    f[1,:] = f[1,:] * c['t']
    f[2,:] = f[2,:] * c['t'] * ( c['t'] - 1.0 )
    f[:,1] = f[:,1] * c['d']
    f[:,2] = f[:,2] * c['d'] * ( c['d'] - 1.0 )

    return f

def expon ( tau, delta, c ):
    """The exponential basis function."""

    import numpy as np
    
    # f[0,0] is n*(tau**t)*(delta**d)*exp(-delta**l)
    # f[i,:] is differentiated i times with respect to tau, and then multiplied by tau**i
    # f[:,j] is differentiated j times with respect to delta, and then multiplied by delta**j

    f = np.full ( (3,3), c['n'] * (tau**c['t']) * (delta**c['d']) * np.exp(-delta**c['l']), dtype=np.float_ )
    f[1,:] = f[1,:] * c['t']
    f[2,:] = f[2,:] * c['t'] * ( c['t'] - 1.0 )
    f[:,1] = f[:,1] * (c['d']-c['l']*delta**c['l'])
    f[:,2] = f[:,2] * ( (c['d']-c['l']*delta**c['l']) * (c['d']-1.0-c['l']*delta**c['l']) - (c['l']**2)*delta**c['l'] )

    return f

def gauss ( tau, delta, c ):
    """The gaussian basis function."""

    import numpy as np

    # f[0,0] is n*(tau**t)*exp(-beta*(tau-gamma)**2)*(delta**d)*exp(-eta*(delta-epsilon)**2)
    # f[i,:] is differentiated i times with respect to tau, and then multiplied by tau**i
    # f[:,j] is differentiated j times with respect to delta, and then multiplied by delta**j

    f = np.full ( (3,3), c['n']*(tau**c['t'])*np.exp(-c['beta']*(tau-c['gamma'])**2)
                      *(delta**c['d'])*np.exp(-c['eta']*(delta-c['epsilon'])**2), dtype=np.float_)
    f[1,:] = f[1,:] * ( c['t'] - 2.0*c['beta']*tau*(tau-c['gamma']) )
    f[2,:] = f[2,:] * ( ( c['t'] - 2.0*c['beta']*tau*(tau-c['gamma']) )**2 - c['t'] - 2*c['beta']*tau**2 )
    f[:,1] = f[:,1] * ( c['d'] - 2.0*c['eta']*delta*(delta-c['epsilon']) )
    f[:,2] = f[:,2] * ( ( c['d'] - 2.0*c['eta']*delta*(delta-c['epsilon']) )**2 - c['d'] - 2*c['eta']*delta**2 )

    return f

def a_res_full ( temp, rho ):
    """Reduced residual free energy and scaled derivatives for the full Lennard-Jones potential."""

    import numpy as np

    # Given temperature and density in LJ units.
    # In a[i,j], index i refers to the tau-derivative and index j to the delta-derivative
    # The derivatives are multiplied by the corresponding powers of tau and delta
    # a[i,:] is differentiated i times with respect to tau, and then multiplied by tau**i
    # a[:,j] is differentiated j times with respect to delta, and then multiplied by delta**j
    

    temp_crit = 1.32 # Critical temperature
    rho_crit  = 0.31 # Critical density

    tau   = temp_crit / temp # Reduced inverse temperature
    delta = rho / rho_crit   # Reduced density

    a = np.zeros ( (3,3), dtype=np.float_ )

    # Coefficients taken from Table 2 of 
    # M Thol, G Rutkai, A Koester, R Lustig, R Span, J Vrabec, J Phys Chem Ref Data 45, 023101 (2016)

    cp = np.empty ( 6, dtype = [ ('n',np.float_), ('t',np.float_), ('d',np.float_) ] )
    cp['n'] = [ 0.005208073, 2.186252, -2.161016, 1.452700, -2.041792, 0.18695286 ]
    cp['t'] = [ 1.000, 0.320, 0.505, 0.672, 0.843, 0.898 ]
    cp['d'] = [ 4.0,   1.0,   1.0,   2.0,   2.0,   3.0   ]
    for c in cp:
        a = a + power ( tau, delta, c )

    ce = np.empty ( 6, dtype = [ ('n',np.float_), ('t',np.float_), ('d',np.float_), ('l',np.float_) ] )
    ce['n'] = [ -0.090988445, -0.49745610, 0.10901431, -0.80055922, -0.56883900, -0.62086250 ]
    ce['t'] = [ 1.294, 2.590, 1.786, 2.770, 1.786, 1.205 ]
    ce['d'] = [ 5.0,   2.0,   2.0,   3.0,   1.0,   1.0   ]
    ce['l'] = [ 1.0,   2.0,   1.0,   2.0,   2.0,   1.0   ]
    for c in ce:
        a = a + expon ( tau, delta, c )

    cg = np.empty ( 11, dtype = [ ('n',np.float_), ('t',np.float_), ('d',np.float_), ('eta',np.float_),
                                  ('beta',np.float_), ('gamma',np.float_), ('epsilon',np.float_)] )
    cg['n']       = [ -1.4667177,  1.8914690,  -0.13837010, -0.38696450, 0.12657020, 0.6057810,
                       1.1791890, -0.47732679, -9.9218575,  -0.57479320, 0.0037729230 ]
    cg['t']       = [ 2.830,  2.548, 4.650, 1.385, 1.460, 1.351, 0.660, 1.496,   1.830, 1.616, 4.970 ]
    cg['d']       = [ 1.0,    1.0,   2.0,   3.0,   3.0,   2.0,   1.0,   2.0,     3.0,   1.0,   1.0   ]
    cg['eta']     = [ 2.067,  1.522, 8.82,  1.722, 0.679, 1.883, 3.925, 2.461,  28.2,   0.753, 0.82  ]
    cg['beta']    = [ 0.625,  0.638, 3.91,  0.156, 0.157, 0.153, 1.16,  1.73,  383.0,   0.112, 0.119 ]
    cg['gamma']   = [ 0.71,   0.86,  1.94,  1.48,  1.49,  1.945, 3.02,  1.11,    1.17,  1.33,  0.24  ]
    cg['epsilon'] = [ 0.2053, 0.409, 0.6,   1.203, 1.829, 1.397, 1.39,  0.539,   0.934, 2.369, 2.43  ]
    for c in cg:
        a = a + gauss ( tau, delta, c )

    return a

def a_res_cutshift ( temp, rho ):
    """Reduced residual free energy and scaled derivatives for Lennard-Jones potential cut-and-shifted at 2.5 sigma."""

    import numpy as np
    
    # Given temperature and density in LJ units.
    # This routine is for the 
    # In a[i,j], index i refers to the tau-derivative and index j to the delta-derivative
    # The derivatives are multiplied by the corresponding powers of tau and delta
    # a[i,:] is differentiated i times with respect to tau, and then multiplied by tau**i
    # a[:,j] is differentiated j times with respect to delta, and then multiplied by delta**j

    temp_crit = 1.086 # Critical temperature
    rho_crit  = 0.319 # Critical density

    tau   = temp_crit / temp # Reduced inverse temperature
    delta = rho / rho_crit   # Reduced density

    a = np.zeros ( (3,3), dtype=np.float_ )

    # Coefficients taken from Table 1 of 
    # M Thol, G Rutkai, R Span, J Vrabec, R Lustig, Int J Thermophys 36, 25 (2015) 

    cp = np.empty ( 6, dtype = [ ('n',np.float_), ('t',np.float_), ('d',np.float_) ] )
    cp['n'] = [ 0.015606084, 1.7917527, -1.9613228, 1.3045604, -1.8117673, 0.15483997 ]
    cp['t'] = [ 1.000, 0.304, 0.583, 0.662, 0.870, 0.870 ]
    cp['d'] = [ 4.0,   1.0,   1.0,   2.0,   2.0,   3.0   ]
    for c in cp:
        a = a + power ( tau, delta, c )

    ce = np.empty ( 6, dtype = [ ('n',np.float_), ('t',np.float_), ('d',np.float_), ('l',np.float_) ] )
    ce['n'] = [ -0.094885204, -0.20092412,  0.11639644, -0.50607364, -0.58422807, -0.47510982 ]
    ce['t'] = [  1.250, 3.000, 1.700, 2.400, 1.960, 1.286 ]
    ce['d'] = [  5.0,   2.0,   2.0,   3.0,   1.0,   1.0   ]
    ce['l'] = [  1.0,   2.0,   1.0,   2.0,   2.0,   1.0   ]
    for c in ce:
        a = a + expon ( tau, delta, c )

    cg = np.empty ( 9, dtype = [ ('n',np.float_), ('t',np.float_), ('d',np.float_), ('eta',np.float_),
                                 ('beta',np.float_), ('gamma',np.float_), ('epsilon',np.float_)] )
    cg['n']       = [  0.0094333106, 0.30444628,  -0.0010820946, -0.099693391, 0.0091193522,
                       0.12970543,   0.023036030, -0.082671073,  -2.2497821 ]
    cg['t']       = [  3.600, 2.080, 5.240, 0.960, 1.360, 1.655, 0.900, 0.860,  3.950 ]
    cg['d']       = [  1.0,   1.0,   2.0,   3.0,   3.0,   2.0,   1.0,   2.0,    3.0   ]
    cg['eta']     = [  4.70,  1.92,  2.70,  1.49,  0.65,  1.73,  3.70,  1.90,  13.2   ]
    cg['beta']    = [ 20.0,   0.77,  0.5,   0.8,   0.4,   0.43,  8.0,   3.3,  114.0   ]
    cg['gamma']   = [  1.0,   0.5,   0.8,   1.5,   0.7,   1.6,   1.3,   0.6,    1.3   ]
    cg['epsilon'] = [  0.55,  0.7,   2.0,   1.14,  1.2,   1.31,  1.14,  0.53,   0.96  ]
    for c in cg:
        a = a + gauss ( tau, delta, c )

    return a
