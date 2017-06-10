#!/usr/bin/env python3
# test_pot_gb.py

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

"""Quadrupole-quadrupole potential and forces."""

import numpy as np
n = 2 # Two-molecule potential
print('test_pot_gb module')
print('Returns potential and force for Gay-Berne')
print(n,'-molecule potential',sep='')

def force ( r, e ):
    """Returns potential pot and numpy arrays f, t of shape (n,3), same as input arguments.

    Demonstrates the calculation of forces from the Gay-Berne potential.
    Written for ease of comparison with the text rather than efficiency!
    """
    from math import isclose

    # Parameters of the Gay-Berne potential                     
    #                                                           
    # The key parameters are                                    
    #    mu, nu  ................ the exponents                 
    #    kappa and kappa' ....... the anisotropies              
    #    kappa is the ratio of intermolecular separations       
    #    sigma_e / sigma_s  i.e. end-to-end / side-by-side      
    #    kappa' is the ratio of well depths                     
    #    epsilon_s / epsilon_e  i.e. side-by-side / end-to-end  
    # The derived parameters are chi and chi'                   
    #    chi = (kappa**2 - 1) / (kappa**2+1)                    
    #    chi' = (z - 1) / (z + 1)
    #    where z = (kappa') ** ( 1 / mu )                       
    #                                                           
    # For convenience kappa' is spelt xappa, chi' is spelt xhi
    # We choose units such that sigma_s = 1.0 and epsilon_0 = 1.0
    # Two of the following three varieties should be commented out

    # Original Gay-Berne-deMiguel potential [J. Chem. Phys, 74, 3316; Mol. Phys. 74, 405 (1991)]
    mu, nu, kappa, xappa = 2, 1, 3.0, 5.0

    # # Luckhurst-Phippen potential [Liq. Cryst., 8, 451 (1990)]
    # mu, nu, kappa, xappa = 1, 2, 3.0, 5.0

    # # Berardi-Zannoni potential [J. Chem. Soc. Faraday Trans., 89, 4069 (1993)]
    # mu, nu, kappa, xappa = 1, 3, 3.0, 5.0

    # Derived parameters
    chi = (kappa**2 - 1.0) / (kappa**2+1.0)
    xhi = (xappa**(1.0/mu) - 1.0) / (xappa**(1.0/mu) + 1.0)

    # Cutoff distance; normally we would use a larger value
    r_cut = 4.0
    
    assert r.shape == (n,3), 'Incorrect shape of r'
    assert e.shape == (n,3), 'Incorrect shape of e'

    # Notation to match appendix
    i = 0
    j = 1

    ei = e[i,:]
    ej = e[j,:]
    assert isclose(np.sum(ei**2),1.0), 'Non-unit vector {} {} {}'.format(*ei)
    assert isclose(np.sum(ej**2),1.0), 'Non-unit vector {} {} {}'.format(*ej)

    rij = r[i,:] - r[j,:]
    rij_mag = np.sqrt( np.sum(rij**2) ) # Magnitude of separation vector
    sij = rij / rij_mag                 # Unit vector
    ci  = np.dot( ei, sij )
    cj  = np.dot( ej, sij )
    cij = np.dot( ei, ej  )
    cp  = ci + cj
    cm  = ci - cj

    # Sigma formula
    cpchi = cp/(1.0+chi*cij)
    cmchi = cm/(1.0-chi*cij)
    sigma = 1.0/np.sqrt(1.0-0.5*chi*(cp*cpchi+cm*cmchi))

    # Epsilon formula
    eps1    = 1.0/np.sqrt(1.0-(chi*cij)**2) # Depends on chi, not xhi
    cpxhi   = cp/(1.0+xhi*cij)
    cmxhi   = cm/(1.0-xhi*cij)
    eps2    = 1.0-0.5*xhi*(cp*cpxhi+cm*cmxhi) # Depends on xhi
    epsilon = (eps1**nu) * (eps2**mu)

    # Potential at rij
    rho      = rij_mag - sigma + 1.0
    rho6     = 1.0 / rho**6
    rho12    = rho6**2
    rhoterm  = 4.0*(rho12 - rho6)                 # Needed for forces and torques
    drhoterm = -24.0 * (2.0 * rho12 - rho6) / rho # Needed for forces and torques
    pot      = epsilon*rhoterm

    # Potential at r_cut
    rho      = r_cut - sigma + 1.0
    rho6     = 1.0 / rho**6
    rho12    = rho6**2
    cutterm  = 4.0*(rho12 - rho6)                 # Needed for cutoff forces and torques
    dcutterm = -24.0 * (2.0 * rho12 - rho6) / rho # Needed for cutoff forces and torques
    pot      = pot - epsilon * cutterm

    # Derivatives of sigma
    prefac    = 0.5*chi*sigma**3
    dsig_dci  = prefac*(cpchi+cmchi)
    dsig_dcj  = prefac*(cpchi-cmchi)
    prefac    = prefac*(0.5*chi)
    dsig_dcij = -prefac*(cpchi**2-cmchi**2)

    # Derivatives of epsilon
    prefac    = -mu*xhi*(eps1**nu)*eps2**(mu-1)
    deps_dci  = prefac*(cpxhi+cmxhi)
    deps_dcj  = prefac*(cpxhi-cmxhi)
    prefac    = prefac*(0.5*xhi)
    deps_dcij = -prefac*(cpxhi**2-cmxhi**2)                           # From derivative of eps2
    deps_dcij = deps_dcij + nu*(chi**2)*(eps1**(nu+2))*(eps2**mu)*cij # From derivative of eps1

    # Derivatives of potential
    dpot_drij = epsilon * drhoterm
    dpot_dci  = rhoterm * deps_dci  - epsilon * drhoterm * dsig_dci
    dpot_dcj  = rhoterm * deps_dcj  - epsilon * drhoterm * dsig_dcj
    dpot_dcij = rhoterm * deps_dcij - epsilon * drhoterm * dsig_dcij

    # Standard formula for forces and torque gradients
    fij = -dpot_drij*sij - dpot_dci*(ei-ci*sij)/rij_mag - dpot_dcj*(ej-cj*sij)/rij_mag
    gi  = dpot_dci*sij  + dpot_dcij*ej
    gj  = dpot_dcj*sij  + dpot_dcij*ei

    # Derivatives of potential at cutoff
    dpot_drij = epsilon * dcutterm
    dpot_dci  = cutterm * deps_dci  - epsilon * dcutterm * dsig_dci
    dpot_dcj  = cutterm * deps_dcj  - epsilon * dcutterm * dsig_dcj
    dpot_dcij = cutterm * deps_dcij - epsilon * dcutterm * dsig_dcij

    # Standard formula for forces and torque gradients (without dpot_drij term)
    fij = fij + dpot_dci*(ei-ci*sij)/rij_mag + dpot_dcj*(ej-cj*sij)/rij_mag
    gi  = gi - ( dpot_dci*sij  + dpot_dcij*ej ) 
    gj  = gj - ( dpot_dcj*sij  + dpot_dcij*ei ) 

    # Final forces and torques
    f      = np.empty_like(r)
    t      = np.empty_like(r)
    f[i,:] = fij
    f[j,:] = -fij
    t[i,:] = -np.cross(ei,gi)
    t[j,:] = -np.cross(ej,gj)

    return pot, f, t
