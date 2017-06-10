! test_pot_gb.f90
! Pair potential, Gay-Berne
MODULE test_pot_module

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             !
  ! to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     !
  ! published by Oxford University Press ("the publishers").                                       !
  !                                                                                                !
  ! LICENCE                                                                                        !
  ! Creative Commons CC0 Public Domain Dedication.                                                 !
  ! To the extent possible under law, the authors have dedicated all copyright and related         !
  ! and neighboring rights to this software to the PUBLIC domain worldwide.                        !
  ! This software is distributed without any warranty.                                             !
  ! You should have received a copy of the CC0 Public Domain Dedication along with this software.  !
  ! If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               !
  !                                                                                                !
  ! DISCLAIMER                                                                                     !
  ! The authors and publishers make no warranties about the software, and disclaim liability       !
  ! for all uses of the software, to the fullest extent permitted by applicable law.               !
  ! The authors and publishers do not recommend use of this software for any purpose.              !
  ! It is made freely available, solely to clarify points made in the text. When using or citing   !
  ! the software, you should not imply endorsement by the authors or publishers.                   !
  !------------------------------------------------------------------------------------------------!

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routine
  PUBLIC :: force

  ! Public data
  INTEGER, PARAMETER, PUBLIC :: n = 2 ! Pair potential

CONTAINS

  SUBROUTINE force  ( r, e, pot, f, t )
    USE maths_module, ONLY : cross_product
    IMPLICIT NONE
    REAL, DIMENSION(:,:),           INTENT(in)  :: r, e
    REAL,                           INTENT(out) :: pot
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out) :: f, t

    ! Parameters of the Gay-Berne potential                     
    !                                                           
    ! The key parameters are                                    
    !    mu, nu  ................ the exponents                 
    !    kappa and kappa' ....... the anisotropies              
    !    kappa is the ratio of intermolecular separations       
    !    sigma_e / sigma_s  i.e. end-to-end / side-by-side      
    !    kappa' is the ratio of well depths                     
    !    epsilon_s / epsilon_e  i.e. side-by-side / end-to-end  
    ! The derived parameters are chi and chi'                   
    !    chi = (kappa**2 - 1) / (kappa**2+1)                    
    !    chi' = (z - 1) / (z + 1)
    !    where z = (kappa') ** ( 1 / mu )                       
    !                                                           
    ! For convenience kappa' is spelt xappa, chi' is spelt xhi
    ! We choose units such that sigma_s = 1.0 and epsilon_0 = 1.0
    ! Two of the following three varieties should be commented out

    ! Original Gay-Berne-deMiguel potential [J. Chem. Phys, 74, 3316; Mol. Phys. 74, 405 (1991)]
    INTEGER, PARAMETER :: mu = 2, nu = 1
    REAL,    PARAMETER :: kappa = 3.0, xappa = 5.0

!!$    ! Luckhurst-Phippen potential [Liq. Cryst., 8, 451 (1990)]
!!$    INTEGER, PARAMETER :: mu = 1, nu = 2
!!$    REAL,    PARAMETER :: kappa = 3.0, xappa = 5.0

!!$    ! Berardi-Zannoni potential [J. Chem. Soc. Faraday Trans., 89, 4069 (1993)]
!!$    INTEGER, PARAMETER :: mu = 1, nu = 3
!!$    REAL,    PARAMETER :: kappa = 3.0, xappa = 5.0

    REAL, PARAMETER :: chi = (kappa**2 - 1.0) / (kappa**2+1.0)
    REAL, PARAMETER :: xhi = (xappa**(1.0/mu) - 1.0) / (xappa**(1.0/mu) + 1.0)

    REAL, DIMENSION(3) :: rij, sij, fij, gi, gj
    REAL               :: rij_mag, rij_sq, ci, cj, cij, cp, cm, prefac
    REAL               :: cpchi, cmchi, sigma
    REAL               :: cpxhi, cmxhi, eps1, eps2, epsilon
    REAL               :: dsig_dci, dsig_dcj, dsig_dcij
    REAL               :: deps_dci, deps_dcj, deps_dcij
    REAL               :: dpot_dci, dpot_dcj, dpot_dcij, dpot_drij
    REAL               :: rho, rho6, rho12, rhoterm, drhoterm, cutterm, dcutterm
    REAL,    PARAMETER :: r_cut = 4.0   ! Normally would use a larger value of r_cut
    REAL,    PARAMETER :: tol = 1.e-6
    INTEGER, PARAMETER :: i = 1, j = 2 ! Notation to match appendix

    ! Routine to demonstrate the calculation of forces and torques from the
    ! Gay-Berne potential, including the spherical cutoff contribution.
    ! Written for ease of comparison with the text, rather than efficiency!

    ! Check dimensions to be sure
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'r shape error', SHAPE(r), 3, n
       STOP 'Error in test_pot_gb'
    END IF
    IF ( ANY ( SHAPE(e) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'e shape error', SHAPE(e), 3, n
       STOP 'Error in test_pot_gb'
    END IF

    ! Check unit vectors
    IF ( ABS(SUM(e(:,i)**2)-1.0) > tol .OR. ABS(SUM(e(:,j)**2)-1.0) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a)'      ) 'Warning, non-unit vectors'
       WRITE ( unit=error_unit, fmt='(4f10.5)' ) e(:,i), SUM(e(:,i)**2)
       WRITE ( unit=error_unit, fmt='(4f10.5)' ) e(:,j), SUM(e(:,j)**2)
    END IF

    rij     = r(:,i) - r(:,j)
    rij_sq  = SUM ( rij**2 )
    rij_mag = SQRT(rij_sq)  ! Magnitude of separation
    sij     = rij / rij_mag ! Unit vector along rij

    ! Orientation-dependent terms
    ci  = DOT_PRODUCT ( sij, e(:,i) )
    cj  = DOT_PRODUCT ( sij, e(:,j) ) 
    cij = DOT_PRODUCT ( e(:,i), e(:,j) ) 
    cp  = ci + cj
    cm  = ci - cj

    ! Sigma formula
    cpchi = cp/(1.0+chi*cij)
    cmchi = cm/(1.0-chi*cij)
    sigma = 1.0/SQRT(1.0-0.5*chi*(cp*cpchi+cm*cmchi))

    ! Epsilon formula
    eps1    = 1.0/SQRT(1.0-(chi*cij)**2) ! Depends on chi, not xhi
    cpxhi   = cp/(1.0+xhi*cij)
    cmxhi   = cm/(1.0-xhi*cij)
    eps2    = 1.0-0.5*xhi*(cp*cpxhi+cm*cmxhi) ! Depends on xhi
    epsilon = (eps1**nu) * (eps2**mu)

    ! Potential at rij
    rho      = rij_mag - sigma + 1.0
    rho6     = 1.0 / rho**6
    rho12    = rho6**2
    rhoterm  = 4.0*(rho12 - rho6)                 ! Needed for forces and torques
    drhoterm = -24.0 * (2.0 * rho12 - rho6) / rho ! Needed for forces and torques
    pot      = epsilon*rhoterm

    ! Potential at r_cut
    rho      = r_cut - sigma + 1.0
    rho6     = 1.0 / rho**6
    rho12    = rho6**2
    cutterm  = 4.0*(rho12 - rho6)                 ! Needed for cutoff forces and torques
    dcutterm = -24.0 * (2.0 * rho12 - rho6) / rho ! Needed for cutoff forces and torques
    pot      = pot - epsilon * cutterm

    IF ( .NOT. PRESENT(f) ) RETURN

    IF ( .NOT. PRESENT(t) ) THEN
       WRITE ( unit=error_unit, fmt='(a)' ) 'Both f and t expected'
       STOP 'Error in test_pot_gb'
    END IF
    IF ( ANY ( SHAPE(f) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'f shape error', SHAPE(f), 3, n
       STOP 'Error in test_pot_gb'
    END IF
    IF ( ANY ( SHAPE(t) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 't shape error', SHAPE(t), 3, n
       STOP 'Error in test_pot_gb'
    END IF

    ! Derivatives of sigma
    prefac    = 0.5*chi*sigma**3
    dsig_dci  = prefac*(cpchi+cmchi)
    dsig_dcj  = prefac*(cpchi-cmchi)
    prefac    = prefac*(0.5*chi)
    dsig_dcij = -prefac*(cpchi**2-cmchi**2)

    ! Derivatives of epsilon
    prefac    = -mu*xhi*(eps1**nu)*eps2**(mu-1)
    deps_dci  = prefac*(cpxhi+cmxhi)
    deps_dcj  = prefac*(cpxhi-cmxhi)
    prefac    = prefac*(0.5*xhi)
    deps_dcij = -prefac*(cpxhi**2-cmxhi**2)                           ! From derivative of eps2
    deps_dcij = deps_dcij + nu*(chi**2)*(eps1**(nu+2))*(eps2**mu)*cij ! From derivative of eps1

    ! Derivatives of potential
    dpot_drij  = epsilon * drhoterm
    dpot_dci   = rhoterm * deps_dci  - epsilon * drhoterm * dsig_dci
    dpot_dcj   = rhoterm * deps_dcj  - epsilon * drhoterm * dsig_dcj
    dpot_dcij  = rhoterm * deps_dcij - epsilon * drhoterm * dsig_dcij

    ! Standard formula for forces and torque gradients
    fij = -dpot_drij*sij - dpot_dci*(e(:,i)-ci*sij)/rij_mag - dpot_dcj*(e(:,j)-cj*sij)/rij_mag
    gi  = dpot_dci*sij  + dpot_dcij*e(:,j)
    gj  = dpot_dcj*sij  + dpot_dcij*e(:,i)

    ! Derivatives of potential at cutoff
    dpot_drij  = epsilon * dcutterm
    dpot_dci   = cutterm * deps_dci  - epsilon * dcutterm * dsig_dci
    dpot_dcj   = cutterm * deps_dcj  - epsilon * dcutterm * dsig_dcj
    dpot_dcij  = cutterm * deps_dcij - epsilon * dcutterm * dsig_dcij

    ! Standard formula for forces and torque gradients (without dpot_drij term)
    fij(:) = fij(:) + dpot_dci*(e(:,i)-ci*sij)/rij_mag + dpot_dcj*(e(:,j)-cj*sij)/rij_mag
    gi(:)  = gi(:) - ( dpot_dci*sij  + dpot_dcij*e(:,j) ) 
    gj(:)  = gj(:) - ( dpot_dcj*sij  + dpot_dcij*e(:,i) ) 

    ! Forces
    f(:,i) = fij
    f(:,j) = -fij

    ! Torques
    t(:,i) = -cross_product ( e(:,i), gi ) ! Torque on i due to j
    t(:,j) = -cross_product ( e(:,j), gj ) ! Torque on j due to i

  END SUBROUTINE force

END MODULE test_pot_module

