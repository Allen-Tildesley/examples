! test_pot_gb.f90
! Pair potential, Gay-Berne
MODULE test_pot_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, force

  INTEGER, PARAMETER :: n = 2 ! pair potential

CONTAINS

  SUBROUTINE force  ( r, e, pot, f, t )
    USE utility_module, ONLY : cross_product
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
    REAL,    PARAMETER :: r_cut = 4.0   ! normally would use a larger value of r_cut
    REAL,    PARAMETER :: tol = 1.e-6
    INTEGER, PARAMETER :: i = 1, j = 2 ! notation to match appendix

    ! Routine to demonstrate the calculation of forces and torques from the
    ! Gay-Berne potential, including the spherical cutoff contribution.
    ! Written for ease of comparison with the text, rather than efficiency!

    ! check dimensions to be sure
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) STOP 'r shape error'
    IF ( ANY ( SHAPE(e) /= [3,n] ) ) STOP 'e shape error'

    ! Check unit vectors
    IF ( ABS(SUM(e(:,i)**2)-1.0) > tol .OR. ABS(SUM(e(:,j)**2)-1.0) > tol ) THEN
       WRITE(*,'(a)') 'Warning, non-unit vectors'
       WRITE(*,'(4f10.5)') e(:,i), SUM(e(:,i)**2)
       WRITE(*,'(4f10.5)') e(:,j), SUM(e(:,j)**2)
    END IF

    rij     = r(:,i) - r(:,j)
    rij_sq  = SUM ( rij**2 )
    rij_mag = SQRT(rij_sq)
    sij     = rij / rij_mag ! unit vector along rij

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
    eps1    = 1.0/SQRT(1.0-(chi*cij)**2) ! depends on chi, not xhi
    cpxhi   = cp/(1.0+xhi*cij)
    cmxhi   = cm/(1.0-xhi*cij)
    eps2    = 1.0-0.5*xhi*(cp*cpxhi+cm*cmxhi) ! depends on xhi
    epsilon = (eps1**nu) * (eps2**mu)

    ! Potential at rij
    rho      = rij_mag - sigma + 1.0
    rho6     = 1.0 / rho**6
    rho12    = rho6**2
    rhoterm  = 4.0*(rho12 - rho6)                 ! needed for forces and torques
    drhoterm = -24.0 * (2.0 * rho12 - rho6) / rho ! needed for forces and torques
    pot      = epsilon*rhoterm

    ! Potential at r_cut
    rho      = r_cut - sigma + 1.0
    rho6     = 1.0 / rho**6
    rho12    = rho6**2
    cutterm  = 4.0*(rho12 - rho6)                 ! needed for cutoff forces and torques
    dcutterm = -24.0 * (2.0 * rho12 - rho6) / rho ! needed for cutoff forces and torques
    pot      = pot - epsilon * cutterm

    IF ( .NOT. PRESENT(f) ) RETURN
    IF ( .NOT. PRESENT(t) ) STOP 'both f and t expected'
    IF ( ANY ( SHAPE(f) /= [3,n] ) ) STOP 'f shape error'
    IF ( ANY ( SHAPE(t) /= [3,n] ) ) STOP 't shape error'

    ! Derivatives of sigma
    prefac   = 0.5*chi*sigma**3
    dsig_dci  = prefac*(cpchi+cmchi)
    dsig_dcj  = prefac*(cpchi-cmchi)
    prefac   = prefac*(0.5*chi)
    dsig_dcij = -prefac*(cpchi**2-cmchi**2)

    ! Derivatives of epsilon
    prefac   = -mu*xhi*(eps1**nu)*eps2**(mu-1)
    deps_dci  = prefac*(cpxhi+cmxhi)
    deps_dcj  = prefac*(cpxhi-cmxhi)
    prefac   = prefac*(0.5*xhi)
    deps_dcij = -prefac*(cpxhi**2-cmxhi**2)                          ! from derivative of eps2
    deps_dcij = deps_dcij + nu*(chi**2)*(eps1**(nu+2))*(eps2**mu)*cij ! from derivative of eps1

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
  
