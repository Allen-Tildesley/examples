! test_pot_qq.f90
! Pair potential, quadrupole-quadrupole, quadrupole moment Q=1
MODULE test_pot_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

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

    REAL, DIMENSION(3) :: rij, sij, fij, gi, gj
    REAL               :: rij_mag, ci, cj, cij
    REAL               :: vij, dvdrij, dvdci, dvdcj, dvdcij ! potential derivatives
    REAL, PARAMETER    :: tol = 1.e-6
    INTEGER, PARAMETER :: i = 1, j = 2, n = 2 ! notation to match appendix

    ! Routine to demonstrate the calculation of forces and torques from the
    ! quadrupole-quadrupole potential
    ! Written for ease of comparison with the text, rather than efficiency!

    ! check dimensions to be sure
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'r shape error', SHAPE(r), 3, n
       STOP 'Error in test_pot_qq'
    END IF
    IF ( ANY ( SHAPE(e) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'e shape error', SHAPE(e), 3, n
       STOP 'Error in test_pot_qq'
    END IF

    ! Check unit vectors
    IF ( ABS(SUM(e(:,i)**2)-1.0) > tol .OR. ABS(SUM(e(:,j)**2)-1.0) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a)'      ) 'Warning, non-unit vectors'
       WRITE ( unit=error_unit, fmt='(4f10.5)' ) e(:,i), SUM(e(:,i)**2)
       WRITE ( unit=error_unit, fmt='(4f10.5)' ) e(:,j), SUM(e(:,j)**2)
    END IF

    rij = r(:,i) - r(:,j)
    rij_mag = SQRT ( SUM(rij**2) ) ! magnitude of separation vector
    sij = rij / rij_mag ! unit vector
    ci  = DOT_PRODUCT ( e(:,i), sij )
    cj  = DOT_PRODUCT ( e(:,j), sij )
    cij = DOT_PRODUCT ( e(:,i), e(:,j) )
    
    ! The quadrupole-quadrupole potential with Q_i = 1, Q_j = 1
    vij = 0.75 * (1.0 - 5.0*ci**2 - 5.0*cj**2 + 2.0*cij**2 &
         & + 35.0*(ci*cj)**2 - 20.0*ci*cj*cij) / rij_mag**5
    pot = vij

    IF ( .NOT. PRESENT(f) ) RETURN

    IF ( .NOT. PRESENT(t) ) THEN
       WRITE ( unit=error_unit, fmt='(a)' ) 'Both f and t expected'
       STOP 'Error in test_pot_qq'
    END IF
    IF ( ANY ( SHAPE(f) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'f shape error', SHAPE(f), 3, n
       STOP 'Error in test_pot_qq'
    END IF
    IF ( ANY ( SHAPE(t) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 't shape error', SHAPE(t), 3, n
       STOP 'Error in test_pot_qq'
    END IF

    ! Forces and torques for dipole-quadrupole potential with mu_i = 1, Q_j = 1
    dvdrij = -5.0 * vij / rij_mag
    dvdci  =  7.5 * (ci*(7.0*cj**2-1.0)-2.0*cj*cij) / rij_mag**5
    dvdcj  =  7.5 * (cj*(7.0*ci**2-1.0)-2.0*ci*cij) / rij_mag**5
    dvdcij = -3.0 * (5.0*ci*cj-cij) / rij_mag**5

    ! Forces
    fij = - dvdrij*sij - dvdci*(e(:,i)-ci*sij)/rij_mag - dvdcj*(e(:,j)-cj*sij)/rij_mag

    ! Torque gradients
    gi = dvdci*sij + dvdcij*e(:,j)
    gj = dvdcj*sij + dvdcij*e(:,i)

    ! Final forces and torques
    f(:,i) = fij
    f(:,j) = -fij
    t(:,i) = -cross_product(e(:,i),gi)
    t(:,j) = -cross_product(e(:,j),gj)

  END SUBROUTINE force

END MODULE test_pot_module
  
