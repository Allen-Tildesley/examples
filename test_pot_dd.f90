! test_pot_dd.f90
! dipole-dipole potential
MODULE test_pot_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, force

  INTEGER, PARAMETER :: n = 2 ! Pair potential

CONTAINS

  SUBROUTINE force  ( r, e, pot, f, t )
    USE utility_module, ONLY : cross_product
    IMPLICIT NONE

    REAL, DIMENSION(:,:),           INTENT(in)  :: r, e
    REAL,                           INTENT(out) :: pot
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out) :: f, t

    REAL, DIMENSION(3) :: rij, sij, fij, gi, gj
    REAL               :: rij_mag, ci, cj, cij
    REAL, PARAMETER    :: tol = 1.e-6
    INTEGER, PARAMETER :: i = 1, j = 2 ! notation to match appendix

    ! Routine to demonstrate the calculation of forces and torques from the
    ! dipole-dipole potential
    ! Written for ease of comparison with the text, rather than efficiency!

    ! check dimensions to be sure
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'r shape error', SHAPE(r), 3, n
       STOP 'Error in test_pot_dd'
    END IF
    IF ( ANY ( SHAPE(e) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'e shape error', SHAPE(e), 3, n
       STOP 'Error in test_pot_dd'
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

    ! The dipole-dipole potential with mu_i = mu_j = 1
    pot = (cij-3.0*ci*cj)/rij_mag**3

    IF ( .NOT. PRESENT(f) ) RETURN

    IF ( .NOT. PRESENT(t) ) THEN
       WRITE ( unit=error_unit, fmt='(a)' ) 'Both f and t expected'
       STOP 'Error in test_pot_dd'
    END IF
    IF ( ANY ( SHAPE(f) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'f shape error', SHAPE(f), 3, n
       STOP 'Error in test_pot_dd'
    END IF
    IF ( ANY ( SHAPE(t) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 't shape error', SHAPE(t), 3, n
       STOP 'Error in test_pot_dd'
    END IF

    ! Forces
    fij = 3.0 * ( (cij-5.0*ci*cj)*sij + cj*e(:,i) + ci*e(:,j) ) / rij_mag**4

    ! Torque gradients
    gi = ( e(:,j) - 3.0*cj*sij ) / rij_mag**3
    gj = ( e(:,i) - 3.0*ci*sij ) / rij_mag**3

    ! Final forces and torques
    f(:,i) = fij
    f(:,j) = -fij
    t(:,i) = -cross_product(e(:,i),gi)
    t(:,j) = -cross_product(e(:,j),gj)

  END SUBROUTINE force

END MODULE test_pot_module
  
