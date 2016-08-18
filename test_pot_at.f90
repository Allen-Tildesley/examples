! test_pot_at.f90
! Axilrod-Teller triple-dipole potential
MODULE test_pot_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, force

  INTEGER, PARAMETER :: n = 3 ! three-body potential

CONTAINS

  SUBROUTINE force  ( r, pot, f )
    IMPLICIT NONE

    REAL, DIMENSION(:,:),           INTENT(in)  :: r
    REAL,                           INTENT(out) :: pot
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out) :: f

    REAL, DIMENSION(3) :: rij, rjk, rki
    REAL               :: rij_sq, rjk_sq, rki_sq, rij_mag, rjk_mag, rki_mag
    REAL               :: rij2, rjk2, rki2, ci, cj, ck, prefac, fac
    INTEGER, PARAMETER :: i = 1, j = 2, k = 3 ! notation to match appendix

    ! Routine to demonstrate the calculation of forces from the
    ! Axilrod-Teller triple-dipole potential
    ! Written for ease of comparison with the text, rather than efficiency!

    ! check dimensions to be sure
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'r shape error', SHAPE(r), 3, n
       STOP 'Error in test_pot_at'
    END IF

    rij = r(:,i) - r(:,j)
    rjk = r(:,j) - r(:,k)
    rki = r(:,k) - r(:,i)
    rij_sq = SUM(rij**2)
    rjk_sq = SUM(rjk**2)
    rki_sq = SUM(rki**2)
    rij2 = 1.0/rij_sq
    rjk2 = 1.0/rjk_sq
    rki2 = 1.0/rki_sq
    rij_mag = SQRT ( rij_sq ) ! magnitude of separation vector
    rjk_mag = SQRT ( rjk_sq ) ! magnitude of separation vector
    rki_mag = SQRT ( rki_sq ) ! magnitude of separation vector
    ci = DOT_PRODUCT ( rki, rij )
    cj = DOT_PRODUCT ( rij, rjk )
    ck = DOT_PRODUCT ( rjk, rki )
    prefac = 1/(rij_mag*rjk_mag*rki_mag)**5

    pot = prefac * ( rij_sq*rjk_sq*rki_sq - 3.0*ci*cj*ck ) ! The triple-dipole potential with strength=nu=1

    IF ( .NOT. PRESENT(f) ) RETURN

    IF ( ANY ( SHAPE(f) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'f shape error', SHAPE(f), 3, n
       STOP 'Error in test_pot_at'
    END IF

    fac = 5.0*(rij_sq*rjk_sq*rki_sq-3.0*ci*cj*ck)
    f(:,i) = prefac * ( fac*(rij2*rij-rki2*rki) &
         & + 3.0*ci*(ck-cj)*rjk + 3.0*cj*ck*(rki-rij) &
         & + 2.0*(rij_sq*rjk_sq*rki-rjk_sq*rki_sq*rij) )
    f(:,j) = prefac * ( fac*(rjk2*rjk-rij2*rij) &
         & + 3.0*cj*(ci-ck)*rki + 3.0*ck*ci*(rij-rjk) &
         & + 2.0*(rjk_sq*rki_sq*rij-rki_sq*rij_sq*rjk) )
    f(:,k) = prefac * ( fac*(rki2*rki-rjk2*rjk) &
         & + 3.0*ck*(cj-ci)*rij + 3.0*ci*cj*(rjk-rki) &
         & + 2.0*(rki_sq*rij_sq*rjk-rij_sq*rjk_sq*rki) )

  END SUBROUTINE force

END MODULE test_pot_module
  
