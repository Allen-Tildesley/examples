! test_pot_at.f90
! Axilrod-Teller triple-dipole potential
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
  INTEGER, PARAMETER, PUBLIC :: n = 3 ! Three-body potential

CONTAINS

  SUBROUTINE force  ( r, pot, f )
    IMPLICIT NONE
    REAL, DIMENSION(:,:),           INTENT(in)  :: r
    REAL,                           INTENT(out) :: pot
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out) :: f

    REAL, DIMENSION(3) :: rij, rjk, rki
    REAL               :: rij_sq, rjk_sq, rki_sq, rij_mag, rjk_mag, rki_mag
    REAL               :: rij2, rjk2, rki2, ci, cj, ck, prefac, fac
    INTEGER, PARAMETER :: i = 1, j = 2, k = 3 ! Notation to match appendix

    ! Routine to demonstrate the calculation of forces from the
    ! Axilrod-Teller triple-dipole potential
    ! Written for ease of comparison with the text, rather than efficiency!

    ! Check dimensions to be sure
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'r shape error', SHAPE(r), 3, n
       STOP 'Error in test_pot_at'
    END IF

    ! Note that we define the separation vectors in a cyclic way
    rij = r(:,i) - r(:,j)
    rjk = r(:,j) - r(:,k)
    rki = r(:,k) - r(:,i)
    rij_sq = SUM(rij**2)
    rjk_sq = SUM(rjk**2)
    rki_sq = SUM(rki**2)
    rij2 = 1.0/rij_sq
    rjk2 = 1.0/rjk_sq
    rki2 = 1.0/rki_sq
    rij_mag = SQRT ( rij_sq )
    rjk_mag = SQRT ( rjk_sq )
    rki_mag = SQRT ( rki_sq )
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

