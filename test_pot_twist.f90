! test_pot_twist.f90
! twist angle, cos(phi) potential
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
  INTEGER, PARAMETER, PUBLIC :: n = 4 ! Four-body potential

CONTAINS

  SUBROUTINE force  ( r, pot, f )
    IMPLICIT NONE
    REAL, DIMENSION(:,:),           INTENT(in)  :: r
    REAL,                           INTENT(out) :: pot
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out) :: f

    ! We choose to make the polymer the minimum length needed for testing
    INTEGER                  :: a, b
    REAL, DIMENSION(3,2:n)   :: d      ! the d vectors
    REAL, DIMENSION(2:n,2:n) :: cc, dd ! the C and D coefficients
    REAL                     :: prefac, fac, fac1, fac2

    ! Routine to demonstrate the calculation of forces from the
    ! polymer angle-twisting potential
    ! Written for ease of comparison with the text, rather than efficiency!

    ! Check dimensions to be sure
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'r shape error', SHAPE(r), 3, n
       STOP 'Error in test_pot_twist'
    END IF

    ! Set up d vectors
    DO a = 2, n
       d(:,a) = r(:,a) - r(:,a-1)
    END DO

    ! Store C coefficients in a matrix
    ! In the general case we would not need to calculate every pair
    ! and also we would make use of the symmetry cc(b,a)=cc(a,b)
    DO a = 2, n
       DO b = 2, n
          cc(a,b) = DOT_PRODUCT ( d(:,a), d(:,b) )
       END DO
    END DO

    ! Store D coefficients in a matrix
    ! In the general case we would not need to calculate every pair
    ! and also we would make use of the symmetry dd(b,a)=dd(a,b)
    DO a = 2, n
       DO b = 2, n
          dd(a,b) = cc(a,a)*cc(b,b) - cc(a,b)**2
       END DO
    END DO

    ! For this test there is just one angle
    a = n

    ! Here is the potential as a function of cos(phi)
    ! For testing we use the simplest form: v= -cos(phi)
    ! The notation matches that used in the appendix

    prefac = 1.0 /  SQRT(dd(a,a-1)*dd(a-1,a-2))
    fac = (cc(a,a-1)*cc(a-1,a-2)-cc(a,a-2)*cc(a-1,a-1))
    pot = prefac * fac ! This is -cos(phi)

    IF ( .NOT. PRESENT(f) ) RETURN

    IF ( ANY ( SHAPE(f) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'f shape error', SHAPE(f), 3, n
       STOP 'Error in test_pot_twist'
    END IF

    ! Here we include the derivative of the potential with respect to cos(phi) in the prefactor
    ! For this simple case it is -1, so the forces are simply gradients of cos(phi) as in the text

    fac1 = fac/dd(a,a-1)
    fac2 = fac/dd(a-1,a-2)
    f(:,a) = -prefac * ( cc(a-1,a-2)*d(:,a-1) - cc(a-1,a-1)*d(:,a-2) &
         & -fac1 * ( cc(a-1,a-1)*d(:,a) - cc(a,a-1)*d(:,a-1) ) )
    f(:,a-1) = -prefac * ( cc(a-1,a-2)*d(:,a) - cc(a-1,a-2)*d(:,a-1) &
         & + cc(a,a-1)*d(:,a-2) + cc(a-1,a-1)*d(:,a-2) - 2.0*cc(a,a-2)*d(:,a-1) &
         & -fac2 * ( cc(a-2,a-2)*d(:,a-1) - cc(a-1,a-2)*d(:,a-2) ) &
         & -fac1 * ( cc(a,a)*d(:,a-1) - cc(a-1,a-1)*d(:,a) - cc(a,a-1)*d(:,a) + cc(a,a-1)*d(:,a-1)) )
    f(:,a-2) = -prefac*( -cc(a-1,a-2)*d(:,a) + cc(a,a-1)*d(:,a-1)  &
         & -cc(a,a-1)*d(:,a-2) - cc(a-1,a-1)*d(:,a) + 2.0*cc(a,a-2)*d(:,a-1) &
         & -fac2 * ( cc(a-1,a-1)*d(:,a-2) - cc(a-2,a-2)*d(:,a-1) - cc(a-1,a-2)*d(:,a-1) + cc(a-1,a-2)*d(:,a-2) ) &
         & -fac1 * ( -cc(a,a)*d(:,a-1) + cc(a,a-1)*d(:,a)))
    f(:,a-3) = -prefac * ( -cc(a,a-1)*d(:,a-1) + cc(a-1,a-1)*d(:,a) &
         & -fac2 * ( -cc(a-1,a-1)*d(:,a-2) + cc(a-1,a-2)*d(:,a-1) ) )

  END SUBROUTINE force

END MODULE test_pot_module

