! test_pot_dq.f90
! dipole-quadrupole potential
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
  INTEGER, PARAMETER, PUBLIC :: n = 2 ! pair potential

CONTAINS

  SUBROUTINE force  ( r, e, pot, f, t )
    USE maths_module, ONLY : cross_product
    IMPLICIT NONE
    REAL, DIMENSION(:,:),           INTENT(in)  :: r, e
    REAL,                           INTENT(out) :: pot
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out) :: f, t

    REAL, DIMENSION(3) :: rij, sij, fij, gi, gj
    REAL               :: rij_mag, ci, cj, cij
    REAL               :: vij_dq, vij_qd, dvdrij, dvdci, dvdcj, dvdcij ! Potential derivatives
    REAL, PARAMETER    :: tol = 1.e-6
    INTEGER, PARAMETER :: i = 1, j = 2 ! notation to match appendix

    ! Routine to demonstrate the calculation of forces and torques from the
    ! dipole-quadrupole potential
    ! Written for ease of comparison with the text, rather than efficiency!

    ! Check dimensions to be sure
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'r shape error', SHAPE(r), 3, n
       STOP 'Error in test_pot_dq'
    END IF
    IF ( ANY ( SHAPE(e) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'e shape error', SHAPE(e), 3, n
       STOP 'Error in test_pot_dq'
    END IF

    ! Check unit vectors
    IF ( ABS(SUM(e(:,i)**2)-1.0) > tol .OR. ABS(SUM(e(:,j)**2)-1.0) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a)'      ) 'Warning, non-unit vectors'
       WRITE ( unit=error_unit, fmt='(4f10.5)' ) e(:,i), SUM(e(:,i)**2)
       WRITE ( unit=error_unit, fmt='(4f10.5)' ) e(:,j), SUM(e(:,j)**2)
    END IF

    rij = r(:,i) - r(:,j)
    rij_mag = SQRT ( SUM(rij**2) ) ! Magnitude of separation vector
    sij = rij / rij_mag            ! Unit vector
    ci  = DOT_PRODUCT ( e(:,i), sij )
    cj  = DOT_PRODUCT ( e(:,j), sij )
    cij = DOT_PRODUCT ( e(:,i), e(:,j) )

    ! We calculate the dipole-quadrupole and quadrupole-dipole contributions separately
    ! for clarity. Either part could be commented out so as to be tested separately
    ! together with the corresponding sections for force and torque below.
    pot = 0.0
    fij = 0.0
    gi  = 0.0
    gj  = 0.0

    ! The dipole-quadrupole potential with mu_i = 1, Q_j = 1
    vij_dq = 1.5 * (ci*(1.0-5.0*cj**2)+2*cj*cij)/rij_mag**4
    pot = pot + vij_dq

    ! The quadrupole-dipole potential with Q_i = 1, mu_j = 1
    vij_qd = -1.5 * (cj*(1.0-5.0*ci**2)+2*ci*cij)/rij_mag**4
    pot = pot + vij_qd

    IF ( .NOT. PRESENT(f) ) RETURN

    IF ( .NOT. PRESENT(t) ) THEN
       WRITE ( unit=error_unit, fmt='(a)' ) 'Both f and t expected'
       STOP 'Error in test_pot_dq'
    END IF
    IF ( ANY ( SHAPE(f) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'f shape error', SHAPE(f), 3, n
       STOP 'Error in test_pot_dq'
    END IF
    IF ( ANY ( SHAPE(t) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 't shape error', SHAPE(t), 3, n
       STOP 'Error in test_pot_dq'
    END IF

    ! Forces and torques for dipole-quadrupole potential with mu_i = 1, Q_j = 1
    dvdrij = -4.0*vij_dq/rij_mag
    dvdci  =  1.5 * (1-5.0*cj**2) / rij_mag**4
    dvdcj  =  3.0 * (cij-5.0*ci*cj) / rij_mag**4
    dvdcij =  3.0 * cj / rij_mag**4

    ! Forces
    fij = fij - dvdrij*sij - dvdci*(e(:,i)-ci*sij)/rij_mag - dvdcj*(e(:,j)-cj*sij)/rij_mag

    ! Torque gradients
    gi = gi + dvdci*sij + dvdcij*e(:,j)
    gj = gj + dvdcj*sij + dvdcij*e(:,i)

    ! Forces and torques for quadrupole-dipole potential with Q_i = 1, mu_j = 1
    dvdrij = -4.0*vij_qd/rij_mag
    dvdci  = -3.0 * (cij-5.0*ci*cj) / rij_mag**4
    dvdcj  = -1.5 * (1-5.0*ci**2) / rij_mag**4
    dvdcij = -3.0 * ci / rij_mag**4

    ! Forces
    fij = fij - dvdrij*sij - dvdci*(e(:,i)-ci*sij)/rij_mag - dvdcj*(e(:,j)-cj*sij)/rij_mag

    ! Torque gradients
    gi = gi + dvdci*sij + dvdcij*e(:,j)
    gj = gj + dvdcj*sij + dvdcij*e(:,i)

    ! Final forces and torques
    f(:,i) = fij
    f(:,j) = -fij
    t(:,i) = -cross_product(e(:,i),gi)
    t(:,j) = -cross_product(e(:,j),gj)

  END SUBROUTINE force

END MODULE test_pot_module

