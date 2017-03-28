! ewald_module.f90
! r-space and k-space parts of Ewald sum for ions
MODULE ewald_module

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <dominic.tildesley@epfl.ch> ("the authors"),                          !
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

  ! References:
  ! Woodcock and Singer, Trans. Faraday Soc. 67, 12, 1971.
  ! de Leeuw et al., Proc. Roy. Soc. A 373, 27, 1980.
  ! Heyes, J. Chem. Phys. 74, 1924, 1981.
  ! see also Fincham, mdions, CCP5 program library.

  ! The self term is subtracted from the k-space contribution
  ! The surface term for simulations in vacuum is not included
  ! A cubic box and unit box length are assumed throughout
  ! No special lists are used

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: pot_r_ewald, pot_k_ewald

  ! Consider k-vectors within a sphere bounded by the cube (-nk,+nk) in each component
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: kfac     ! k-dept quantities (k_sq_max)
  INTEGER,                         SAVE :: k_sq_max ! =nk**2

CONTAINS

  FUNCTION pot_r_ewald ( n, r, q, kappa ) RESULT ( pot )
    IMPLICIT NONE
    REAL                              :: pot   ! Returns r-space part of potential energy
    INTEGER,              INTENT(in)  :: n     ! Number of atoms
    REAL, DIMENSION(3,n), INTENT(in)  :: r     ! Positions
    REAL, DIMENSION(n),   INTENT(in)  :: q     ! Charges
    REAL,                 INTENT(in)  :: kappa ! Ewald width parameter

    INTEGER            :: i, j
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_mag, vij

    pot = 0.0

    DO i = 1, n - 1 ! Begin outer loop
       DO j = i + 1, n ! Begin inner loop

          rij(:)  = r(:,i) - r(:,j)           ! Separation vector (box=1 units)
          rij(:)  = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries
          rij_mag = SQRT ( SUM ( rij**2 ) )   ! Magnitude of separation

          vij = q(i) * q(j) * ERFC ( kappa * rij_mag ) / rij_mag ! Screened Coulomb term
          pot = pot + vij

       END DO ! End inner loop
    END DO ! End outer loop

  END FUNCTION pot_r_ewald

  FUNCTION pot_k_ewald ( nk, n, r, q, kappa ) RESULT ( pot )
    IMPLICIT NONE
    REAL                              :: pot   ! Returns k-space part of potential energy
    INTEGER,              INTENT(in)  :: nk    ! Determines number of wave-vectors
    INTEGER,              INTENT(in)  :: n     ! Number of atoms
    REAL, DIMENSION(3,n), INTENT(in)  :: r     ! Positions
    REAL, DIMENSION(n),   INTENT(in)  :: q     ! Charges
    REAL,                 INTENT(in)  :: kappa ! Ewald width parameter

    INTEGER                      :: k, kx, ky, kz, k_sq
    REAL                         :: b, factor, kr_sq
    COMPLEX, DIMENSION(n,  0:nk) :: eikx
    COMPLEX, DIMENSION(n,-nk:nk) :: eiky, eikz
    COMPLEX                      :: term

    LOGICAL, SAVE :: first_call = .TRUE.

    REAL, PARAMETER :: pi = 4.0*ATAN(1.0), twopi = 2.0*pi, twopi_sq = twopi**2

    IF ( first_call ) THEN ! Precalculation of expressions at the start

       b = 1.0 / 4.0 / kappa / kappa
       k_sq_max = nk**2            ! Store this in module data
       ALLOCATE ( kfac(k_sq_max) ) ! Allocate module data array

       ! Lazy triple loop, which over-writes the same values of ksq several times
       ! and leaves some locations empty (those which are not the sum of three squared integers)
       ! We are only interested in the value of k**2, so skip all negative components
       DO kx = 0, nk
          DO ky = 0, nk
             DO kz = 0, nk

                k_sq = kx**2 + ky**2 + kz**2

                IF ( ( k_sq <= k_sq_max ) .AND. ( k_sq /= 0 ) ) THEN ! Test to ensure within range
                   kr_sq      = twopi_sq * REAL ( k_sq )           ! k**2 in real units
                   kfac(k_sq) = twopi * EXP ( -b * kr_sq ) / kr_sq ! Stored expression for later use
                END IF ! End test to ensure within range

             END DO
          END DO
       END DO
       ! End lazy triple loop

       first_call = .FALSE.

    END IF ! End of precalculation of expressions at the start

    ! Double-check value on later calls
    IF ( k_sq_max /= nk**2 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'nk error ', k_sq_max, nk**2
       STOP 'Error in pot_k_ewald'
    END IF

    ! Construct exp(ik.r) for all ions and k-vectors

    ! Calculate kx, ky, kz = 0, 1 explicitly
    eikx(:,0) = (1.0, 0.0)
    eiky(:,0) = (1.0, 0.0)
    eikz(:,0) = (1.0, 0.0)

    eikx(:,1) = CMPLX ( COS ( twopi * r(1,:) ), SIN ( twopi * r(1,:) ) )
    eiky(:,1) = CMPLX ( COS ( twopi * r(2,:) ), SIN ( twopi * r(2,:) ) )
    eikz(:,1) = CMPLX ( COS ( twopi * r(3,:) ), SIN ( twopi * r(3,:) ) )

    ! Calculate remaining positive kx, ky and kz by recurrence
    DO k = 2, nk
       eikx(:,k) = eikx(:,k-1) * eikx(:,1)
       eiky(:,k) = eiky(:,k-1) * eiky(:,1)
       eikz(:,k) = eikz(:,k-1) * eikz(:,1)
    END DO

    ! Negative k values are complex conjugates of positive ones
    ! We do not need negative values of kx
    eiky(:,-nk:-1) = CONJG ( eiky(:,nk:1:-1) )
    eikz(:,-nk:-1) = CONJG ( eikz(:,nk:1:-1) )

    pot = 0.0

    DO kx = 0, nk ! Begin outer loop over non-negative kx

       IF ( kx == 0 ) THEN
          factor = 1.0
       ELSE
          factor = 2.0 ! Accounts for skipping negative kx
       END IF

       DO ky = -nk, nk ! Begin inner loop over all ky

          DO kz = -nk, nk ! Begin inner loop over all kz

             k_sq = kx**2 + ky**2 + kz**2

             IF ( ( k_sq <= k_sq_max ) .AND. ( k_sq /= 0 ) ) THEN ! Test to ensure within range

                term = SUM ( q(:) * eikx(:,kx) * eiky(:,ky) * eikz(:,kz) ) ! Sum over all ions
                pot  = pot + factor * kfac(k_sq) * REAL ( CONJG ( term ) * term )

             END IF ! End test to ensure within range

          END DO ! End inner loop over all kz

       END DO ! End inner loop over all ky

    END DO ! End outer loop over non-negative kx

    ! Subtract self part of k-space sum
    pot = pot - kappa * SUM ( q(:)**2 ) / SQRT(pi)

  END FUNCTION pot_k_ewald

END MODULE ewald_module
