! ewald_module.f90
! r-space and k-space parts of Ewald sum for ions
MODULE ewald_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pot_r_ewald, pot_k_ewald
  
  ! References:
  ! Woodcock and Singer, Trans. Faraday Soc. 67, 12, 1971.
  ! de Leeuw et al., Proc. Roy. Soc. A 373, 27, 1980.
  ! Heyes, J. Chem. Phys. 74, 1924, 1981.
  ! see also Fincham, mdions, CCP5 program library.

  ! The self term is subtracted from the k-space contribution
  ! The surface term for simulations in vacuum is not included
  ! A cubic box and unit box length are assumed throughout
  ! No special lists are used

CONTAINS

  FUNCTION pot_r_ewald ( n, r, q, kappa ) RESULT ( pot ) ! r-space part of potential energy
    REAL                              :: pot   ! Result
    INTEGER,              INTENT(in)  :: n     ! Number of atoms
    REAL, DIMENSION(3,n), INTENT(in)  :: r     ! positions
    REAL, DIMENSION(n),   INTENT(in)  :: q     ! charges
    REAL,                 INTENT(in)  :: kappa ! Ewald width parameter

    INTEGER            :: i, j
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_mag, vij

    pot = 0.0

    DO i = 1, n - 1 ! Begin outer loop
       DO j = i + 1, n ! Begin inner loop

          rij(:) = r(:,i) - r(:,j)
          rij(:) = rij(:) - ANINT ( rij(:) )

          rij_mag = SQRT ( SUM ( rij**2 ) )
          vij     = q(i) * q(j) * ERFC ( kappa * rij_mag ) / rij_mag
          pot     = pot + vij

       END DO ! End inner loop
    END DO ! End outer loop

  END FUNCTION pot_r_ewald

  FUNCTION pot_k_ewald ( nk, n, r, q, kappa ) RESULT ( pot ) ! k-space part of potential energy
    REAL                              :: pot   ! Result
    INTEGER,              INTENT(in)  :: nk    ! Determines number of wave-vectors
    INTEGER,              INTENT(in)  :: n     ! Number of atoms
    REAL, DIMENSION(3,n), INTENT(in)  :: r     ! positions
    REAL, DIMENSION(n),   INTENT(in)  :: q     ! charges
    REAL,                 INTENT(in)  :: kappa ! Ewald width parameter

    ! Consider k-vectors within a sphere bounded by the cube (-nk,+nk) in each component
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: kfac     ! k-dept quantities (k_sq_max)
    INTEGER,                         SAVE :: k_sq_max ! =nk**2

    INTEGER                        :: k, kx, ky, kz, k_sq
    REAL                           :: b, factor, kr_sq
    COMPLEX, DIMENSION(1:n,  0:nk) :: eikx
    COMPLEX, DIMENSION(1:n,-nk:nk) :: eiky, eikz
    COMPLEX                        :: term

    LOGICAL, SAVE :: first_call = .TRUE.

    REAL, PARAMETER :: pi = 4.0*ATAN(1.0), twopi = 2.0*pi

    IF ( first_call ) THEN ! Precalculation of expressions at the start

       b = 1.0 / 4.0 / kappa / kappa
       k_sq_max = nk**2
       ALLOCATE ( kfac(k_sq_max) )

       ! Lazy triple loop, which over-writes the same values of ksq several times
       ! and leaves some locations empty (those which are not the sum of three squares)
       ! We are only interested in the value of k**2, so skip all negative components
       DO kx = 0, nk
          DO ky = 0, nk
             DO kz = 0, nk
                k_sq = kx**2 + ky**2 + kz**2
                IF ( ( k_sq <= k_sq_max ) .AND. ( k_sq /= 0 ) ) THEN
                   kr_sq = ( twopi**2 ) * REAL ( k_sq )            ! k**2 in real units
                   kfac(k_sq) = twopi * EXP ( -b * kr_sq ) / kr_sq ! Stored expression for later use
                END IF
             END DO
          END DO
       END DO

       first_call = .FALSE.

    END IF ! End of precalculation of expressions at the start

    ! Double-check value on later calls
    IF ( k_sq_max /= nk**2 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'nk error ', k_sq_max, nk**2
       STOP 'Error in pot_k_ewald'
    END IF

    ! construct EXP(ik.r) for all ions and k-vectors **

    ! calculate kx, ky, kz = 0, 1 explicitly
    eikx(:,0) = (1.0, 0.0)
    eiky(:,0) = (1.0, 0.0)
    eikz(:,0) = (1.0, 0.0)

    eikx(:,1) = CMPLX ( COS ( twopi * r(1,:) ), SIN ( twopi * r(1,:) ) )
    eiky(:,1) = CMPLX ( COS ( twopi * r(2,:) ), SIN ( twopi * r(2,:) ) )
    eikz(:,1) = CMPLX ( COS ( twopi * r(3,:) ), SIN ( twopi * r(3,:) ) )

    ! calculate remaining kx, ky and kz by recurrence
    DO k = 2, nk
       eikx(:,k) = eikx(:,k-1) * eikx(:,1)
       eiky(:,k) = eiky(:,k-1) * eiky(:,1)
       eikz(:,k) = eikz(:,k-1) * eikz(:,1)
    END DO

    ! Negative k values are complex conjugates of positive ones
    eiky(:,-k:-1) = CONJG ( eiky(:,k:1:-1) )
    eikz(:,-k:-1) = CONJG ( eikz(:,k:1:-1) )

    pot = 0.0

    DO kx = 0, nk ! begin outer loop over non-negative kx

       IF ( kx == 0 ) THEN
          factor = 1.0
       ELSE
          factor = 2.0 ! accounts for skipping negative kx
       END IF

       DO ky = -nk, nk ! Begin inner loop over all ky

          DO kz = -nk, nk ! Begin inner loop over all kz

             k_sq = kx**2 + ky**2 + kz**2

             IF ( ( k_sq <= k_sq_max ) .AND. ( k_sq /= 0 ) ) THEN
                term = SUM ( q(:) * eikx(:,kx) * eiky(:,ky) * eikz(:,kz) )
                pot = pot + factor * kfac(k_sq) * REAL ( CONJG ( term ) * term )
             END IF

          END DO ! End inner loop over all kz
       END DO ! End inner loop over all ky
    END DO ! End outer loop over non-negative kx

    ! subtract self part of k-space sum
    pot = pot - kappa * SUM ( q(:)**2 ) / SQRT(pi)

  END FUNCTION pot_k_ewald

END MODULE ewald_module
