! initialize_module.f90
! Routines to initialize configurations and velocities
MODULE initialize_module

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

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: allocate_arrays, deallocate_arrays
  PUBLIC :: initialize_lattice, initialize_random, initialize_velocities

  ! Public data
  INTEGER,                           PUBLIC :: n      ! Number of atoms
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r      ! Positions (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v      ! Velocities (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: e      ! Orientations (3,n) or (0:3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: w      ! Angular velocities (3,n)

CONTAINS

  SUBROUTINE allocate_arrays ( quaternions )
    IMPLICIT NONE
    LOGICAL, INTENT(in) :: quaternions

    ALLOCATE ( r(3,n), v(3,n), w(3,n) )

    IF ( quaternions ) THEN
       ALLOCATE ( e(0:3,n) )
    ELSE
       ALLOCATE ( e(3,n) )
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, v, w, e )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE initialize_lattice ( box, length, soft )
    IMPLICIT NONE
    REAL,    INTENT(in) :: box    ! Simulation box length
    REAL,    INTENT(in) :: length ! Molecule length
    LOGICAL, INTENT(in) :: soft   ! Flag for soft interactions (no overlap check)

    ! Sets up the fcc lattice: four molecules per unit cell
    ! For atoms, for which length=0.0, the e-coordinates will be ignored
    ! For linear molecules, the orientations comply with the alpha-fcc pattern
    ! For nonlinear molecules, the 0-element is set to zero

    REAL, DIMENSION(3,4), PARAMETER :: r_fcc = RESHAPE ( [ &
         & 0.25, 0.25, 0.25, &
         & 0.25, 0.75, 0.75, &
         & 0.75, 0.75, 0.25, &
         & 0.75, 0.25, 0.75 ], [3,4] ) ! Positions in unit cell

    REAL, DIMENSION(3,4), PARAMETER :: e_fcc = RESHAPE ( SQRT(1.0/3.0) * [ &
         &  1.0,  1.0,  1.0,  &
         &  1.0, -1.0, -1.0,  &
         & -1.0,  1.0, -1.0,  &
         & -1.0, -1.0,  1.0  ], [3,4] ) ! Orientations in unit cell

    REAL    :: cell, box2
    INTEGER :: nc, ix, iy, iz, a, i

    WRITE ( unit=output_unit, fmt='(a)' ) 'Close-packed fcc lattice positions'

    nc = NINT ( REAL(n/4) ** (1.0/3.0) )
    IF ( n /= 4 * nc ** 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'n, nc mismatch ', n, 4 * nc ** 3
       STOP 'Error in initialize_lattice'
    END IF

    cell = box / REAL(nc) ! Unit cell
    box2 = box / 2.0      ! Half box length

    IF ( LBOUND(e,dim=1) == 0 ) e(0,:) = 0.0 ! For quaternions

    i = 0

    ! Begin triple loop over unit cell indices
    DO iz = 0, nc-1
       DO iy = 0, nc-1
          DO ix = 0, nc-1

             DO a = 1, 4 ! Begin loop over atoms in unit cell
                i = i + 1
                r(1:3,i) = r_fcc(:,a) + REAL ( [ ix, iy, iz ] ) ! In range 0 .. real(nc)
                r(1:3,i) = r(1:3,i) * cell                      ! In range 0 .. box
                r(1:3,i) = r(1:3,i) - box2                      ! In range -box/2 .. box/2
                e(1:3,i) = e_fcc(:,a)
                IF ( .NOT. soft ) THEN
                   IF ( overlap ( i, 1, i-1, box, length ) ) THEN
                      WRITE ( unit=error_unit, fmt='(a)' ) 'Density too high'
                      STOP 'Error in initialize_lattice'
                   END IF
                END IF
             END DO ! End loop over atoms in unit cell

          END DO
       END DO
    END DO
    ! End triple loop over unit cell indices

  END SUBROUTINE initialize_lattice

  SUBROUTINE initialize_random ( box, length, soft )
    USE maths_module, ONLY : random_vector, random_quaternion
    IMPLICIT NONE
    REAL, INTENT(in)    :: box    ! Simulation box length
    REAL,    INTENT(in) :: length ! Molecule length
    LOGICAL, INTENT(in) :: soft   ! Flag for soft interactions (no overlap check)

    ! Places atoms at random positions
    ! Unlikely to be useful, unless the interaction potential is soft
    ! or the density rather low

    INTEGER, PARAMETER :: iter_max = 1000 ! Max random placement iterations
    INTEGER :: i, iter

    WRITE ( unit=output_unit, fmt='(a)' ) 'Random positions'

    DO i = 1, n

       iter = 0
       DO ! Loop until non-overlapping position found

          CALL RANDOM_NUMBER ( r(:,i) )   ! All in range (0,1)
          r(:,i) = ( r(:,i) - 0.5 ) * box ! In range -box/2 .. box/2
          IF ( LBOUND(e,dim=1) == 0 ) THEN
             e(:,i) = random_quaternion ( )
          ELSE
             e(:,i) = random_vector ( )
          END IF

          IF ( soft ) EXIT ! No overlap test
          IF ( .NOT. overlap ( i, 1, i-1, box, length ) ) EXIT

          iter = iter + 1
          IF ( iter > iter_max ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations ', iter, iter_max
             STOP 'Error in initialize_random'
          END IF
       END DO ! End loop until non-overlapping position found

    END DO

  END SUBROUTINE initialize_random

  SUBROUTINE initialize_velocities ( temperature, inertia )
    USE maths_module, ONLY : random_perpendicular_vector, random_normals
    IMPLICIT NONE
    REAL, INTENT(in) :: temperature ! Reduced temperature
    REAL, INTENT(in) :: inertia     ! Reduced moment of inertia

    ! Chooses linear velocities from Maxwell-Boltzmann (Gaussian) distribution
    ! We set the total momentum to zero
    ! We assume unit molecular mass

    ! For linear molecules we choose the direction of the angular velocity
    ! randomly but perpendicular to the molecular axis.
    ! The square of the magnitude of the angular velocity
    ! is chosen from an exponential distribution
    ! For nonlinear molecules we choose all three components of angular velocity
    ! from a Gaussian distribution, assuming equal moments of inertia
    ! There is no attempt to set the total angular momentum to zero

    REAL               :: factor
    REAL, DIMENSION(3) :: v_cm
    REAL    ::  w_sq, w_sq_mean, w_std_dev, zeta
    INTEGER ::  i

    WRITE ( unit=output_unit, fmt='(a,t40,2f15.6)' ) 'Velocities at temperature, inertia', temperature, inertia

    ! Linear velocities

    CALL random_normals ( 0.0, 1.0, v ) ! Unit normal random numbers

    v_cm(:) = SUM ( v(:,:), dim=2 ) / REAL ( n )            ! Compute centre of mass velocity
    v(:,:)  = v(:,:) - SPREAD ( v_cm(:), dim=2, ncopies=n ) ! Set net momentum to zero
    factor  = SQRT ( REAL(3*n-3)*temperature / SUM (v**2) ) ! Sqrt of ratio of kinetic energies
    v       = factor * v

    ! Angular velocities

    IF ( LBOUND(e,dim=1) == 0 ) THEN ! Nonlinear molecule, treat as spherical top

       w_std_dev = SQRT(temperature/inertia)
       CALL random_normals ( 0.0, w_std_dev, w )

    ELSE ! Linear molecule

       w_sq_mean = 2.0 * temperature / inertia

       DO i = 1, n ! Begin loop over molecules
          w(:,i) = random_perpendicular_vector ( e(:,i) ) ! Set direction of the angular velocity
          CALL RANDOM_NUMBER ( zeta )
          w_sq   = - w_sq_mean * LOG ( zeta ) ! Squared magnitude of angular velocity
          w(:,i) = w(:,i) * SQRT ( w_sq )
       END DO ! End loop over molecules

    END IF

  END SUBROUTINE initialize_velocities

  FUNCTION overlap ( i, j1, j2, box, ell )
    IMPLICIT NONE
    LOGICAL             :: overlap ! Returns a flag indicating overlap with any j<i
    INTEGER, INTENT(in) :: i       ! Index of atom/molecule to be checked
    INTEGER, INTENT(in) :: j1      ! First j index
    INTEGER, INTENT(in) :: j2      ! Last j index
    REAL,    INTENT(in) :: box     ! Simulation box length for PBC
    REAL,    INTENT(in) :: ell     ! Molecule length

    ! This routine checks for overlaps of atoms (ell=0) or spherocylinders (ell>0)

    INTEGER            :: j
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq, rei, rej, eij, sij_sq
    REAL               :: sin_sq, ci, cj, ai, aj, di, dj, ell2
    REAL, PARAMETER    :: tol = 1.0e-6

    overlap = .FALSE.
    IF ( j2 < j1 ) RETURN

    IF ( ABS ( ell ) < tol ) THEN ! Handle spherical case separately

       DO j = j1, j2
          rij    = r(:,i) - r(:,j)
          rij    = rij - ANINT ( rij / box ) * box
          rij_sq = SUM ( rij**2 )

          IF ( rij_sq < 1.0 ) THEN
             overlap = .TRUE.
             RETURN
          END IF

       END DO

    ELSE ! Handle non-spherical case

       ell2   = ell / 2.0    ! Half the line segment length

       DO j = j1, j2

          rij    = r(:,i) - r(:,j)
          rij    = rij - ANINT ( rij / box ) * box
          rij_sq = SUM ( rij**2 )
          rei    = DOT_PRODUCT ( rij, e(:,i) )
          rej    = DOT_PRODUCT ( rij, e(:,j) )
          eij    = DOT_PRODUCT ( e(:,i), e(:,j ) )

          sin_sq = 1.0 - eij**2 ! Squared sine of angle between line segments

          IF ( sin_sq < tol ) THEN ! Guard against nearly-parallel lines
             ci = -rei
             cj =  rej
          ELSE
             ci = ( - rei + eij * rej ) / sin_sq
             cj = (   rej - eij * rei ) / sin_sq
          ENDIF

          ai = ABS ( ci )
          aj = ABS ( cj )
          IF ( ai > ell2 ) ci = SIGN ( ell2, ci )
          IF ( aj > ell2 ) cj = SIGN ( ell2, cj )

          IF ( ai > aj ) THEN
             cj =  rej + ci * eij
          ELSE
             ci = -rei + cj * eij
          ENDIF

          ai = ABS ( ci )
          aj = ABS ( cj )
          IF ( ai > ell2 ) ci = SIGN ( ell2, ci )
          IF ( aj > ell2 ) cj = SIGN ( ell2, cj )

          di =  2.0 * rei + ci - cj * eij
          dj = -2.0 * rej + cj - ci * eij

          sij_sq = rij_sq + ci * di + cj * dj ! Squared distance between line segments

          IF ( sij_sq < 1.0 ) THEN
             overlap = .TRUE.
             RETURN
          END IF
       END DO
    END IF

  END FUNCTION overlap

END MODULE initialize_module
