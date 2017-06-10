! initialize_module.f90
! Routines to initialize configurations and velocities
MODULE initialize_module

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

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: allocate_arrays, deallocate_arrays
  PUBLIC :: fcc_positions, ran_positions, ran_velocities
  PUBLIC :: chain_positions, chain_velocities

  ! Public data
  INTEGER,                           PUBLIC :: n ! Number of atoms
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v ! Velocities (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: e ! Orientations (3,n) or (0:3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: w ! Angular velocities (3,n)

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

  SUBROUTINE fcc_positions ( box, length, soft )
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
       STOP 'Error in fcc_positions'
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
                      STOP 'Error in fcc_positions'
                   END IF
                END IF
             END DO ! End loop over atoms in unit cell

          END DO
       END DO
    END DO
    ! End triple loop over unit cell indices

  END SUBROUTINE fcc_positions

  SUBROUTINE ran_positions ( box, length, soft )
    USE maths_module, ONLY : random_vector, random_quaternion
    IMPLICIT NONE
    REAL, INTENT(in)    :: box    ! Simulation box length
    REAL,    INTENT(in) :: length ! Molecule length
    LOGICAL, INTENT(in) :: soft   ! Flag for soft interactions (no overlap check)

    ! Places molecules at random positions
    ! Unlikely to be useful, unless the interaction potential is soft
    ! or the density rather low
    ! For atoms, for which length=0.0, the e-coordinates will be ignored

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
             STOP 'Error in ran_positions'
          END IF
       END DO ! End loop until non-overlapping position found

    END DO

  END SUBROUTINE ran_positions

  SUBROUTINE ran_velocities ( temperature, inertia )
    USE maths_module, ONLY : random_perpendicular_vector, random_normals
    IMPLICIT NONE
    REAL, INTENT(in) :: temperature ! Reduced temperature
    REAL, INTENT(in) :: inertia     ! Reduced moment of inertia

    ! Chooses translational velocities from Maxwell-Boltzmann (Gaussian) distribution
    ! We set the total momentum to zero
    ! We assume unit molecular mass

    ! For linear molecules we choose the direction of the angular velocity
    ! randomly but perpendicular to the molecular axis.
    ! The square of the magnitude of the angular velocity
    ! is chosen from an exponential distribution
    ! For nonlinear molecules we choose all three components of angular velocity
    ! from a Gaussian distribution, assuming equal moments of inertia
    ! There is no attempt to set the total angular momentum to zero
    ! For atoms, the w array is set here, but ignored later

    REAL               :: factor
    REAL, DIMENSION(3) :: v_cm
    REAL               :: w_sq, w_sq_mean, w_std_dev, zeta
    INTEGER            :: i
    REAL, PARAMETER    :: tol = 1.0e-6

    WRITE ( unit=output_unit, fmt='(a,t40,2f15.6)' ) 'Velocities at temperature, inertia', temperature, inertia

    ! Translational velocities

    CALL random_normals ( 0.0, 1.0, v ) ! Unit normal random numbers

    v_cm(:) = SUM ( v(:,:), dim=2 ) / REAL ( n )            ! Compute centre of mass velocity
    v(:,:)  = v(:,:) - SPREAD ( v_cm(:), dim=2, ncopies=n ) ! Set net momentum to zero
    factor  = SQRT ( REAL(3*n-3)*temperature / SUM (v**2) ) ! Sqrt of ratio of kinetic energies
    v       = factor * v

    ! Angular velocities (will be ignored in the case of atoms)

    IF ( inertia < tol ) THEN ! Should never happen
       WRITE ( unit=error_unit, fmt='(a,f15.6)' ) 'Error, inertia = ', inertia
       STOP 'Error in ran_velocities'
    END IF

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

  END SUBROUTINE ran_velocities

  SUBROUTINE chain_positions ( bond, soft )
    USE maths_module, ONLY : random_vector
    IMPLICIT NONE
    REAL,    INTENT(in) :: bond ! Chain bond length
    LOGICAL, INTENT(in) :: soft ! Flag for soft interactions (no overlap check)

    ! Chooses chain positions randomly, at desired bond length, avoiding overlap

    REAL               :: diff_sq
    REAL, DIMENSION(3) :: r_cm
    INTEGER            :: i, iter

    REAL,    PARAMETER :: tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 1000

    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Chain, randomly oriented bonds = ', bond

    r(:,1) = [0.0,0.0,0.0]        ! First atom at origin
    r(:,2) = bond*random_vector() ! Second atom at random position (bond length away)

    DO i = 3, n ! Begin loop over atom indices

       iter = 0
       DO ! Loop until non-overlapping position found (N must not be too large!)
          r(:,i) = r(:,i-1) + bond*random_vector() ! Subsequent atoms randomly placed (bond length away)

          IF ( soft ) EXIT ! No overlap test
          IF ( .NOT. chain_overlap ( i, 1, i-2 ) ) EXIT ! Check all so far except bonded neighbour

          iter = iter + 1
          IF ( iter > iter_max ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations ', iter, iter_max
             STOP 'Error in chain_positions'
          END IF

       END DO ! End loop until non-overlapping position found

    END DO ! End loop over atom indices

    r_cm(:) = SUM ( r(:,:), dim=2 ) / REAL(n)               ! Compute centre of mass positions
    r(:,:)  = r(:,:) - SPREAD ( r_cm(:), dim=2, ncopies=n ) ! Shift centre of mass to the origin

    DO i = 1, n-1 ! Loop to confirm bond lengths
       diff_sq = SUM ( (r(:,i)-r(:,i+1))**2 ) - bond**2
       IF ( ABS(diff_sq) > tol ) WRITE ( unit=error_unit, fmt='(a,2i15,f15.8)' ) 'Bond length warning ', i, i+1, diff_sq
    END DO ! End loop to confirm bond lengths

  END SUBROUTINE chain_positions

  SUBROUTINE chain_velocities ( temperature, constraints )
    USE maths_module, ONLY : random_normals, cross_product, solve, outer_product
    IMPLICIT NONE
    REAL,    INTENT(in) :: temperature ! reduced temperature
    LOGICAL, INTENT(in) :: constraints ! Option to constrain velocities relative to bonds

    ! Chooses velocities from Maxwell-Boltzmann (Gaussian) distribution
    ! For simplicity, we just pick each atom velocity randomly and
    ! apply bond constraints (if required) afterwards
    ! In between, we take steps to remove linear and angular momentum
    ! since the configuration is intended to be used in MD simulations without periodic boundaries
    ! in which case both these quantities are conserved
    ! NB there is at present no check for a singular inertia tensor in the angular momentum fix!
    ! We assume centre of mass is already at the origin
    ! We assume unit molecular mass and employ Lennard-Jones units

    REAL                 :: temp
    REAL, DIMENSION(3)   :: v_cm, r_cm, ang_mom, ang_vel
    REAL, DIMENSION(3,3) :: inertia
    INTEGER              :: i, xyz
    REAL, PARAMETER      :: tol = 1.e-6

    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Chain velocities at temperature', temperature

    ! Confirm centre-of-mass is at origin
    r_cm(:) = SUM ( r(:,:), dim=2 ) / REAL ( n ) ! Compute centre of mass
    IF ( ANY ( ABS ( r_cm ) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.8)' ) 'Centre of mass error', r_cm
       STOP 'Error in initialize_chain_velocities'
    END IF

    CALL random_normals ( 0.0, SQRT(temperature), v ) ! Choose 3N random velocities

    ! Compute and remove total momentum
    v_cm(:) = SUM ( v(:,:), dim=2 ) / REAL ( n )            ! Compute centre of mass velocity
    v(:,:)  = v(:,:) - SPREAD ( v_cm(:), dim=2, ncopies=n ) ! Set net momentum to zero

    ! Compute total angular momentum and moment of inertia tensor
    ang_mom = 0.0
    inertia = 0.0
    DO i = 1, n
       ang_mom = ang_mom + cross_product ( r(:,i), v(:,i) )
       inertia = inertia - outer_product ( r(:,i), r(:,i) )
       FORALL ( xyz=1:3 ) inertia(xyz,xyz) = inertia(xyz,xyz) + DOT_PRODUCT ( r(:,i), r(:,i) )
    END DO

    ! Solve linear system to get angular velocity
    ang_vel = solve ( inertia, ang_mom )

    ! Remove angular momentum
    DO i = 1, n
       v(:,i) = v(:,i) - cross_product ( ang_vel, r(:,i) )
    END DO

    IF ( constraints ) THEN

       ! Apply bond constraints (which should not introduce linear or angular momentum)
       WRITE ( unit=output_unit, fmt='(a)' ) 'Applying velocity constraints relative to bonds'
       CALL rattle_b

       ! Scale velocities to get correct temperature
       ! Number of degrees of freedom is 3*n - (n-1) bonds - 6 for angular and linear momentum
       temp   = SUM(v(:,:)**2) / REAL ( 3*n - (n-1) - 6 )
       v(:,:) = v(:,:) * SQRT ( temperature / temp )

    ELSE

       ! Scale velocities to get correct temperature
       ! Number of degrees of freedom is 3*n - 6 for angular and linear momentum
       temp   = SUM(v(:,:)**2) / REAL ( 3*n - 6 )
       v(:,:) = v(:,:) * SQRT ( temperature / temp )

    END IF

    ! Final check on angular and linear momenta
    v_cm    = 0.0
    ang_mom = 0.0
    DO i = 1, n
       v_cm    = v_cm + v(:,i)
       ang_mom = ang_mom + cross_product ( r(:,i), v(:,i) )
    END DO

    IF ( ANY ( ABS ( v_cm ) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.8)' ) 'Linear momentum error', v_cm
       STOP 'Error in chain_velocities'
    END IF

    IF ( ANY ( ABS ( ang_mom ) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.8)' ) 'Angular momentum error', ang_mom
       STOP 'Error in chain_velocities'
    END IF

  END SUBROUTINE chain_velocities

  SUBROUTINE rattle_b ! A version of velocity Verlet constraint algorithm
    IMPLICIT NONE

    ! This subroutine iteratively adjusts the velocities stored in the array v
    ! to satisfy the time derivatives of the bond constraints

    REAL, DIMENSION(3) :: rij, vij, dv
    LOGICAL            :: done
    REAL               :: g
    INTEGER            :: i, j, iter
    REAL,    PARAMETER :: tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 500
    LOGICAL, DIMENSION(:), ALLOCATABLE :: move, moved

    ALLOCATE ( move(n), moved(n) )

    iter     = 0
    done     = .FALSE.
    moved(:) = .TRUE.

    DO ! Iterative loop until done

       IF ( done ) EXIT

       done    = .TRUE.
       move(:) = .FALSE.

       DO i = 1, n-1 ! Loop over constraints
          j = i + 1 ! Partner atom for this constraint

          IF ( moved(i) .OR. moved(j) ) THEN ! Test whether need to re-examine ij
             vij = v(:,i) - v(:,j)
             rij = r(:,i) - r(:,j)

             ! In the following formulae, inverse masses are all unity
             g = -0.5*DOT_PRODUCT ( rij, vij ) / DOT_PRODUCT ( rij, rij )

             IF ( ABS ( g ) > tol ) THEN ! Test whether constraint already satisfied

                dv      = rij * g     ! Velocity adjustment
                v(:,i)  = v(:,i) + dv ! Adjust velocity i
                v(:,j)  = v(:,j) - dv ! Adjust velocity j
                move(i) = .TRUE.      ! Flag that we moved i
                move(j) = .TRUE.      ! Flag that we moved j
                done    = .FALSE.     ! Flag that we moved something

             END IF ! End test whether constraint already satisfied

          END IF ! End test whether need to re-examine ij

       END DO ! End loop over constraints

       ! Prepare for next iteration
       moved(:) = move(:)
       iter     = iter + 1
       IF ( iter > iter_max ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations', iter, iter_max
          STOP 'Error in rattle_b'
       END IF

    END DO ! End iterative loop until done

    DEALLOCATE ( move, moved )

  END SUBROUTINE rattle_b

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

  FUNCTION chain_overlap ( i, j1, j2 ) RESULT (overlap )
    IMPLICIT NONE
    LOGICAL             :: overlap ! Returns a flag indicating overlap with any j<i
    INTEGER, INTENT(in) :: i       ! Index of atom/molecule to be checked
    INTEGER, INTENT(in) :: j1      ! First j index
    INTEGER, INTENT(in) :: j2      ! Last j index

    ! This routine checks for overlaps of atoms
    ! NO box, NO periodic boundary conditions

    INTEGER            :: j
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq

    overlap = .FALSE.
    IF ( j2 < j1 ) RETURN

    DO j = j1, j2
       rij    = r(:,i) - r(:,j)
       rij_sq = SUM ( rij**2 )

       IF ( rij_sq < 1.0 ) THEN
          overlap = .TRUE.
          RETURN
       END IF

    END DO

  END FUNCTION chain_overlap

END MODULE initialize_module
