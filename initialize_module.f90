! initialize_module.f90
! Routines to initialize configurations and velocities
MODULE initialize_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: allocate_arrays, deallocate_arrays
  PUBLIC :: initialize_positions_lattice, initialize_orientations_lattice
  PUBLIC :: initialize_positions_random, initialize_orientations_random
  PUBLIC :: initialize_velocities, initialize_angular_velocities
  PUBLIC :: initialize_chain_lattice, initialize_chain_random, initialize_chain_velocities

  ! Public data
  INTEGER,                           PUBLIC :: n ! Number of atoms
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v ! Velocities (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: e ! Orientations (3,n) or (0:3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: w ! Angular velocities (3,n)

  ! Private data
  LOGICAL, DIMENSION(:), ALLOCATABLE :: move, moved
  
CONTAINS

  SUBROUTINE allocate_arrays ( quaternions )
    LOGICAL, INTENT(in) :: quaternions

    ALLOCATE ( r(3,n), v(3,n), w(3,n), move(n), moved(n) )

    IF ( quaternions ) THEN
       ALLOCATE ( e(0:3,n) )
    ELSE
       ALLOCATE ( e(3,n) )
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays

    DEALLOCATE ( r, v, w, e, move, moved )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE initialize_positions_lattice

    ! Sets up the fcc lattice: four molecules per unit cell
    ! Initially, employ a unit cell of unit length
    ! Afterwards, rescale into simulation box of unit length centred at the origin

    REAL, DIMENSION(3,4), PARAMETER :: r_fcc = RESHAPE ( [ &
         & 0.25, 0.25, 0.25, &
         & 0.25, 0.75, 0.75, &
         & 0.75, 0.75, 0.25, &
         & 0.75, 0.25, 0.75 ], [3,4] ) ! Positions in unit cell

    REAL, DIMENSION(3) :: r_cm
    INTEGER            :: nc, ix, iy, iz, a, i

    WRITE ( unit=output_unit, fmt='(a)' ) 'Close-packed lattice positions'

    nc = NINT ( REAL(n/4) ** (1.0/3.0) )
    IF ( n /= 4 * nc ** 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'n, nc mismatch ', n, 4 * nc ** 3
       STOP 'Error in initialize_positions_lattice'
    END IF

    IF ( .NOT. ALLOCATED ( r ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array r is not allocated'
       STOP 'Error in initialize_positions_lattice'
    END IF

    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'Error in shape of r', SHAPE(r), 3, n
       STOP 'Error in initialize_positions_lattice'
    END IF

    i = 0

    ! Begin triple loop over unit cell indices
    DO iz = 0, nc-1
       DO iy = 0, nc-1
          DO ix = 0, nc-1

             DO a = 1, 4 ! Begin loop over atoms in unit cell
                i = i + 1
                r(:,i) = r_fcc(:,a) + REAL ( [ ix, iy, iz ] )
             END DO ! End loop over atoms in unit cell

          END DO
       END DO
    END DO
    ! End triple loop over unit cell indices

    r(:,:)  = r(:,:)  / REAL ( nc )                         ! Scale positions into unit cell
    r_cm(:) = SUM ( r, dim=2 ) / REAL(n)                    ! Compute centre of mass position
    r(:,:)  = r(:,:) - SPREAD ( r_cm(:), dim=2, ncopies=n ) ! Shift centre of mass to the origin

  END SUBROUTINE initialize_positions_lattice
 
  SUBROUTINE initialize_positions_random

    ! Places atoms at random positions
    ! Unlikely to be useful, unless the interaction potential is soft
    ! Simulation box is a unit cube centred at the origin

    REAL, DIMENSION(3) :: r_cm

    WRITE ( unit=output_unit, fmt='(a)' ) 'Random positions'

    IF ( .NOT. ALLOCATED ( r ) ) THEN
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array r is not allocated'
       STOP 'Error in initialize_positions_random'
    END IF

    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'Error in shape of r', SHAPE(r), 3, n
       STOP 'Error in initialize_positions_random'
    END IF

    CALL RANDOM_NUMBER ( r(:,:) )                           ! All in range (0,1)
    r_cm(:) = SUM ( r(:,:), dim=2 ) / REAL(n)               ! Compute centre of mass position
    r(:,:)  = r(:,:) - SPREAD ( r_cm(:), dim=2, ncopies=n ) ! Shift centre of mass to the origin
    
  END SUBROUTINE initialize_positions_random
 
  SUBROUTINE initialize_orientations_lattice

    ! Sets up the alpha-fcc lattice (4 per unit cell) for linear molecules
    ! Sets up perfectly aligned configuration for nonlinear molecules

    REAL, DIMENSION(3,4), PARAMETER :: e_fcc = RESHAPE ( SQRT(1.0/3.0) * [ &
         &  1.0,  1.0,  1.0,  &
         &  1.0, -1.0, -1.0,  &
         & -1.0,  1.0, -1.0,  &
         & -1.0, -1.0,  1.0  ], [3,4] ) ! Orientations in unit cell

    INTEGER :: nc, k

    WRITE ( unit=output_unit, fmt='(a)' ) 'Regular lattice of orientations'

    nc = NINT ( REAL(n/4) ** (1.0/3.0) )
    IF ( n /= 4 * nc ** 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'n, nc mismatch ', n, 4 * nc ** 3
       STOP 'Error in initialize_orientations_lattice'
    END IF

    IF ( .NOT. ALLOCATED ( e ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array e is not allocated'
       STOP 'Error in initialize_orientations_lattice'
    END IF

    IF ( LBOUND(e,dim=1) /= 0 .AND. LBOUND(e,dim=1) /= 1 ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Array e lower bound mismatch ', LBOUND(e,dim=1)
       STOP 'Error in initialize_orientations_lattice'
    END IF

    IF ( UBOUND(e,dim=1) /= 3 ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Array e upper bound mismatch ', UBOUND(e,dim=1)
       STOP 'Error in initialize_orientations_lattice'
    END IF

    IF ( SIZE(e,dim=2) /= n   ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array e bounds mismatch ', SIZE(e,dim=2), n
       STOP 'Error in initialize_orientations_lattice'
    END IF

    IF ( LBOUND(e,dim=1) == 0 ) THEN

       DO k = 1, n
          e(:,k) = [1.0,0.0,0.0,0.0] ! For nonlinear molecules simply set all orientations the same
       END DO

    ELSE

       DO k = 0, n-4, 4 ! Loop over unit cells (4 molecules per unit cell)
          e(:,k+1:k+4) = e_fcc(:,1:4) ! Copy unit cell orientations
       END DO ! End loop over unit cells

    END IF

  END SUBROUTINE initialize_orientations_lattice

  SUBROUTINE initialize_orientations_random
    USE maths_module, ONLY : random_quaternion, random_vector

    ! Sets up random orientations for linear or nonlinear molecules

    INTEGER :: i

    WRITE ( unit=output_unit, fmt='(a)' ) 'Random orientations'

    IF ( .NOT. ALLOCATED ( e ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array e is not allocated'
       STOP 'Error in initialize_orientations_random'
    END IF

    IF ( LBOUND(e,dim=1) /= 0 .AND. LBOUND(e,dim=1) /= 1 ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Array e lower bound mismatch ', LBOUND(e,dim=1)
       STOP 'Error in initialize_orientations_random'
    END IF

    IF ( UBOUND(e,dim=1) /= 3 ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Array e upper bound mismatch ', UBOUND(e,dim=1)
       STOP 'Error in initialize_orientations_random'
    END IF

    IF ( SIZE(e,dim=2) /= n   ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array e bounds mismatch ', SIZE(e,dim=2), n
       STOP 'Error in initialize_orientations_random'
    END IF

    IF ( LBOUND(e,dim=1) == 0 ) THEN

       DO i = 1, n
          e(:,i) = random_quaternion ( )
       END DO

    ELSE

       DO i = 1, n
          e(:,i) = random_vector ( )
       END DO

    END IF

  END SUBROUTINE initialize_orientations_random

  SUBROUTINE initialize_velocities ( temperature )
    USE maths_module, ONLY : random_normals
    REAL, INTENT(in) :: temperature ! Reduced temperature

    ! Chooses velocities from Maxwell-Boltzmann (Gaussian) distribution
    ! We set the total momentum to zero afterwards
    ! We assume unit molecular mass and employ simulation (e.g. Lennard-Jones) units
    ! property                  units
    ! energy                    epsilon ( = 1 )
    ! molecular mass            m ( = 1 )
    ! velocity v                sqrt(epsilon/m)
    ! angular velocity w        sqrt(epsilon/m*sigma**2)
    ! moment of inertia         m*sigma**2

    REAL               :: v_rms ! Root-mean-square velocity component
    REAL, DIMENSION(3) :: v_cm  ! Centre of mass velocity

    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Velocities at temperature', temperature

    IF ( .NOT. ALLOCATED ( v ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array v is not allocated'
       STOP 'Error in initialize_velocities'
    END IF

    IF ( ANY ( SHAPE(v) /= [3,n] ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'Error in shape of v', SHAPE(v), 3, n
       STOP 'Error in initialize_velocities'
    END IF

    v_rms = SQRT ( temperature )
    CALL random_normals ( 0.0, v_rms, v )

    v_cm(:) = SUM ( v(:,:), dim=2 ) / REAL ( n )            ! Compute centre of mass velocity
    v(:,:)  = v(:,:) - SPREAD ( v_cm(:), dim=2, ncopies=n ) ! Set net momentum to zero

  END SUBROUTINE initialize_velocities

  SUBROUTINE initialize_angular_velocities ( temperature, inertia )
    USE maths_module, ONLY : random_perpendicular_vector, random_normals
    REAL, INTENT(in) :: temperature ! Reduced temperature
    REAL, INTENT(in) :: inertia     ! Reduced moment of inertia

    ! For linear molecules we choose the direction of the angular velocity
    ! randomly but perpendicular to the molecular axis.
    ! The square of the magnitude of the angular velocity
    ! is chosen from an exponential distribution
    ! For nonlinear molecules we choose all three components of angular velocity
    ! from a Gaussian distribution, assuming equal moments of inertia
    ! There is no attempt to set the total angular momentum to zero

    REAL    ::  w_sq, w_sq_mean, w_std_dev,  zeta
    INTEGER ::  i

    WRITE ( unit=output_unit, fmt='(a,t40,2f15.6)' ) 'Angular velocities at temperature, inertia', temperature, inertia

    IF ( .NOT. ALLOCATED ( w ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array w is not allocated'
       STOP 'Error in initialize_angular_velocities'
    END IF

    IF ( ANY ( SHAPE(w) /= [3,n] ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'Error in shape of w', SHAPE(r), 3, n
       STOP 'Error in initialize_angular_velocities'
    END IF

    IF ( .NOT. ALLOCATED ( e ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array e is not allocated'
       STOP 'Error in initialize_angular_velocities'
    END IF

    IF ( SIZE(e,dim=2) /= n ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array e bounds mismatch ', SIZE(e,dim=2), n
       STOP 'Error in initialize_angular_velocities'
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

  END SUBROUTINE initialize_angular_velocities

  SUBROUTINE initialize_chain_lattice

    ! Sets up the fcc lattice, four molecules per unit cell
    ! Unit cell is a unit cube, nearest neighbour distance is sqrt(0.5)
    ! Results are then scaled to give unit bond length

    ! The aim of the slightly complicated triple loop over unit cells
    ! with various options for atom positions within the unit cell
    ! is to ensure that successive atoms are exactly one bond length apart

    REAL, DIMENSION(3,4), PARAMETER :: r_fcc = RESHAPE ( [ &
         & 0.25, 0.25, 0.25, &
         & 0.25, 0.75, 0.75, &
         & 0.75, 0.75, 0.25, &
         & 0.75, 0.25, 0.75 ],[3,4] ) ! Positions in unit cell

    ! Pick out a particular sequence of atoms within cell, according to x_direction and y_direction
    INTEGER, DIMENSION(4,2,2), PARAMETER :: atoms_mid = RESHAPE ( [ 1,4,2,3, 4,1,3,2, 2,3,1,4, 3,2,4,1 ], [4,2,2] )
    INTEGER, DIMENSION(4,2,2), PARAMETER :: atoms_inp = RESHAPE ( [ 1,4,2,3, 3,4,1,2, 1,2,3,4, 3,2,4,1 ], [4,2,2] )
    INTEGER, DIMENSION(4,2,2), PARAMETER :: atoms_out = RESHAPE ( [ 1,2,3,4, 4,1,3,2, 2,3,1,4, 3,4,1,2 ], [4,2,2] )

    INTEGER, DIMENSION(4)     :: atoms
    REAL                      :: r_sq
    REAL,    DIMENSION(3)     :: r_cm
    INTEGER                   :: nc, ix, iy, iz
    INTEGER                   :: a, i, j, x_direction, y_direction, plane_count
    INTEGER, DIMENSION(2)     :: istart, istop, istep
    REAL, PARAMETER           :: tol = 1.0e-9

    WRITE ( unit=output_unit, fmt='(a)' ) 'Chain, close-packed atoms surrounded by vacuum'

    nc = NINT ( ( REAL(n)/4.0 )**(1.0/3.0) )
    IF ( n /= 4 * nc ** 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'n, nc mismatch ', n, 4 * nc ** 3
       STOP 'Error in initialize_chain_lattice'
    END IF

    IF ( .NOT. ALLOCATED ( r ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array r is not allocated'
       STOP 'Error in initialize_chain_lattice'
    END IF

    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'Error in shape of r', SHAPE(r), 3, n
       STOP 'Error in initialize_chain_lattice'
    END IF

    ! Forward and reverse options to traverse cells in xy plane
    istart = [ 1, nc ]
    istop  = [ nc, 1 ]
    istep  = [ 1, -1 ]

    ! Begin triple loop over unit cell indices
    i = 0
    DO iz = 1, nc
       y_direction = MOD(iz,2) + 1 ! Alternate scans across xy plane
       plane_count = 0
       DO iy = istart(y_direction), istop(y_direction), istep(y_direction)
          x_direction = MOD(iz+iy,2) + 1 ! Alternate scans along x rows
          DO ix = istart(x_direction), istop(x_direction), istep(x_direction)

             plane_count = plane_count + 1
             IF ( plane_count == 1 ) THEN
                atoms = atoms_inp(:,x_direction,y_direction) ! First corner of plane
             ELSE IF ( plane_count == nc**2 ) THEN
                atoms = atoms_out(:,x_direction,y_direction) ! Last corner of plane
             ELSE
                atoms = atoms_mid(:,x_direction,y_direction) ! Any other cell in plane
             END IF

             DO a = 1, 4 ! Begin loop over atoms in unit cell
                i = i+1
                r(:,i) = r_fcc(:,atoms(a)) + REAL ( [ ix-1, iy-1, iz-1 ] )
             END DO ! End loop over atoms in unit cell

          END DO
       END DO
    END DO
    ! End triple loop over unit cell indices

    r_cm(:) = SUM ( r(:,:), dim=2 ) / REAL(n)               ! Compute centre of mass position
    r(:,:)  = r(:,:) - SPREAD ( r_cm(:), dim=2, ncopies=n ) ! Shift centre of box to the origin
    r(:,:)  = r(:,:) / SQRT(0.5)                            ! Scale to give unit bond length

    DO i = 1, n-1 ! Loop to confirm unit bond lengths
       r_sq = SUM ( (r(:,i)-r(:,i+1))**2 )
       IF ( ABS(r_sq-1.0) > tol ) WRITE ( unit=error_unit, fmt='(a,2i15,f15.8)' ) 'Bond length warning ', i, i+1, r_sq
    END DO ! End loop to confirm unit bond lengths

    ! Double loop to confirm no overlaps
    DO i = 1, n-2
       DO j = i+2, n
          r_sq = SUM((r(:,i)-r(:,j))**2)
          IF ( ABS(r_sq) < 1.0 - tol ) WRITE ( unit=error_unit, fmt='(a,2i15,f15.8)' ) 'Overlap warning ', i, j, r_sq
       END DO
    END DO
    ! End double loop to confirm no overlaps

  END SUBROUTINE initialize_chain_lattice

  SUBROUTINE initialize_chain_random
    USE maths_module, ONLY : random_vector

    ! Chooses chain positions randomly, at unit bond length, avoiding overlap

    REAL               :: r_sq
    REAL, DIMENSION(3) :: r_cm
    INTEGER            :: i, j, iter
    LOGICAL            :: overlap

    REAL,    PARAMETER :: tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 1000

    WRITE ( unit=output_unit, fmt='(a)' ) 'Chain, bonds randomly oriented, avoiding overlaps'

    IF ( .NOT. ALLOCATED ( r ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array r is not allocated'
       STOP 'Error in initialize_chain_random'
    END IF
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'Error in shape of r', SHAPE(r), 3, n
       STOP 'Error in initialize_chain_random'
    END IF

    r(:,1) = [0.0,0.0,0.0]   ! First atom at origin
    r(:,2) = random_vector() ! Second atom at random position (unit bond length away)

    DO i = 3, n ! Begin loop over atom indices

       iter = 0
       DO ! Loop until non-overlapping position found (N must not be too large!)
          r(:,i) = r(:,i-1) + random_vector() ! Subsequent atoms randomly placed (unit bond length)
          overlap = .FALSE.

          DO j = 1, i-2 ! Loop to check for overlap
             overlap = SUM ( (r(:,i)-r(:,j))**2 ) < 1.0
             IF ( overlap ) EXIT
          END DO ! End loop to check for overlap

          IF ( .NOT. overlap ) EXIT

          iter = iter + 1
          IF ( iter > iter_max ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations ', iter, iter_max
             STOP 'Error in initialize_chain_random'
          END IF

       END DO ! End loop until non-overlapping position found

    END DO ! End loop over atom indices

    r_cm(:) = SUM ( r(:,:), dim=2 ) / REAL(n)               ! Compute centre of mass positions
    r(:,:)  = r(:,:) - SPREAD ( r_cm(:), dim=2, ncopies=n ) ! Shift centre of mass to the origin

    DO i = 1, n-1 ! Loop to confirm unit bond lengths
       r_sq = SUM ( (r(:,i)-r(:,i+1))**2 )
       IF ( ABS(r_sq-1.0) > tol ) WRITE ( unit=error_unit, fmt='(a,2i15,f15.8)' ) 'Bond length warning ', i, i+1, r_sq
    END DO ! End loop to confirm unit bond lengths

    ! Double loop to confirm no overlaps
    DO i = 1, n-2
       DO j = i+2, n
          r_sq = SUM((r(:,i)-r(:,j))**2)
          IF ( ABS(r_sq) < 1.0 - tol ) WRITE ( unit=error_unit, fmt='(a,2i15,f15.8)' ) 'Overlap warning ', i, j, r_sq
       END DO
    END DO
    ! End double loop to confirm no overlaps

  END SUBROUTINE initialize_chain_random

  SUBROUTINE initialize_chain_velocities ( temperature )
    USE maths_module, ONLY : random_normals, cross_product, solve, outer_product
    REAL, INTENT(in) :: temperature ! reduced temperature

    ! Chooses velocities from Maxwell-Boltzmann (Gaussian) distribution
    ! For simplicity, we just pick each atom velocity randomly and
    ! apply bond constraints afterwards
    ! In between, we take steps to remove linear and angular momentum
    ! since the configuration may be used in MD simulations without periodic boundaries
    ! in which case both these quantities are conserved
    ! We assume centre of mass is already at the origin
    ! We assume unit molecular mass and employ Lennard-Jones units
    ! property                  units
    ! energy                    epsilon ( = 1 )
    ! molecular mass            m ( = 1 )
    ! velocity v                sqrt(epsilon/m)

    REAL                 :: temp
    REAL, DIMENSION(3)   :: v_cm, ang_mom, ang_vel
    REAL, DIMENSION(3,3) :: inertia
    INTEGER              :: i, xyz
    REAL, PARAMETER      :: tol = 1.e-6

    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Chain velocities at temperature', temperature

    IF ( .NOT. ALLOCATED ( v ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a)' ) 'Array v is not allocated'
       STOP 'Error in initialize_chain_velocities'
    END IF

    IF ( ANY ( SHAPE(v) /= [3,n] ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'Error in shape of v', SHAPE(v), 3, n
       STOP 'Error in initialize_chain_velocities'
    END IF

    CALL random_normals ( 0.0, SQRT(temperature), v ) ! Choose 3N random velocities

    ! Compute and remove total momentum
    v_cm(:) = SUM ( v(:,:), dim=2 ) / REAL ( n )            ! Compute centre of mass velocity
    v(:,:)  = v(:,:) - SPREAD ( v_cm(:), dim=2, ncopies=n ) ! Set net momentum to zero

    ! Compute total angular momentum and moment of inertia tensor
    ang_mom = 0.0
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

    ! Apply bond constraints (which should not introduce linear or angular momentum)
    CALL rattle_b

    ! Scale velocities to get correct temperature
    ! Number of degrees of freedom is 3*n - (n-1) bonds - 6 for angular and linear momentum
    temp   = SUM(v(:,:)**2) / REAL ( 3*n - (n-1) - 6 )
    v(:,:) = v(:,:) * SQRT ( temperature / temp )

    ! Check angular and linear momenta
    v_cm    = 0.0
    ang_mom = 0.0
    DO i = 1, n
       v_cm    = v_cm + v(:,i)
       ang_mom = ang_mom + cross_product ( r(:,i), v(:,i) )
    END DO

    IF ( ANY ( ABS ( v_cm ) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.8)' ) 'Linear momentum error', v_cm
       STOP 'Error in initialize_chain_velocities'
    END IF

    IF ( ANY ( ABS ( ang_mom ) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.8)' ) 'Angular momentum error', ang_mom
       STOP 'Error in initialize_chain_velocities'
    END IF

  END SUBROUTINE initialize_chain_velocities

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

  END SUBROUTINE rattle_b

END MODULE initialize_module
