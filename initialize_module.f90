! initialize_module.f90
! Routines to initialize configurations and velocities
MODULE initialize_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_positions, initialize_orientations, initialize_velocities, initialize_angular_velocities
  PUBLIC :: initialize_chain_lattice, initialize_chain_random, initialize_chain_velocities
  PUBLIC :: n, r, v, e, w

  INTEGER                           :: n
  REAL, DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,:)
  REAL, DIMENSION(:,:), ALLOCATABLE :: v ! velocities (3,:)
  REAL, DIMENSION(:,:), ALLOCATABLE :: e ! orientations (3,:)
  REAL, DIMENSION(:,:), ALLOCATABLE :: w ! angular velocities (3,:)

CONTAINS

  SUBROUTINE initialize_positions ( nc )
    INTEGER, INTENT(in) :: nc ! number of unit cells in each coordinate direction

    ! Sets up the fcc lattice
    ! Four molecules per unit cell
    ! simulation box is a unit cube centred at the origin

    REAL, DIMENSION(3,4), PARAMETER :: r0 = RESHAPE ( [ &
         & 0.25, 0.25, 0.25,  0.75, 0.75, 0.25,  &
         & 0.25, 0.75, 0.75,  0.75, 0.25, 0.75 ],[3,4] ) ! positions in unit cell
    REAL, DIMENSION(3) :: rcm ! centre of mass
    INTEGER     ix, iy, iz, a, i

    WRITE(*,'(''Initializing positions on fcc lattice'')')
    
    IF ( n /= 4 * nc ** 3 ) STOP 'n, nc mismatch in init_positions'

    ALLOCATE ( r(3,n) )

    i = 0

    ! Begin triple loop over unit cell indices
    DO iz = 0, nc-1
       DO iy = 0, nc-1
          DO ix = 0, nc-1

             DO a = 1, 4 ! Begin loop over atoms in unit cell
                i = i + 1
                r(:,i) = r0(:,a) + REAL ( [ ix, iy, iz ] )
             END DO ! End loop over atoms in unit cell

          END DO
       END DO
    END DO
    ! End triple loop over unit cell indices

    r(:,:) = r(:,:)  / REAL ( nc ) ! Scale positions into unit cell
    rcm = SUM ( r, dim=2 ) / real(n)
    r(:,:) = r(:,:) - SPREAD ( rcm(:), dim=2, ncopies=n ) ! shift centre of mass to the origin
  END SUBROUTINE initialize_positions

  SUBROUTINE initialize_orientations ( nc )
    INTEGER, INTENT(in) :: nc ! number of unit cells in each coordinate direction

    ! Sets up the alpha-fcc lattice for linear molecules
    ! Four molecules per unit cell

    REAL, PARAMETER :: rroot3 = 1.0 / SQRT ( 3.0 )

    REAL, DIMENSION(3,4), PARAMETER :: e0 = RESHAPE (  rroot3*[ &
         &  1.0,  1.0,  1.0,    1.0, -1.0, -1.0,  &
         & -1.0,  1.0, -1.0,   -1.0, -1.0,  1.0 ],[3,4] ) ! orientations in unit cell

    INTEGER     k

    WRITE(*,'(''Initializing orientations on alpha-fcc lattice'')')

    IF ( n /= 4 * nc ** 3 ) STOP 'n, nc mismatch in initialize_orientations'

    ALLOCATE ( e(3,n) )

    ! Begin loop over unit cells (4 molecules per unit cell)
    DO k = 0, n-4, 4
       e(:,k+1:k+4) = e0(:,1:4) ! copy unit cell orientations
    END DO
    ! End loop over unit cells

  END SUBROUTINE initialize_orientations

  SUBROUTINE initialize_velocities ( temperature )
    USE utility_module, ONLY : random_normal
    REAL, INTENT(in) :: temperature ! reduced temperature

    ! Chooses velocities from Maxwell-Boltzmann (Gaussian) distribution
    ! We set the total momentum to zero afterwards
    ! we assume unit molecular mass and employ Lennard-Jones units
    ! property                  units
    ! energy                    epsilon ( = 1 )
    ! molecular mass            m ( = 1 )
    ! velocity v                sqrt(epsilon/m)
    ! angular velocity w        sqrt(epsilon/m*sigma**2)
    ! moment of inertia         m*sigma**2

    REAL               :: v_rms ! root-mean-square velocity component
    REAL, DIMENSION(3) :: v_sum ! total momentum
    INTEGER            :: i, k

    WRITE(*,'(''Initializing velocities at temperature'',t40,f15.5)') temperature

    ALLOCATE ( v(3,n) )

    v_rms = SQRT ( temperature )

    DO i = 1, n
       DO k = 1, 3
          v(k,i) = random_normal ( 0.0, v_rms )
       END DO
    END DO

    ! Compute and remove total momentum
    v_sum(:) = SUM ( v(:,:), dim=2 )
    v_sum(:) = v_sum(:) / REAL ( n )
    DO k = 1, 3
       v(k,:) = v(k,:) - v_sum(k)
    END DO

  END SUBROUTINE initialize_velocities

  SUBROUTINE initialize_angular_velocities ( temperature, inertia )
    USE utility_module, ONLY : random_perpendicular_vector
    REAL, INTENT(in) :: temperature ! reduced temperature
    REAL, INTENT(in) :: inertia     ! reduced moment of inertia

    ! this routine is specific to linear molecules
    ! it chooses the direction of the angular velocity randomly
    ! but perpendicular to the molecular axis.
    ! the square of the magnitude of the angular velocity
    ! is chosen from an exponential distribution
    ! there is no attempt to set the total angular momentum to zero

    REAL     ::   w_sq, w_sq_mean, zeta
    INTEGER   ::  i

    WRITE(*,'(''Initializing angular velocities at temperature'',t40,f15.5)') temperature
    WRITE(*,'(''Moment of inertia'',                             t40,f15.5)') inertia

    ALLOCATE ( w(3,n) )
    IF ( .NOT. ALLOCATED ( e ) ) STOP 'Allocation error in init_angular_velocities'
    IF ( SIZE(e,dim=2) /= n ) STOP 'Array mismatch in init_angular_velocities'

    w_sq_mean = 2.0 * temperature / inertia

    DO i = 1, n ! Begin loop over molecules
       w(:,i) = random_perpendicular_vector ( e(:,i) ) ! set direction of the angular velocity
       CALL RANDOM_NUMBER ( zeta )
       w_sq   = - w_sq_mean * LOG ( zeta ) ! squared magnitude of angular velocity
       w(:,i) = w(:,i) * SQRT ( w_sq )
    END DO

  END SUBROUTINE initialize_angular_velocities

  SUBROUTINE initialize_chain_lattice

    ! Sets up the fcc lattice, four molecules per unit cell
    ! unit cell is a unit cube, nearest neighbour distance is sqrt(0.5)
    ! results are then scaled to give unit bond length

    REAL, DIMENSION(3,4), PARAMETER :: r0 = RESHAPE ( [ &
         & 0.25, 0.25, 0.25, &
         & 0.25, 0.75, 0.75, &
         & 0.75, 0.75, 0.25, &
         & 0.75, 0.25, 0.75 ],[3,4] ) ! positions in unit cell

    INTEGER, DIMENSION(4,2,2) :: atoms_mid, atoms_inp, atoms_out
    INTEGER, DIMENSION(4)     :: atoms
    REAL                      :: rsq
    REAL,    dimension(3)     :: rcm ! centre of mass
    INTEGER                   :: nc, ix, iy, iz ! unit cell indices
    INTEGER                   :: a, i, j, x_direction, y_direction, plane_count
    INTEGER, DIMENSION(2)     :: istart, istop, istep
    REAL, PARAMETER           :: tol = 1.0e-9

    WRITE(*,'(''Initializing positions on fcc lattice'')')

    nc = NINT ( ( REAL(n)/4.0 )**(1.0/3.0) )
    IF ( n /= 4 * nc ** 3 ) STOP 'n, nc mismatch in initialize_chain_lattice'

    ALLOCATE ( r(3,n) )

    ! define sequences of atoms through unit cells
    atoms_inp(:,1,1) = [1,4,2,3]
    atoms_mid(:,1,1) = [1,4,2,3]
    atoms_out(:,1,1) = [1,2,3,4]
    atoms_inp(:,2,1) = [3,4,1,2]
    atoms_mid(:,2,1) = [4,1,3,2]
    atoms_out(:,2,1) = [4,1,3,2]
    atoms_inp(:,1,2) = [1,2,3,4]
    atoms_mid(:,1,2) = [2,3,1,4]
    atoms_out(:,1,2) = [2,3,1,4]
    atoms_inp(:,2,2) = [3,2,4,1]
    atoms_mid(:,2,2) = [3,2,4,1]
    atoms_out(:,2,2) = [3,4,1,2]

    ! forward and reverse options to traverse cells in xy plane
    istart = [ 1, nc ]
    istop  = [ nc, 1 ]
    istep  = [ 1, -1 ]

    ! Begin triple loop over unit cell indices
    i = 0
    DO iz = 1, nc
       y_direction = MOD(iz,2) + 1 ! alternate scans across xy plane
       plane_count = 0
       DO iy = istart(y_direction), istop(y_direction), istep(y_direction)
          x_direction = MOD(iz+iy,2) + 1 ! alternate scans across x rows
          DO ix = istart(x_direction), istop(x_direction), istep(x_direction)
             plane_count = plane_count + 1
             IF ( plane_count == 1 ) THEN
                atoms = atoms_inp(:,x_direction,y_direction) ! first corner of plane
             ELSE IF ( plane_count == nc**2 ) THEN
                atoms = atoms_out(:,x_direction,y_direction) ! last corner of plane
             ELSE
                atoms = atoms_mid(:,x_direction,y_direction) ! any other cell in plane
             END IF

             DO a = 1, 4 ! Begin loop over atoms in unit cell
                i = i+1
                r(:,i) = r0(:,atoms(a)) + REAL ( [ ix-1, iy-1, iz-1 ] )
             END DO ! End loop over atoms in unit cell

          END DO
       END DO
    END DO
    ! End triple loop over unit cell indices

    rcm = SUM ( r, dim=2 ) / real(n)
    r(:,:) = r(:,:) - SPREAD ( rcm, dim=2, ncopies = n )  ! shift centre of box to the origin
    r(:,:) = r(:,:) / SQRT(0.5)     ! scale to give unit bond length

    ! Double loop to confirm unit bond lengths and no overlaps
    DO i = 1, n-1
       DO j = i+1, n
          rsq =  SUM((r(:,i)-r(:,j))**2)
          IF ( j == i+1 ) THEN
             IF ( ABS(rsq-1.0) > tol ) PRINT *, 'Bond length warning ', i, i+1, rsq
          ELSE
             IF ( ABS(rsq) < 1.0 - tol ) PRINT *, 'Overlap warning ', i, i+1, rsq
          END IF
       END DO
    END DO

  END SUBROUTINE initialize_chain_lattice

  SUBROUTINE initialize_chain_random
    USE utility_module, ONLY : random_bond => random_orientation_vector_alt1

    ! Chooses chain positions randomly, at unit bond length, avoiding overlap

    REAL               :: rsq
    REAL, DIMENSION(3) :: d, rcm
    INTEGER            :: i, j, iter
    REAL, PARAMETER    :: tol = 1.0e-9
    INTEGER, parameter :: iter_max = 1000
    LOGICAL            :: overlap

    WRITE(*,'(''Initializing positions randomly'')')

    ALLOCATE ( r(3,n) )

    ! First atom at origin
    r(:,1) = [0.0,0.0,0.0]

    ! Second atom at random position
    CALL random_bond ( d )
    r(:,2) = r(:,1) + d

    DO i = 3, n ! Begin loop over atom indices

       iter = 0
       DO ! Loop until non-overlapping position found (n must not be too large!)
          CALL random_bond ( d )
          r(:,i) = r(:,i-1) + d
          overlap = .FALSE.

          DO j = 1, i-2 ! Loop to check for overlap
             overlap = SUM ( (r(:,i)-r(:,j))**2 ) < 1.0
             IF ( overlap ) EXIT
          END DO ! End loop to check for overlap

          IF ( .NOT. overlap ) EXIT

          iter = iter + 1
          IF ( iter > iter_max ) STOP 'too many iterations in initialize_chain_random'
          
       END DO ! End loop until non-overlapping position found

    END DO ! End loop over atom indices

    rcm = SUM ( r, dim=2 ) / REAL(n)
    r(:,:) = r(:,:) - SPREAD ( rcm, dim=2, ncopies = n )  ! shift centre of box to the origin

    ! Double loop to confirm unit bond lengths and no overlaps
    DO i = 1, n-1
       DO j = i+1, n
          rsq =  SUM((r(:,i)-r(:,j))**2)
          IF ( j == i+1 ) THEN
             IF ( ABS(rsq-1.0) > tol ) PRINT *, 'Bond length warning ', i, j, rsq
          ELSE
             IF ( ABS(rsq) < 1.0 - tol ) PRINT *, 'Overlap warning ', i, j, rsq
          END IF
       END DO
    END DO

  END SUBROUTINE initialize_chain_random

  SUBROUTINE initialize_chain_velocities ( temperature )
    USE utility_module, ONLY : random_normal, cross_product, random_perpendicular_vector
    REAL, INTENT(in) :: temperature ! reduced temperature

    ! Chooses velocities from Maxwell-Boltzmann (Gaussian) distribution
    ! For simplicity, we just pick each atom velocity perpendicular to bonds
    ! which makes it completely unsuitable as a general-purpose thermostat
    ! We assume there is no need to apply periodic boundary conditions
    ! We set the total momentum to zero afterwards
    ! we assume unit molecular mass and employ Lennard-Jones units
    ! property                  units
    ! energy                    epsilon ( = 1 )
    ! molecular mass            m ( = 1 )
    ! velocity v                sqrt(epsilon/m)

    REAL               :: temp, v1, v2
    REAL, DIMENSION(3) :: v_vec ! total momentum
    INTEGER            :: i, k

    WRITE(*,'(''Initializing velocities at temperature'',t40,f15.5)') temperature

    ALLOCATE ( v(3,n) )

    DO i = 1, n
       IF ( i == 1 ) THEN
          v_vec = random_perpendicular_vector ( r(:,1)-r(:,2) )
          v1 = random_normal ( 0.0, 1.0 )
          v2 = random_normal ( 0.0, 1.0 )
          v(:,i) = v_vec * SQRT ( v1**2 + v2**2 )
       ELSE IF ( i == n ) THEN
          v_vec = random_perpendicular_vector ( r(:,n)-r(:,n-1) )
          v1 = random_normal ( 0.0, 1.0 )
          v2 = random_normal ( 0.0, 1.0 )
          v(:,i) = v_vec * SQRT ( v1**2 + v2**2 )
       ELSE
          v_vec = cross_product ( r(:,i)-r(:,i-1), r(:,i)-r(:,i+1) )
          temp = SQRT ( SUM ( v_vec**2) )
          v_vec = v_vec / temp ! unit vector
          v(:,i) = v_vec * random_normal ( 0.0, 1.0 )
       END IF
    END DO

    ! Compute and remove total momentum
    v_vec(:) = SUM ( v(:,:), dim=2 )
    v_vec(:) = v_vec(:) / REAL ( n )
    DO k = 1, 3
       v(k,:) = v(k,:) - v_vec(k)
    END DO

    ! number of degrees of freedom is 3*(n-1) - (n-1) constraints
    temp = SUM(v**2) / REAL ( 2*(n-1) )
    v = v * SQRT ( temperature / temp )

  END SUBROUTINE initialize_chain_velocities

END MODULE initialize_module
