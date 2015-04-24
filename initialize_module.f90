! initialize_module.f90
! Routines to initialize configurations and velocities
MODULE initialize_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_positions, initialize_orientations, initialize_velocities, initialize_angular_velocities
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

    INTEGER     ix, iy, iz, a, m

    WRITE(*,'(''Initializing positions on fcc lattice'')')
    
    IF ( n /= 4 * nc ** 3 ) STOP 'n, nc mismatch in init_positions'

    ALLOCATE ( r(3,n) )

    m = 0

    ! Begin triple loop over unit cell indices
    DO iz = 0, nc-1
       DO iy = 0, nc-1
          DO ix = 0, nc-1

             DO a = 1, 4 ! Begin loop over atoms in unit cell
                r(:,a+m) = r0(:,a) + REAL ( [ ix, iy, iz ] )
             END DO ! End loop over atoms in unit cell

             m = m + 4

          END DO
       END DO
    END DO
    ! End triple loop over unit cell indices

    r(:,:) = r(:,:)  / REAL ( nc ) ! Scale positions into unit box
    r(:,:) = r(:,:) - 0.5          ! shift centre of box to the origin
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

    IF ( n /= 4 * nc ** 3 ) STOP 'n, nc mismatch in init_orientations'

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

END MODULE initialize_module
