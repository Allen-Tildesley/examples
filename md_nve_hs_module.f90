! md_nve_hs_module.f90
! Collisions and overlap for MD of hard spheres
MODULE md_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, r, v, coltime, partner, lt, ne, gt
  PUBLIC :: allocate_arrays, deallocate_arrays, update, overlap, collide

  INTEGER                              :: n       ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r, v    ! positions, velocities (3,n)
  REAL,    DIMENSION(:),   ALLOCATABLE :: coltime ! time to next collision (n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: partner ! collision partner (n)

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options

CONTAINS

  SUBROUTINE allocate_arrays
    ALLOCATE ( r(3,n), v(3,n), coltime(n), partner(n) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, v, coltime, partner )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE update ( i, j_range, box ) ! updates collision details for atom i
    INTEGER, INTENT(in) :: i, j_range   ! atom and set of collision partners
    REAL,    INTENT(in) :: box          ! simulation box length

    INTEGER            :: j, j1, j2
    REAL, DIMENSION(3) :: rij, vij
    REAL               :: rijsq, vijsq, bij, tij, discr

    SELECT CASE ( j_range )
    CASE ( lt ) ! j < i
       j1 = 1
       j2 = i-1
    CASE ( gt ) ! j > i
       j1 = i+1
       j2 = n
    CASE ( ne ) ! j /= i
       j1 = 1
       j2 = n
    END SELECT

    coltime(i) = HUGE(1.0)

    DO j = j1, j2

       IF ( i == j ) CYCLE

       rij(:) = r(:,i) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) )
       rij(:) = rij(:) * box ! now in sigma=1 units
       vij(:) = v(:,i) - v(:,j)
       bij  = DOT_PRODUCT ( rij, vij )

       IF ( bij < 0.0 ) THEN

          rijsq = SUM ( rij**2 )
          vijsq = SUM ( vij**2 )
          discr = bij ** 2 - vijsq * ( rijsq - 1.0 ) ! sigma**2 = 1.0

          IF ( discr > 0.0 ) THEN

             tij = ( -bij - SQRT ( discr ) ) / vijsq

             IF ( tij < coltime(i) ) THEN

                coltime(i) = tij
                partner(i) = j

             END IF

          END IF

       END IF

    END DO

  END SUBROUTINE update

  FUNCTION overlap ( box )      ! tests configuration for pair overlaps
    LOGICAL          :: overlap ! function result
    REAL, INTENT(in) :: box     ! simulation box length

    INTEGER            :: i, j
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq, box_sq, rij_mag
    REAL,    PARAMETER :: tol = 1.0e-4 

    overlap = .FALSE.
    box_sq  = box**2

    DO i = 1, n - 1
       DO j = i + 1, n

          rij(:) = r(:,i) - r(:,j)
          rij(:) = rij(:) - ANINT ( rij(:) )
          rij_sq = SUM ( rij**2 )  ! squared distance
          rij_sq = rij_sq * box_sq ! now in sigma=1 units

          IF ( rij_sq < 1.0 ) THEN
             rij_mag = SQRT(rij_sq)
             WRITE ( unit=error_unit, fmt='(a,2i5,f15.8)' ) 'Warning: i,j,rij = ', i, j, rij_mag
             IF ( ( 1.0 - rij_mag ) > tol ) overlap = .TRUE.
          END IF

       END DO
    END DO

  END FUNCTION overlap

  SUBROUTINE collide ( i, j, box, virial ) ! collision dynamics
    INTEGER, INTENT(in)  :: i, j           ! colliding atoms, assumed to be in contact
    REAL,    INTENT(in)  :: box            ! simulation box length
    REAL,    INTENT(out) :: virial         ! collision contribution to pressure

    REAL, DIMENSION(3) :: rij, vij
    REAL               :: factor

    rij(:) = r(:,i) - r(:,j)
    rij(:) = rij(:) - ANINT ( rij(:) ) ! separation vector
    rij(:) = rij(:) * box              ! now in sigma=1 units
    vij(:) = v(:,i) - v(:,j)           ! relative velocity

    factor = DOT_PRODUCT ( rij, vij )
    vij = - factor * rij

    v(:,i) = v(:,i) + vij
    v(:,j) = v(:,j) - vij
    virial = DOT_PRODUCT ( vij, rij ) / 3.0
  END SUBROUTINE collide

END MODULE md_module
