! md_nve_hs_module.f90
! Collisions and overlap for MD of hard spheres
MODULE md_nve_hs_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, v, coltime, partner, lt, ne, gt
  PUBLIC :: initialize, finalize, update, overlap, collide

  INTEGER                              :: n       ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r, v    ! positions, velocities (3,n)
  REAL,    DIMENSION(:),   ALLOCATABLE :: coltime ! time to next collision (n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: partner ! collision partner (n)

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options

CONTAINS

  SUBROUTINE initialize
    ALLOCATE ( r(3,n), v(3,n), coltime(n), partner(n) )
  END SUBROUTINE initialize

  SUBROUTINE finalize
    DEALLOCATE ( r, v, coltime, partner )
  END SUBROUTINE finalize

  SUBROUTINE update ( i, j_range, sigma_sq ) ! updates collision details for atom i
    INTEGER, INTENT(in) :: i, j_range
    REAL,    INTENT(in) :: sigma_sq

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
       vij(:) = v(:,i) - v(:,j)
       bij  = DOT_PRODUCT ( rij, vij )

       IF ( bij < 0.0 ) THEN

          rijsq = SUM ( rij**2 )
          vijsq = SUM ( vij**2 )
          discr = bij ** 2 - vijsq * ( rijsq - sigma_sq )

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

  FUNCTION overlap ( sigma_sq ) ! tests configuration for pair overlaps
    LOGICAL          :: overlap  ! function result
    REAL, INTENT(in) :: sigma_sq ! particle diameter squared

    INTEGER            :: i, j
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq, rij_mag
    REAL,    PARAMETER :: tol = 1.0e-4 

    overlap  = .FALSE.

    DO i = 1, n - 1
       DO j = i + 1, n

          rij(:) = r(:,i) - r(:,j)
          rij(:) = rij(:) - ANINT ( rij(:) )
          rij_sq = SUM ( rij**2 )

          IF ( rij_sq < sigma_sq ) THEN
             rij_mag = SQRT ( rij_sq / sigma_sq )
             WRITE ( unit=error_unit, fmt='(a,2i5,f15.8)' ) 'Warning: i,j,rij/sigma = ', i, j, rij_mag
             IF ( ( 1.0 - rij_mag ) > tol ) overlap = .TRUE.
          END IF

       END DO
    END DO

  END FUNCTION overlap

  SUBROUTINE collide ( i, j, sigma_sq, virial ) ! collision dynamics
    INTEGER, INTENT(in)  :: i, j
    REAL,    INTENT(in)  :: sigma_sq
    REAL,    INTENT(out) :: virial

    ! it is assumed that i and j are in contact
    ! the routine also computes the collisional virial

    REAL, DIMENSION(3) :: rij, vij
    REAL :: factor

    rij(:) = r(:,i) - r(:,j)
    rij(:) = rij(:) - ANINT ( rij(:) )
    vij(:) = v(:,i) - v(:,j)

    factor = DOT_PRODUCT ( rij, vij ) / sigma_sq
    vij = - factor * rij

    v(:,i) = v(:,i) + vij
    v(:,j) = v(:,j) - vij
    virial = DOT_PRODUCT ( vij, rij ) / 3.0
  END SUBROUTINE collide

END MODULE md_nve_hs_module
