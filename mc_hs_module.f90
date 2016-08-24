! mc_hs_module.f90
! Overlap routines for MC simulation, hard spheres
MODULE mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, lt, ne, gt
  PUBLIC :: model_description, allocate_arrays, deallocate_arrays
  PUBLIC :: overlap_1, overlap, n_overlap_1, n_overlap

  INTEGER                             :: n ! number of atoms
  REAL,   DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,n)

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options
  REAL,    PARAMETER :: sigma = 1.0             ! hard-sphere diameter (unit of length)

CONTAINS

  SUBROUTINE model_description ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Hard sphere potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ', sigma    
  END SUBROUTINE model_description
  
  SUBROUTINE allocate_arrays
    ALLOCATE ( r(3,n) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r )
  END SUBROUTINE deallocate_arrays

  FUNCTION overlap ( box )
    LOGICAL             :: overlap ! shows if an overlap was detected
    REAL,    INTENT(in) :: box     ! simulation box length

    INTEGER :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in overlap'
    END IF

    overlap  = .FALSE.

    DO i = 1, n - 1
       IF ( overlap_1 ( r(:,i), i, gt, box ) ) THEN
          overlap = .TRUE.
          EXIT ! jump out of loop
       END IF
    END DO

  END FUNCTION overlap

  FUNCTION overlap_1 ( ri, i, j_range, box ) RESULT ( overlap )
    LOGICAL                        :: overlap    ! shows if an overlap was detected
    REAL, DIMENSION(3), INTENT(in) :: ri         ! coordinates of atom of interest
    INTEGER,            INTENT(in) :: i, j_range ! index, and partner index range
    REAL,               INTENT(in) :: box        ! simulation box length

    ! Detects overlap of atom in ri
    ! with j/=i, j>i, or j<i depending on j_range
    ! It is assumed that r is in units where box = 1

    INTEGER            :: j, j1, j2
    REAL               :: box_sq, rij_sq
    REAL, DIMENSION(3) :: rij

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in overlap_1'
    END IF

    box_sq = box**2

    overlap = .FALSE.

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

    DO j = j1, j2

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )
       rij_sq = rij_sq * box_sq ! now in sigma=1 units

       IF ( rij_sq < 1.0 ) THEN
          overlap = .TRUE.
          EXIT ! jump out of loop
       END IF

    END DO

  END FUNCTION overlap_1

  FUNCTION n_overlap ( box )
    INTEGER             :: n_overlap ! counts overlaps
    REAL,    INTENT(in) :: box       ! simulation box length

    INTEGER :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in n_overlap'
    END IF

    n_overlap  = 0

    DO i = 1, n - 1
       n_overlap = n_overlap + n_overlap_1 ( r(:,i), i, gt, box )
    END DO

  END FUNCTION n_overlap

  FUNCTION n_overlap_1 ( ri, i, j_range, box ) RESULT ( n_overlap )
    INTEGER                        :: n_overlap  ! counts overlaps
    REAL, DIMENSION(3), INTENT(in) :: ri         ! coordinates of atom of interest
    INTEGER,            INTENT(in) :: i, j_range ! index, and partner index range
    REAL,               INTENT(in) :: box        ! simulation box length

    ! Counts overlaps of atom in ri
    ! with j/=i, j>i, or j<i depending on j_range
    ! It is assumed that r is in units where box = 1

    INTEGER            :: j, j1, j2
    REAL               :: box_sq, rij_sq
    REAL, DIMENSION(3) :: rij

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in n_overlap_1'
    END IF

    box_sq = box**2

    n_overlap = 0

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

    DO j = j1, j2

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )
       rij_sq = rij_sq * box_sq ! now in sigma=1 units

       IF ( rij_sq < 1.0 ) n_overlap = n_overlap + 1

    END DO

  END FUNCTION n_overlap_1

END MODULE mc_module
