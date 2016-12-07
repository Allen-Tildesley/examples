! mc_hs_module.f90
! Overlap routines for MC simulation, hard spheres
MODULE mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: overlap_1, overlap, n_overlap

  ! Public data
  INTEGER,                             PUBLIC :: n ! Number of atoms
  REAL,   DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n)

  ! Private data
  INTEGER, PARAMETER :: lt = -1, gt = 1 ! j-range options

CONTAINS

  SUBROUTINE introduction

    WRITE ( unit=output_unit, fmt='(a)' ) 'Hard sphere potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'   
    WRITE ( unit=output_unit, fmt='(a)' ) 'Energy, kT = 1'   

  END SUBROUTINE introduction

  SUBROUTINE conclusion

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays
    ALLOCATE ( r(3,n) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r )
  END SUBROUTINE deallocate_arrays

  FUNCTION overlap ( box )
    LOGICAL            :: overlap ! Shows if an overlap was detected
    REAL,   INTENT(in) :: box     ! Simulation box length

    ! Actual calculation is performed by function overlap_1

    INTEGER :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in overlap'
    END IF

    DO i = 1, n - 1
       IF ( overlap_1 ( r(:,i), i, box, gt ) ) THEN
          overlap = .TRUE. ! Overlap detected
          RETURN           ! Return immediately
       END IF
    END DO

    overlap  = .FALSE. ! No overlaps detected

  END FUNCTION overlap

  FUNCTION overlap_1 ( ri, i, box, j_range ) RESULT ( overlap )
    LOGICAL                        :: overlap ! Shows if an overlap was detected
    REAL, DIMENSION(3), INTENT(in) :: ri      ! Coordinates of atom of interest
    INTEGER,            INTENT(in) :: i       ! Index of atom of interest
    REAL,               INTENT(in) :: box     ! Simulation box length
    INTEGER, OPTIONAL,  INTENT(in) :: j_range ! Optional partner index range

    ! Detects overlap of atom in ri
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i
    ! It is assumed that r is in units where box = 1

    INTEGER            :: j, j1, j2
    REAL               :: box_sq, rij_sq
    REAL, DIMENSION(3) :: rij

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in overlap_1'
    END IF

    IF ( PRESENT ( j_range ) ) THEN
       SELECT CASE ( j_range )
       CASE ( lt ) ! j < i
          j1 = 1
          j2 = i-1
       CASE ( gt ) ! j > i
          j1 = i+1
          j2 = n
       CASE default ! should never happen
          WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error ', j_range
          STOP 'Impossible error in overlap_1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    box_sq = box**2

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( i == j ) CYCLE ! Skip self

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )
       rij_sq = rij_sq * box_sq ! now in sigma=1 units

       IF ( rij_sq < 1.0 ) THEN
          overlap = .TRUE. ! Overlap detected
          RETURN           ! Return immediately
       END IF

    END DO ! End loop over selected range of partners

    overlap = .FALSE. ! No overlaps detected

  END FUNCTION overlap_1

  FUNCTION n_overlap ( box )
    INTEGER             :: n_overlap ! Counts overlaps
    REAL,    INTENT(in) :: box       ! Simulation box length

    ! This routine is used in the calculation of pressure
    ! Actual calculation is performed by function n_overlap_1

    INTEGER :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in n_overlap'
    END IF

    n_overlap  = 0

    DO i = 1, n - 1
       n_overlap = n_overlap + n_overlap_1 ( r(:,i), i, box, gt )
    END DO

  END FUNCTION n_overlap

  FUNCTION n_overlap_1 ( ri, i, box, j_range ) RESULT ( n_overlap )
    INTEGER                        :: n_overlap ! Counts overlaps
    REAL, DIMENSION(3), INTENT(in) :: ri        ! Coordinates of atom of interest
    INTEGER,            INTENT(in) :: i         ! Index of atom of interest
    REAL,               INTENT(in) :: box       ! Simulation box length
    INTEGER, OPTIONAL,  INTENT(in) :: j_range   ! Optional partner index range

    ! Counts overlaps of atom in ri
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i
    ! It is assumed that r is in units where box = 1
    ! This routine is used in the calculation of pressure

    INTEGER            :: j, j1, j2
    REAL               :: box_sq, rij_sq
    REAL, DIMENSION(3) :: rij

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in n_overlap_1'
    END IF

    IF ( PRESENT ( j_range ) ) THEN
       SELECT CASE ( j_range )
       CASE ( lt ) ! j < i
          j1 = 1
          j2 = i-1
       CASE ( gt ) ! j > i
          j1 = i+1
          j2 = n
       CASE default ! should never happen
          WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error ', j_range
          STOP 'Impossible error in n_overlap_1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    box_sq    = box**2
    n_overlap = 0

    DO j = j1, j2

       IF ( i == j ) CYCLE ! Skip self

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )
       rij_sq = rij_sq * box_sq ! Now in sigma=1 units

       IF ( rij_sq < 1.0 ) n_overlap = n_overlap + 1

    END DO

  END FUNCTION n_overlap_1

END MODULE mc_module
