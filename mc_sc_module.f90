! mc_sc_module.f90
! Overlap routines for MC simulation, hard spherocylinders
MODULE mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, r, e
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: overlap_1, overlap, n_overlap

  INTEGER                             :: n ! number of atoms
  REAL,   DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,n)
  REAL,   DIMENSION(:,:), ALLOCATABLE :: e ! orientations (3,n)

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range

CONTAINS

  SUBROUTINE introduction ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)' ) 'Hard spherocylinder potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'  
    WRITE ( unit=output_unit, fmt='(a)' ) 'Energy, kT = 1'   
  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output
    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'
  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays
    ALLOCATE ( r(3,n), e(3,n) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, e )
  END SUBROUTINE deallocate_arrays
  
  FUNCTION overlap ( box, length )
    LOGICAL             :: overlap ! Shows if an overlap was detected
    REAL,    INTENT(in) :: box     ! Simulation box length
    REAL,    INTENT(in) :: length  ! Cylinder length

    ! Actual calculation is performed by function overlap_1

    INTEGER :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in overlap'
    END IF
    IF ( n > SIZE(e,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for e', n, SIZE(e,dim=2)
       STOP 'Error in overlap'
    END IF

    DO i = 1, n - 1
       IF ( overlap_1 ( r(:,i), e(:,i), i, box, length, gt ) ) THEN
          overlap = .TRUE. ! Overlap detected
          return           ! Return immediately
       END IF
    END DO

    overlap  = .FALSE. ! No overlaps detected

  END FUNCTION overlap

  FUNCTION overlap_1 ( ri, ei, i, box, length, j_range ) RESULT ( overlap )
    LOGICAL                        :: overlap ! Shows if an overlap was detected
    REAL, DIMENSION(3), INTENT(in) :: ri      ! Coordinates of molecule of interest
    REAL, DIMENSION(3), INTENT(in) :: ei      ! Orientation (vector) of molecule of interest
    INTEGER,            INTENT(in) :: i       ! Index of molecule of interest
    REAL,               INTENT(in) :: box     ! Simulation box length
    REAL,               INTENT(in) :: length  ! Cylinder length
    INTEGER, OPTIONAL,  INTENT(in) :: j_range ! Optional partner index range

    ! Detects overlap of molecule in ri/ei
    ! The coordinates in ri and ei are not necessarily identical with those in r(:,i) and e(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i
    ! It is assumed that r is in units where box = 1

    INTEGER            :: j, j1, j2
    REAL               :: box_sq, range, range_box_sq, rij_sq, rei, rej, eij, sij_sq
    REAL, DIMENSION(3) :: rij

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in overlap_1'
    END IF
    IF ( n > SIZE(e,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for e', n, SIZE(e,dim=2)
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

    box_sq       = box**2
    range        = 1.0 + length         ! centre-centre interaction range
    range_box_sq = ( range / box ) ** 2 ! squared range in box=1 units

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )
       IF ( rij_sq > range_box_sq ) CYCLE ! No possibility of overlap

       rij_sq = rij_sq * box_sq ! now in sigma=1 units
       rij    = rij    * box    ! now in sigma=1 units
       rei = DOT_PRODUCT ( rij, ei     )
       rej = DOT_PRODUCT ( rij, e(:,j) )
       eij = DOT_PRODUCT ( ei,  e(:,j) )

       sij_sq = sc_dist_sq ( rij_sq, rei, rej, eij, length )
       IF ( sij_sq < 1.0 ) THEN
          overlap = .TRUE. ! Overlap detected
          return           ! Return immediately
       END IF

    END DO ! End loop over selected range of partners

    overlap = .FALSE. ! No overlaps detected

  END FUNCTION overlap_1

  FUNCTION n_overlap ( box, length )
    INTEGER             :: n_overlap ! Returns number of overlaps
    REAL,    INTENT(in) :: box       ! simulation box length
    REAL,    INTENT(in) :: length    ! cylinder length

    ! This routine is used in the calculation of pressure
    ! Actual calculation is performed by function n_overlap_1

    INTEGER :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in n_overlap'
    END IF
    IF ( n > SIZE(e,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for e', n, SIZE(e,dim=2)
       STOP 'Error in n_overlap'
    END IF

    n_overlap  = 0

    DO i = 1, n - 1
       n_overlap = n_overlap + n_overlap_1 ( r(:,i), e(:,i), i, box, length, gt )
    END DO

  END FUNCTION n_overlap

  FUNCTION n_overlap_1 ( ri, ei, i, box, length, j_range ) RESULT ( n_overlap )
    INTEGER                        :: n_overlap ! Returns number of overlaps
    REAL, DIMENSION(3), INTENT(in) :: ri        ! Coordinates of molecule of interest
    REAL, DIMENSION(3), INTENT(in) :: ei        ! Orientation (vector) of molecule of interest
    INTEGER,            INTENT(in) :: i         ! Index of molecule of interest
    REAL,               INTENT(in) :: box       ! Simulation box length
    REAL,               INTENT(in) :: length    ! Cylinder length
    INTEGER, OPTIONAL,  INTENT(in) :: j_range   ! Optional partner index range

    ! Counts overlaps of molecule in ri/ei
    ! The coordinates in ri and ei are not necessarily identical with those in r(:,i) and e(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i
    ! It is assumed that r is in units where box = 1
    ! This routine is used in the calculation of pressure

    INTEGER            :: j, j1, j2
    REAL               :: box_sq, range, range_box_sq, rij_sq, rei, rej, eij, sij_sq
    REAL, DIMENSION(3) :: rij

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in n_overlap_1'
    END IF
    IF ( n > SIZE(e,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for e', n, SIZE(e,dim=2)
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

    box_sq       = box**2
    range        = 1.0 + length         ! centre-centre interaction range
    range_box_sq = ( range / box ) ** 2 ! squared range in box=1 units

    n_overlap = 0

    DO j = j1, j2

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )
       IF ( rij_sq > range_box_sq ) CYCLE

       rij_sq = rij_sq * box_sq ! now in sigma=1 units
       rij    = rij    * box    ! now in sigma=1 units
       rei = DOT_PRODUCT ( rij, ei     )
       rej = DOT_PRODUCT ( rij, e(:,j) )
       eij = DOT_PRODUCT ( ei,  e(:,j) )

       sij_sq = sc_dist_sq ( rij_sq, rei, rej, eij, length )

       IF ( sij_sq < 1.0 ) n_overlap = n_overlap + 1

    END DO

  END FUNCTION n_overlap_1

  FUNCTION sc_dist_sq ( rij_sq, rei, rej, eij, length ) RESULT ( sij_sq )
    REAL             :: sij_sq                        ! squared distance between line segments
    REAL, INTENT(in) :: rij_sq, rei, rej, eij, length ! geometric parameters

    REAL            :: sin_sq, ci, cj, ai, aj, di, dj, length2
    REAL, PARAMETER :: tol = 0.000001

    sin_sq  = 1.0 - eij**2
    length2 = length

    IF ( sin_sq < tol ) THEN
       ci = -rei
       cj =  rej
    ELSE
       ci = ( - rei + eij * rej ) / sin_sq
       cj = (   rej - eij * rei ) / sin_sq
    ENDIF
    ai = ABS ( ci )
    aj = ABS ( cj )
    IF ( ai > length2 ) ci = SIGN ( length2, ci )
    IF ( aj > length2 ) cj = SIGN ( length2, cj )
    IF ( ai > aj ) THEN
       cj =  rej + ci * eij
    ELSE
       ci = -rei + cj * eij
    ENDIF
    ai = ABS ( ci )
    aj = ABS ( cj )
    IF ( ai > length2 ) ci = SIGN ( length2, ci )
    IF ( aj > length2 ) cj = SIGN ( length2, cj )
    di =  2.0 * rei + ci - cj * eij
    dj = -2.0 * rej + cj - ci * eij
    sij_sq = rij_sq + ci * di + cj * dj
  END FUNCTION sc_dist_sq

END MODULE mc_module
