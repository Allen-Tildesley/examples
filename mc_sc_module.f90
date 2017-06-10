! mc_sc_module.f90
! Overlap routines for MC simulation, hard spherocylinders
MODULE mc_module

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
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: overlap_1, overlap, n_overlap

  ! Public data
  INTEGER,                             PUBLIC :: n ! Number of atoms
  REAL,   DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n)
  REAL,   DIMENSION(:,:), ALLOCATABLE, PUBLIC :: e ! Orientations (3,n)

  ! Private data
  REAL,    PARAMETER :: pi     = 4.0*atan(1.0)
  REAL,    PARAMETER :: length = 5.0                            ! Cylinder length L (in units where D=1)
  REAL,    PARAMETER :: vmol   = pi * ( 0.25*length + 1.0/6.0 ) ! Spherocylinder volume

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range

CONTAINS

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)'       ) 'Hard spherocylinder potential'
    WRITE ( unit=output_unit, fmt='(a,f15.6)' ) 'Spherocylinder L/D ratio',    length
    WRITE ( unit=output_unit, fmt='(a,f15.6)' ) 'Spherocylinder volume/D**3',  vmol
    WRITE ( unit=output_unit, fmt='(a)'       ) 'Diameter, D = 1'  
    WRITE ( unit=output_unit, fmt='(a)'       ) 'Energy, kT = 1'   

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays
    IMPLICIT NONE

    ALLOCATE ( r(3,n), e(3,n) )

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, e )

  END SUBROUTINE deallocate_arrays

  FUNCTION overlap ( box )
    IMPLICIT NONE
    LOGICAL             :: overlap ! Shows if an overlap was detected
    REAL,    INTENT(in) :: box     ! Simulation box length

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
       IF ( overlap_1 ( r(:,i), e(:,i), i, box, gt ) ) THEN
          overlap = .TRUE. ! Overlap detected
          RETURN           ! Return immediately
       END IF
    END DO

    overlap  = .FALSE. ! No overlaps detected

  END FUNCTION overlap

  FUNCTION overlap_1 ( ri, ei, i, box, j_range ) RESULT ( overlap )
    IMPLICIT NONE
    LOGICAL                        :: overlap ! Shows if an overlap was detected
    REAL, DIMENSION(3), INTENT(in) :: ri      ! Coordinates of molecule of interest
    REAL, DIMENSION(3), INTENT(in) :: ei      ! Orientation (vector) of molecule of interest
    INTEGER,            INTENT(in) :: i       ! Index of molecule of interest
    REAL,               INTENT(in) :: box     ! Simulation box length
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

    range = 1.0 + length ! Centre-centre interaction range (D=1 units)

    IF ( range > 0.5*box ) THEN ! Check for box being too small
       WRITE ( unit=error_unit, fmt='(a,2f15.6)' ) 'Box too small', box, range
       STOP 'Error in overlap_1'
    END IF

    range_box_sq = ( range / box ) ** 2 ! Squared range in box=1 units
    box_sq       = box**2               ! Squared box length

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( i == j ) CYCLE ! Skip self

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )
       IF ( rij_sq > range_box_sq ) CYCLE ! No possibility of overlap

       rij_sq = rij_sq * box_sq ! Now in D=1 units
       rij    = rij    * box    ! Now in D=1 units
       rei = DOT_PRODUCT ( rij, ei     )
       rej = DOT_PRODUCT ( rij, e(:,j) )
       eij = DOT_PRODUCT ( ei,  e(:,j) )

       sij_sq = dist_sq ( rij_sq, rei, rej, eij, length ) ! Squared distance between line segments
       IF ( sij_sq < 1.0 ) THEN
          overlap = .TRUE. ! Overlap detected
          RETURN           ! Return immediately
       END IF

    END DO ! End loop over selected range of partners

    overlap = .FALSE. ! No overlaps detected

  END FUNCTION overlap_1

  FUNCTION n_overlap ( box )
    IMPLICIT NONE
    INTEGER             :: n_overlap ! Returns number of overlaps
    REAL,    INTENT(in) :: box       ! Simulation box length

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
       n_overlap = n_overlap + n_overlap_1 ( r(:,i), e(:,i), i, box, gt )
    END DO

  END FUNCTION n_overlap

  FUNCTION n_overlap_1 ( ri, ei, i, box, j_range ) RESULT ( n_overlap )
    IMPLICIT NONE
    INTEGER                        :: n_overlap ! Returns number of overlaps
    REAL, DIMENSION(3), INTENT(in) :: ri        ! Coordinates of molecule of interest
    REAL, DIMENSION(3), INTENT(in) :: ei        ! Orientation (vector) of molecule of interest
    INTEGER,            INTENT(in) :: i         ! Index of molecule of interest
    REAL,               INTENT(in) :: box       ! Simulation box length
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

    range = 1.0 + length ! Centre-centre interaction range (D=1 units)

    IF ( range > 0.5*box ) THEN ! Check for box being too small
       WRITE ( unit=error_unit, fmt='(a,2f15.6)' ) 'Box too small', box, range
       STOP 'Error in n_overlap_1'
    END IF

    range_box_sq = ( range / box ) ** 2 ! Squared range in box=1 units
    box_sq       = box**2               ! Squared box length

    n_overlap = 0

    DO j = j1, j2

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )
       IF ( rij_sq > range_box_sq ) CYCLE ! No possibility of overlap

       rij_sq = rij_sq * box_sq ! Now in D=1 units
       rij    = rij    * box    ! Now in D=1 units
       rei = DOT_PRODUCT ( rij, ei     )
       rej = DOT_PRODUCT ( rij, e(:,j) )
       eij = DOT_PRODUCT ( ei,  e(:,j) )

       sij_sq = dist_sq ( rij_sq, rei, rej, eij, length ) ! Squared distance between line segments

       IF ( sij_sq < 1.0 ) n_overlap = n_overlap + 1

    END DO

  END FUNCTION n_overlap_1

  FUNCTION dist_sq ( rij_sq, rei, rej, eij, ell ) RESULT ( sij_sq )
    IMPLICIT NONE
    REAL             :: sij_sq ! Returns squared distance between line segments
    REAL, INTENT(in) :: rij_sq ! Squared centre-centre distance
    REAL, INTENT(in) :: rei    ! Scalar product rij.ei where ei is unit vector along i
    REAL, INTENT(in) :: rej    ! Scalar product rij.ej where ej is unit vector along j
    REAL, INTENT(in) :: eij    ! Scalar product ei.ej
    REAL, INTENT(in) :: ell    ! Line segment length

    REAL            :: sin_sq, ci, cj, ai, aj, di, dj, ell2
    REAL, PARAMETER :: tol = 1.e-6

    sin_sq = 1.0 - eij**2 ! Squared sine of angle between line segments
    ell2   = ell / 2.0    ! Half the line segment length

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

    sij_sq = rij_sq + ci * di + cj * dj

  END FUNCTION dist_sq

END MODULE mc_module
