! mc_sc_module.f90
! Overlap routines for MC simulation, hard spherocylinders
MODULE mc_sc_module

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, e, lt, ne, gt
  PUBLIC :: initialize, finalize, overlap_1, overlap, n_overlap_1, n_overlap

  INTEGER                             :: n ! number of atoms
  REAL,   DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,:)
  REAL,   DIMENSION(:,:), ALLOCATABLE :: e ! orientations (3,:)

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options

CONTAINS

  subroutine initialize
    ALLOCATE ( r(3,n), e(3,n) )
  end subroutine initialize

  subroutine finalize
    DEALLOCATE ( r, e )
  end subroutine finalize
  
  FUNCTION overlap ( sigma, length )
    LOGICAL             :: overlap ! shows if an overlap was detected
    REAL,    INTENT(in) :: sigma   ! cylinder diameter
    REAL,    INTENT(in) :: length  ! cylinder length

    INTEGER :: i

    IF ( n > SIZE(r,dim=2) ) STOP 'Array bounds error for r in overlap'
    IF ( n > SIZE(e,dim=2) ) STOP 'Array bounds error for e in overlap'

    overlap  = .FALSE.

    DO i = 1, n - 1
       IF ( overlap_1 ( r(:,i), e(:,i), i, gt, sigma, length ) ) THEN
          overlap = .TRUE.
          EXIT ! jump out of loop
       END IF
    END DO

  END FUNCTION overlap

  FUNCTION overlap_1 ( ri, ei, i, j_range, sigma, length ) RESULT ( overlap )
    LOGICAL                        :: overlap    ! shows if an overlap was detected
    REAL, DIMENSION(3), INTENT(in) :: ri         ! coordinates of molecule of interest
    REAL, DIMENSION(3), INTENT(in) :: ei         ! orientation of molecule of interest
    INTEGER,            INTENT(in) :: i, j_range ! index, and partner index range
    REAL,               INTENT(in) :: sigma      ! cylinder diameter
    REAL,               INTENT(in) :: length     ! cylinder length

    ! Detects overlap of atom in ri
    ! with j/=i, j>i, or j<i depending on j_range
    ! It is assumed that r, sigma are in units where box = 1

    INTEGER            :: j, j1, j2
    REAL               :: sigma_sq, range_sq, rij_sq, rei, rej, eij
    REAL, DIMENSION(3) :: rij

    IF ( n > SIZE(r,dim=2) ) STOP 'Array bounds error for r in energy'
    IF ( n > SIZE(e,dim=2) ) STOP 'Array bounds error for e in energy'

    sigma_sq = sigma**2
    range_sq = ( sigma + length ) ** 2

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
       IF ( rij_sq > range_sq ) CYCLE

       rei = DOT_PRODUCT ( rij, ei     )
       rej = DOT_PRODUCT ( rij, e(:,j) )
       eij = DOT_PRODUCT ( ei,  e(:,j) )

       rij_sq = sc_dist_sq ( rij_sq, rei, rej, eij, length )
       IF ( rij_sq < sigma_sq ) THEN
          overlap = .TRUE.
          EXIT ! jump out of loop
       END IF

    END DO

  END FUNCTION overlap_1

  FUNCTION n_overlap ( sigma, length )
    INTEGER             :: n_overlap ! counts overlaps
    REAL,    INTENT(in) :: sigma   ! cylinder diameter
    REAL,    INTENT(in) :: length  ! cylinder length

    INTEGER :: i

    IF ( n > SIZE(r,dim=2) ) STOP 'Array bounds error for r in energy'

    n_overlap  = 0

    DO i = 1, n - 1
       n_overlap = n_overlap + n_overlap_1 ( r(:,i), e(:,i), i, gt, sigma, length )
    END DO

  END FUNCTION n_overlap

  FUNCTION n_overlap_1 ( ri, ei, i, j_range, sigma, length ) RESULT ( n_overlap )
    INTEGER                        :: n_overlap  ! counts overlaps
    REAL, DIMENSION(3), INTENT(in) :: ri         ! coordinates of molecule of interest
    REAL, DIMENSION(3), INTENT(in) :: ei         ! orientation of molecule of interest
    INTEGER,            INTENT(in) :: i, j_range ! index, and partner index range
    REAL,               INTENT(in) :: sigma      ! cylinder diameter
    REAL,               INTENT(in) :: length     ! cylinder length

    ! Counts overlaps of atom in ri
    ! with j/=i, j>i, or j<i depending on j_range
    ! It is assumed that r, sigma are in units where box = 1

    INTEGER            :: j, j1, j2
    REAL               :: sigma_sq, range_sq, rij_sq, rei, rej, eij
    REAL, DIMENSION(3) :: rij

    IF ( n > SIZE(r,dim=2) ) STOP 'Array bounds error for r in energy'

    sigma_sq = sigma**2
    range_sq = ( sigma + length ) ** 2

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
       IF ( rij_sq > range_sq ) CYCLE

       rei = DOT_PRODUCT ( rij, ei     )
       rej = DOT_PRODUCT ( rij, e(:,j) )
       eij = DOT_PRODUCT ( ei,  e(:,j) )

       rij_sq = sc_dist_sq ( rij_sq, rei, rej, eij, length )

       IF ( rij_sq < sigma_sq ) n_overlap = n_overlap + 1

    END DO

  END FUNCTION n_overlap_1

  FUNCTION sc_dist_sq ( rij_sq, rei, rej, eij, length ) RESULT ( sij_sq )
    REAL             :: sij_sq                        ! squared distance between line segments
    REAL, INTENT(in) :: rij_sq, rei, rej, eij, length ! geometric parameters

    REAL :: sin_sq, ci, cj, ai, aj, di, dj, length2
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

END MODULE mc_sc_module
