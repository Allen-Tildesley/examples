! verlet_list_module.f90 (used by md_lj_vl_module.f90)
! Molecular dynamics simulation
MODULE verlet_list_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_list, finalize_list, make_list
  PUBLIC :: point, list

  INTEGER                              :: nl     ! size of list
  REAL                                 :: r_list ! list range parameter
  REAL                                 :: r_skin ! list skin parameter
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r_save ! saved positions for list (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: dr     ! displacements (3,n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: point  ! index to neighbour list (n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: list   ! Verlet neighbour list (nl)

CONTAINS

  SUBROUTINE initialize_list ( n, r_cut )
    INTEGER, INTENT(in) :: n
    REAL,    INTENT(in) :: r_cut

    REAL :: r_list_factor

    REAL, PARAMETER :: pi = 4.0*ATAN(1.0)
    NAMELIST /list_parameters/ r_list_factor

    ! Sensible default for r_list_factor
    r_list_factor = 1.2
    READ(*,nml=list_parameters)
    WRITE(*,'(''Verlet list factor = '',t40,f15.5)') r_list_factor
    IF ( r_list_factor <= 1.0 ) STOP 'r_list_factor must be > 1'
    r_list = r_cut * r_list_factor
    IF ( r_list > 0.5   ) STOP 'r_list too large'
    r_skin = r_list - r_cut
    WRITE(*,'(''Verlet list range (box units) = '',t40,f15.5)') r_list
    WRITE(*,'(''Verlet list skin  (box units) = '',t40,f15.5)') r_skin

    ! Estimate list size based on density + 10 per cent
    nl = CEILING ( 1.1*(4.0*pi/6.0)*(r_list**3)*REAL(n**2) )
    WRITE(*,'(''Verlet list size = '',t40,i15)') nl
    ALLOCATE ( r_save(3,n), dr(3,n), point(n), list(nl) )
  END SUBROUTINE initialize_list

  SUBROUTINE finalize_list
    DEALLOCATE ( r_save, dr, point, list )
  END SUBROUTINE finalize_list
  
  SUBROUTINE resize_list ! reallocates list array, somewhat larger
    
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp
    INTEGER                            :: nl_new

    nl_new = CEILING ( 1.25 * REAL(nl) )
    WRITE(*,'(''Warning: reallocating Verlet list array'')')

    ALLOCATE ( tmp(nl_new) ) ! new size for list
    tmp(1:nl) = list(:)      ! copy elements across

    CALL move_ALLOC ( tmp, list )
    nl = SIZE(list)
    WRITE(*,'(''New Verlet list size = '',t40,i15)') nl

  END SUBROUTINE resize_list

  SUBROUTINE make_list ( n, r )
    INTEGER,                 INTENT(in) :: n
    REAL,    DIMENSION(3,n), INTENT(in) :: r

    INTEGER            :: i, j, k
    REAL               :: r_list_sq, rij_sq, dr_sq_max
    REAL, DIMENSION(3) :: rij

    LOGICAL, SAVE :: first_call = .TRUE.

    IF ( .NOT. first_call ) THEN
       dr = r - r_save                         ! displacement since last list update
       dr = dr - ANINT ( dr )                  ! periodic boundaries in box=1
       dr_sq_max = MAXVAL ( SUM(dr**2,dim=1) ) ! squared maximum displacement
       IF ( 4.0*dr_sq_max < r_skin ** 2 )  RETURN ! no need to make list
    END IF

    first_call = .FALSE.

    k = 0
    r_list_sq = r_list ** 2

    DO i = 1, n - 1 ! Begin outer loop over atoms

       point(i) = k + 1

       DO j = i + 1, n ! Begin inner loop over partner atoms

          rij(:) = r(:,i) - r(:,j)
          rij(:) = rij(:) - ANINT ( rij(:) )
          rij_sq = SUM ( rij**2 )

          IF ( rij_sq < r_list_sq ) THEN

             k = k + 1
             IF ( k > nl ) CALL resize_list
             list(k) = j

          END IF

       END DO ! End inner loop over partner atoms

    END DO ! End outer loop over atoms

    point(n) = k + 1

    r_save = r

  END SUBROUTINE make_list

END MODULE verlet_list_module
