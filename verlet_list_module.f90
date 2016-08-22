! verlet_list_module.f90
! Verlet list handling routines for MD simulation
MODULE verlet_list_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: initialize_list, finalize_list, make_list
  PUBLIC :: point, list

  ! The initialize_list routine reads the value of r_list_factor
  ! from standard input using a NAMELIST nml_list
  ! Leave namelist empty to accept supplied default
  ! It is assumed that all positions and displacements are divided by box
  ! r_list_box is set to r_cut_box*r_list_factor
  
  INTEGER                                         :: nl         ! size of list
  REAL                                            :: r_list_box ! list range parameter / box length
  REAL                                            :: r_skin_box ! list skin parameter / box_length
  REAL,    DIMENSION(:,:), ALLOCATABLE            :: r_save     ! saved positions for list (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE            :: dr         ! displacements (3,n)
  INTEGER, DIMENSION(:),   ALLOCATABLE, PROTECTED :: point      ! index to neighbour list (n)
  INTEGER, DIMENSION(:),   ALLOCATABLE, PROTECTED :: list       ! Verlet neighbour list (nl)

CONTAINS

  SUBROUTINE initialize_list ( n, r_cut_box )
    INTEGER, INTENT(in) :: n         ! number of particles
    REAL,    INTENT(in) :: r_cut_box ! r_cut / box

    REAL    :: r_list_factor
    INTEGER :: ioerr

    REAL, PARAMETER :: pi = 4.0*ATAN(1.0)
    NAMELIST /nml_list/ r_list_factor

        WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Verlet list based on r_cut/box =', r_cut_box

    ! Sensible default for r_list_factor
    r_list_factor = 1.2
    READ ( unit=input_unit, nml=nml_list, iostat=ioerr ) ! namelist input
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error reading namelist nml_list from standard input', ioerr
       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
       STOP 'Error in initialize_list'
    END IF
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Verlet list factor = ', r_list_factor
    IF ( r_list_factor <= 1.0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_list_factor must be > 1', r_list_factor
       STOP 'Error in initialize_list'
    END IF
    r_list_box = r_cut_box * r_list_factor
    IF ( r_list_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_list/box too large', r_list_box
       STOP 'Error in initialize_list'
    END IF
    r_skin_box = r_list_box - r_cut_box
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Verlet list range (box units) = ', r_list_box
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Verlet list skin  (box units) = ', r_skin_box

    ! Estimate list size based on density + 10 per cent
    nl = CEILING ( 1.1*(4.0*pi/3.0)*(r_list_box**3)*REAL(n*(n-1)) / 2.0 )
    WRITE ( unit=output_unit, fmt='(a,t40,i15)') 'Verlet list size = ', nl
    ALLOCATE ( r_save(3,n), dr(3,n), point(n), list(nl) )
  END SUBROUTINE initialize_list

  SUBROUTINE finalize_list
    DEALLOCATE ( r_save, dr, point, list )
  END SUBROUTINE finalize_list
  
  SUBROUTINE resize_list ! reallocates list array, somewhat larger
    
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp
    INTEGER                            :: nl_new

    nl_new = CEILING ( 1.25 * REAL(nl) )
    WRITE( unit=output_unit, fmt='(a)', advance='no' ) 'Warning: new Verlet list array size = '

    ALLOCATE ( tmp(nl_new) ) ! new size for list
    tmp(1:nl) = list(:)      ! copy elements across

    CALL MOVE_ALLOC ( tmp, list )
    nl = SIZE(list)
    WRITE( unit=error_unit, fmt='(t60,i15)') nl

  END SUBROUTINE resize_list

  SUBROUTINE make_list ( n, r )
    INTEGER,                 INTENT(in) :: n
    REAL,    DIMENSION(3,n), INTENT(in) :: r

    INTEGER            :: i, j, k
    REAL               :: r_list_box_sq, rij_sq, dr_sq_max
    REAL, DIMENSION(3) :: rij

    LOGICAL, SAVE :: first_call = .TRUE.

    IF ( .NOT. first_call ) THEN
       dr = r - r_save                             ! displacement since last list update
       dr = dr - ANINT ( dr )                      ! periodic boundaries in box=1 units
       dr_sq_max = MAXVAL ( SUM(dr**2,dim=1) )     ! squared maximum displacement
       IF ( 4.0*dr_sq_max < r_skin_box**2 ) RETURN ! no need to make list
    END IF

    first_call = .FALSE.

    k = 0
    r_list_box_sq = r_list_box ** 2

    DO i = 1, n - 1 ! Begin outer loop over atoms

       point(i) = k + 1

       DO j = i + 1, n ! Begin inner loop over partner atoms

          rij(:) = r(:,i) - r(:,j)
          rij(:) = rij(:) - ANINT ( rij(:) )
          rij_sq = SUM ( rij**2 )

          IF ( rij_sq < r_list_box_sq ) THEN

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
