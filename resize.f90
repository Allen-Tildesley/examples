
  SUBROUTINE resize 
    IMPLICIT NONE

    ! Reallocates r array, twice as large
    ! This is employed by mc_zvt_lj, grand canonical ensemble

    REAL, DIMENSION(:,:), ALLOCATABLE :: tmp
    INTEGER                           :: n_old, n_new

    n_old = SIZE(r,dim=2)
    n_new = 2*n_old
    WRITE( unit=output_unit, fmt='(a,i10,a,i10)' ) 'Reallocating r from old ', n_old, ' to ', n_new

    ALLOCATE ( tmp(3,n_new) ) ! New size for r
    tmp(:,1:n_old) = r(:,:)   ! Copy elements across

    CALL MOVE_ALLOC ( tmp, r )

  END SUBROUTINE resize
