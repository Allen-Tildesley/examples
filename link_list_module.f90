! link_list_module.f90 (used by md_lj_ll_module.f90)
! Monte Carlo or molecular dynamics simulation
MODULE link_list_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_list, finalize_list, make_list, check_list
  PUBLIC :: move_in_list, create_in_list, destroy_in_list, c_index
  PUBLIC :: nc, head, list, c

  INTEGER,                                PROTECTED :: nc     ! dimensions of head array, assume cubic box
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE, PROTECTED :: head   ! head(0:nc-1,0:nc-1,0:nc-1), assume cubic box
  INTEGER, DIMENSION(:),     ALLOCATABLE, PROTECTED :: list   ! list(n)
  INTEGER, DIMENSION(:,:),   ALLOCATABLE, PROTECTED :: c      ! c(3,n) 3D cell index of each atom

CONTAINS

  SUBROUTINE initialize_list ( n, r_cut )
    INTEGER, INTENT(in) :: n
    REAL,    INTENT(in) :: r_cut ! We assume that r_cut (in box=1 units) never changes

    WRITE(*,'(''Link list initialization'')')
    WRITE(*,'(''Link cells based on r_cut ='',t40,f15.5)') r_cut
    nc = FLOOR ( 1.0 / r_cut ) ! number of cells
    IF ( nc < 3 ) STOP 'System is too small to use link cells'

    ALLOCATE ( list(n), c(3,n) )
    ALLOCATE ( head(0:nc-1,0:nc-1,0:nc-1) )

  END SUBROUTINE initialize_list

  SUBROUTINE finalize_list
    DEALLOCATE ( list, c )
    DEALLOCATE ( head )
  END SUBROUTINE finalize_list

  SUBROUTINE make_list ( n, r )
    INTEGER,                 INTENT(in) :: n ! number of atoms
    REAL,    DIMENSION(3,n), INTENT(in) :: r ! atom coordinates

    ! Local variables
    INTEGER :: i

    head(:,:,:) = 0

    ! Allocate to cells
    DO i = 1, n 
       c(:,i)  = c_index ( r(:,i) )    ! Index function
       CALL create_in_list ( i, c(:,i) )
    END DO

  END SUBROUTINE make_list

  FUNCTION c_index ( ri ) RESULT ( ci )
    INTEGER, DIMENSION(3)             :: ci ! 3D cell index in range 0..nc-1
    REAL,    DIMENSION(3), INTENT(in) :: ri ! position in box = 1 units

    ! We do not want to do any periodic imaging here, so as to cope with 
    ! Lees-Edwards boundaries as well as normal ones
    IF ( ANY ( ABS(ri) > 0.5 ) ) STOP 'atom not in main box in c_index'
    ci(:) = FLOOR ( ( ri(:) + 0.5 ) * REAL(nc) ) ! The index formula

    ! Guard against small chance of roundoff error
    WHERE ( ci(:) < 0    ) ci(:) = 0
    WHERE ( ci(:) > nc-1 ) ci(:) = nc-1

  END FUNCTION c_index
  
  SUBROUTINE move_in_list ( i, ci ) ! Move i from current cell to ci
    INTEGER,               INTENT(in) :: i
    INTEGER, DIMENSION(3), INTENT(in) :: ci

    IF ( ALL ( ci(:) == c(:,i) ) ) RETURN ! no need to do anything

    CALL destroy_in_list ( i, c(:,i) ) ! Remove from old cell
    CALL create_in_list  ( i, ci(:)  ) ! Add to new cell

  END SUBROUTINE move_in_list

  SUBROUTINE create_in_list ( i, ci ) ! Create i in cell ci
    INTEGER,               INTENT(in) :: i
    INTEGER, DIMENSION(3), INTENT(in) :: ci

    list(i)                 = head(ci(1),ci(2),ci(3))
    head(ci(1),ci(2),ci(3)) = i 
    c(:,i)                  = ci(:)

  END SUBROUTINE create_in_list

  SUBROUTINE destroy_in_list ( i, ci )      ! Destroy i in cell ci
    INTEGER,               INTENT(in) :: i
    INTEGER, DIMENSION(3), INTENT(in) :: ci

    INTEGER :: this, next

    this = head(ci(1),ci(2),ci(3))

    IF ( this == i ) THEN ! i is the head atom in that cell

       head(ci(1),ci(2),ci(3)) = list(i) ! simply point head at next atom

    ELSE

       DO ! traverse link-list, keeping previous index
          next = list(this) ! find the next one

          IF ( next == i ) THEN ! found our atom, just link over it
             list(this) = list(i)
             EXIT ! leave the loop

          ELSE IF ( next == 0 ) THEN ! this should never happen
             STOP 'could not find particle in its cell'

          ELSE ! move on to the next
             this = next
          END IF

       END DO

    END IF

  END SUBROUTINE destroy_in_list

  SUBROUTINE check_list ( n, r ) ! Routine to check consistency of cell lists
    INTEGER,                 INTENT(in) :: n
    REAL,    DIMENSION(3,n), INTENT(in) :: r

    INTEGER, DIMENSION(3) :: ci
    INTEGER               :: ci1, ci2, ci3, i

    ! First check that the cells are correct
    DO i = 1, n 
       ci = c_index ( r(:,i) ) ! Index function
       IF ( ANY ( ci(:) /= c(:,i) ) ) STOP 'check_list index error 1'
    END DO

    ! Next check that atoms in cells really are in those cells
    DO ci1 = 0, nc-1
       DO ci2 = 0, nc-1
          DO ci3 = 0, nc-1
             ci = [ ci1, ci2, ci3 ]
             i = head(ci1,ci2,ci3)
             DO
                IF ( i == 0 ) EXIT
                IF ( ANY ( ci(:) /= c(:,i) ) ) STOP 'check_list index error 2'
                i = list(i)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE check_list

END MODULE link_list_module
