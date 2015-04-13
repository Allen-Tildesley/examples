! link_list_module.f90 (used by md_lj_ll_module.f90)
! Monte Carlo or molecular dynamics simulation
MODULE link_list_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_list, finalize_list, make_list, update_list, check_list, get_neighbours
  PUBLIC :: nj, j_list, nc, head, list, c, dc, nk

  ! Set up vectors to each cell in neighbourhood of 3x3x3 cells in cubic lattice
  INTEGER, PARAMETER :: nk = 13 
  INTEGER, DIMENSION(3,-nk:nk), PARAMETER :: dc = RESHAPE( [ &
       &   -1,-1,-1,    0,-1,-1,    1,-1,-1, &
       &   -1, 1,-1,    0, 1,-1,    1, 1,-1, &
       &   -1, 0,-1,    1, 0,-1,    0, 0,-1, &
       &    0,-1, 0,    1,-1, 0,   -1,-1, 0, &
       &   -1, 0, 0,    0, 0, 0,    1, 0, 0, &
       &    1, 1, 0,   -1, 1, 0,    0, 1, 0, &
       &    0, 0, 1,   -1, 0, 1,    1, 0, 1, &
       &   -1,-1, 1,    0,-1, 1,    1,-1, 1, &
       &   -1, 1, 1,    0, 1, 1,    1, 1, 1    ], [ 3, 2*nk+1 ] )
  INTEGER,                                PROTECTED :: nc     ! dimensions of head array, assume cubic box
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE, PROTECTED :: head   ! head(0:nc-1,0:nc-1,0:nc-1), assume cubic box
  INTEGER, DIMENSION(:),     ALLOCATABLE, PROTECTED :: list   ! list(n)
  INTEGER, DIMENSION(:,:),   ALLOCATABLE, PROTECTED :: c      ! c(3,n) 3D cell index of each atom
  INTEGER, DIMENSION(:),     ALLOCATABLE, PROTECTED :: j_list ! j_list(n) list of j-partners
  INTEGER,                                PROTECTED :: nj     ! number of j-partners

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options

CONTAINS

  SUBROUTINE initialize_list ( n, r_cut )
    INTEGER, INTENT(in) :: n
    REAL,    INTENT(in) :: r_cut ! We assume that r_cut (in box=1 units) never changes

    WRITE(*,'(''Link list initialization'')')
    WRITE(*,'(''Link cells based on r_cut ='',t40,f15.5)') r_cut
    nc = FLOOR ( 1.0 / r_cut ) ! number of cells
    IF ( nc < 3 ) STOP 'System is too small to use link cells'

    ALLOCATE ( list(n), j_list(n), c(3,n) )
    ALLOCATE ( head(0:nc-1,0:nc-1,0:nc-1) )

  END SUBROUTINE initialize_list

  SUBROUTINE finalize_list
    DEALLOCATE ( list, j_list, c )
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
       c(:,i)  = c_index ( r(:,i) )         ! Index function
       list(i) = head(c(1,i),c(2,i),c(3,i)) ! store old head
       head(c(1,i),c(2,i),c(3,i)) = i       ! replace with current atom index
    END DO

  END SUBROUTINE make_list

  FUNCTION c_index ( ri ) RESULT ( ci )
    INTEGER, DIMENSION(3)             :: ci ! 3D cell index in range 0..nc-1
    REAL,    DIMENSION(3), INTENT(in) :: ri ! position in box = 1 units

    ci(:) = FLOOR ( ( 0.5 + ( ri(:) - ANINT(ri(:)) ) ) * REAL(nc) ) ! The index formula

    WHERE ( ci(:) < 0    ) ci(:) = 0    ! guard against small chance of roundoff error
    WHERE ( ci(:) > nc-1 ) ci(:) = nc-1 ! guard against small chance of roundoff error

  END FUNCTION c_index
  
  SUBROUTINE update_list ( i, ri ) ! Routine to update list structures when particle i has been moved
    INTEGER,               INTENT(in) :: i
    REAL,    DIMENSION(3), INTENT(in) :: ri

    INTEGER               :: this, next
    INTEGER, DIMENSION(3) :: ci

    ci(:) = c_index ( ri(:) ) ! Index function

    IF ( ALL ( ci(:) == c(:,i) ) ) RETURN ! no need to do anything

    ! Find i and remove it from c(:,i)
    this = head(c(1,i),c(2,i),c(3,i))
    IF ( this == i ) THEN ! i is the head atom in that cell

       head(c(1,i),c(2,i),c(3,i)) = list(i) ! simply point head at next atom

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
    
    ! Add i to the new cell
    c(:,i)  = ci(:)
    list(i) = head(c(1,i),c(2,i),c(3,i))
    head(c(1,i),c(2,i),c(3,i)) = i 

  END SUBROUTINE update_list

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

  SUBROUTINE get_neighbours ( i, op )

    ! Arguments
    INTEGER, intent(in) :: i  ! particle whose neighbours are required
    INTEGER, intent(in) :: op ! ne or gt, determining the range of neighbours

    ! Local variables
    INTEGER :: k1, k2, k, j
    INTEGER, DIMENSION(3) :: cj

    SELECT CASE ( op )
    CASE ( ne ) ! check every other atom in all cells
       k1 = -nk
       k2 =  nk
    CASE ( gt ) ! check half neighbour cells and j downlist from i in current cell
       k1 = 0
       k2 = nk
    CASE default
       STOP 'This should never happen'
    END SELECT

    nj = 0 ! Will store number of neighbours found

    DO k = k1, k2 ! Begin loop over neighbouring cells

       cj(:) = MODULO ( c(:,i) + dc(:,k), nc )

       IF ( k == 0 .AND. op == gt ) THEN
          j = list(i) ! check down-list from i in i-cell
       ELSE
          j = head(cj(1),cj(2),cj(3)) ! check entire j-cell
       END IF

       DO ! Begin loop over j atoms in list

          IF ( j == 0 ) EXIT
          
          IF ( j /= i ) THEN
             nj         = nj + 1 ! increment count of j atoms
             j_list(nj) = j      ! store new j atom
          END IF
          j = list(j) ! Next atom in j cell

       ENDDO ! End loop over j atoms in list

    ENDDO ! End loop over neighbouring cells 

  END SUBROUTINE get_neighbours

END MODULE link_list_module
