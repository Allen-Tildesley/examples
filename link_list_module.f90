! link_list_module.f90
! Link list handling routines for MC or MD simulation
MODULE link_list_module
  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: initialize_list, finalize_list, make_list, check_list
  PUBLIC :: move_in_list, create_in_list, destroy_in_list, c_index
  PUBLIC :: sc, head, list, c

  INTEGER,                                PROTECTED :: sc   ! dimensions of head array, assume cubic box
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE, PROTECTED :: head ! head(0:sc-1,0:sc-1,0:sc-1), assume cubic box
  INTEGER, DIMENSION(:),     ALLOCATABLE, PROTECTED :: list ! list(n)
  INTEGER, DIMENSION(:,:),   ALLOCATABLE, PROTECTED :: c    ! c(3,n) 3D cell index of each atom

CONTAINS

  SUBROUTINE initialize_list ( n, r_cut_box ) ! Routine to allocate list arrays
    INTEGER, INTENT(in) :: n         ! Number of particles
    REAL,    INTENT(in) :: r_cut_box ! rcut/box, assume never changes

    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Link cells based on r_cut/box =', r_cut_box

    sc = FLOOR ( 1.0 / r_cut_box ) ! Number of cells in each dimension
    IF ( sc < 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)') 'System is too small to use link cells', sc
       STOP 'Error in initialize_list'
    END IF

    ALLOCATE ( list(n), c(3,n) )
    ALLOCATE ( head(0:sc-1,0:sc-1,0:sc-1) )

  END SUBROUTINE initialize_list

  SUBROUTINE finalize_list ! Routine to deallocate list arrays
    DEALLOCATE ( list, c )
    DEALLOCATE ( head )
  END SUBROUTINE finalize_list

  SUBROUTINE make_list ( n, r ) ! Routine to make list
    INTEGER,                 INTENT(in) :: n ! Number of atoms
    REAL,    DIMENSION(3,n), INTENT(in) :: r ! Atom coordinates

    INTEGER :: i

    head(:,:,:) = 0

    DO i = 1, n ! Loop over all atoms
       c(:,i) = c_index ( r(:,i) )       ! Index function. allocating i to cell
       CALL create_in_list ( i, c(:,i) ) ! This does the work of adding i to list
    END DO ! End loop over all atoms

  END SUBROUTINE make_list

  FUNCTION c_index ( ri ) RESULT ( ci )
    INTEGER, DIMENSION(3)             :: ci ! Returns 3D cell index in range 0..sc-1, calculated from
    REAL,    DIMENSION(3), INTENT(in) :: ri ! position in box = 1 units

    ! We do not want to do any periodic imaging here, so as to cope with 
    ! Lees-Edwards boundaries as well as normal ones
    ! But we must check that ri is within bounds
    IF ( ANY ( ABS(ri) > 0.5 ) ) THEN ! Should never happen
       WRITE ( unit=error_unit, fmt='(a,3f15.5)') 'Atom not in main box', ri
       STOP 'Error in c_index'
    END IF
    
    ci(:) = FLOOR ( ( ri(:) + 0.5 ) * REAL(sc) ) ! The index formula

    ! Guard against small chance of roundoff error
    WHERE ( ci(:) < 0    ) ci(:) = 0
    WHERE ( ci(:) > sc-1 ) ci(:) = sc-1

  END FUNCTION c_index

  SUBROUTINE create_in_list ( i, ci ) ! Routine to create i in cell ci
    INTEGER,               INTENT(in) :: i  ! index of atom
    INTEGER, DIMENSION(3), INTENT(in) :: ci ! 3D index of cell in which i lies

    list(i)                 = head(ci(1),ci(2),ci(3)) ! transfer old head to list
    head(ci(1),ci(2),ci(3)) = i                       ! i becomes new head for this list
    c(:,i)                  = ci(:)                   ! store 3D index in array

  END SUBROUTINE create_in_list

  SUBROUTINE destroy_in_list ( i, ci ) ! Routine to destroy i in cell ci
    INTEGER,               INTENT(in) :: i  ! index of atom
    INTEGER, DIMENSION(3), INTENT(in) :: ci ! 3D index of cell in which i lies

    INTEGER :: this, next

    this = head(ci(1),ci(2),ci(3)) ! Locate head of list corresponding to cell

    IF ( this == i ) THEN ! i is the head atom in that cell

       head(ci(1),ci(2),ci(3)) = list(i) ! Simply point head at next atom, we're done

    ELSE ! i lies further down the list

       DO ! Loop traversing link-list
          next = list(this) ! Look ahead to the next entry

          IF ( next == i ) THEN ! Found our atom, just link over it
             list(this) = list(i)
             EXIT ! Leave the loop

          ELSE IF ( next == 0 ) THEN ! This should never happen
             WRITE ( unit=error_unit, fmt='(a,4i15)') 'Could not find particle in its cell', i, ci
             STOP 'Error in destroy_in_list'

          ELSE ! Move on to the next
             this = next ! Keep this index for next iteration
          END IF

       END DO ! End loop traversing link-list

    END IF

  END SUBROUTINE destroy_in_list
  
  SUBROUTINE move_in_list ( i, ci ) ! Routine to move i from current cell to ci
    INTEGER,               INTENT(in) :: i
    INTEGER, DIMENSION(3), INTENT(in) :: ci

    IF ( ALL ( ci(:) == c(:,i) ) ) RETURN ! no need to do anything

    CALL destroy_in_list ( i, c(:,i) ) ! Remove from old cell
    CALL create_in_list  ( i, ci(:)  ) ! Add to new cell

  END SUBROUTINE move_in_list

  SUBROUTINE check_list ( n, r ) ! Routine to check consistency of cell lists
    INTEGER,                 INTENT(in) :: n ! Number of atoms
    REAL,    DIMENSION(3,n), INTENT(in) :: r ! Atom positions

    INTEGER, DIMENSION(3) :: ci
    INTEGER               :: ci1, ci2, ci3, i

    DO i = 1, n ! Loop to check that each atom's cell is correct
       ci = c_index ( r(:,i) ) ! Index function
       IF ( ANY ( ci(:) /= c(:,i) ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,7i10)') 'Inconsistency 1 found:', i, ci, c(:,i)
          STOP 'Error in check_list'
       END IF
    END DO ! End loop to check that each atom's cell is correct

    ! Triple loop over cells
    DO ci1 = 0, sc-1
       DO ci2 = 0, sc-1
          DO ci3 = 0, sc-1

             ci = [ ci1, ci2, ci3 ] ! Store 3D cell index in array
             i = head(ci1,ci2,ci3)  ! Locate head of list corresponding to cell

             DO ! Loop traversing link list
                IF ( i == 0 ) EXIT ! Reached end of list

                IF ( ANY ( ci(:) /= c(:,i) ) ) THEN
                   WRITE ( unit=error_unit, fmt='(a,7i10)') 'Inconsistency 2 found:', i, ci, c(:,i)
                   STOP 'Error in check_list'
                END IF

                i = list(i) ! Move on to next list entry
             END DO ! End loop traversing link list

          END DO
       END DO
    END DO
    ! End triple loop over cells

  END SUBROUTINE check_list

END MODULE link_list_module
