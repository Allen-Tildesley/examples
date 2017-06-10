! link_list_module.f90
! Link list handling routines for MC or MD simulation
MODULE link_list_module

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
  PUBLIC :: initialize_list, finalize_list, make_list, check_list
  PUBLIC :: move_in_list, create_in_list, destroy_in_list, c_index, neighbours

  ! Public (protected) data
  INTEGER,                                PROTECTED, SAVE, PUBLIC :: sc     ! dimensions of head array
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE, PROTECTED, SAVE, PUBLIC :: head   ! head(0:sc-1,0:sc-1,0:sc-1)
  INTEGER, DIMENSION(:),     ALLOCATABLE, PROTECTED, SAVE, PUBLIC :: list   ! list(n)
  INTEGER, DIMENSION(:,:),   ALLOCATABLE, PROTECTED, SAVE, PUBLIC :: c      ! c(3,n) 3D cell index of each atom

CONTAINS

  SUBROUTINE initialize_list ( n, r_cut_box ) ! Routine to allocate list arrays
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n         ! Number of particles
    REAL,    INTENT(in) :: r_cut_box ! rcut/box, assume never changes

    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Link cells based on r_cut/box =', r_cut_box

    sc = FLOOR ( 1.0 / r_cut_box ) ! Number of cells in each dimension
    IF ( sc < 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)') 'System is too small to use link cells', sc
       STOP 'Error in initialize_list'
    END IF

    ALLOCATE ( list(n), c(3,n) )
    ALLOCATE ( head(0:sc-1,0:sc-1,0:sc-1) )

  END SUBROUTINE initialize_list

  SUBROUTINE finalize_list ! Routine to deallocate list arrays
    IMPLICIT NONE

    DEALLOCATE ( list, c )
    DEALLOCATE ( head )

  END SUBROUTINE finalize_list

  SUBROUTINE make_list ( n, r ) ! Routine to make list
    IMPLICIT NONE
    INTEGER,                 INTENT(in) :: n ! Number of atoms
    REAL,    DIMENSION(3,n), INTENT(in) :: r ! Atom coordinates

    INTEGER :: i
    
    head(:,:,:) = 0

    DO i = 1, n ! Loop over all atoms
       c(:,i) = c_index ( r(:,i) )       ! Index function allocating atom i to cell
       CALL create_in_list ( i, c(:,i) ) ! This does the work of adding atom i to list
    END DO ! End loop over all atoms

  END SUBROUTINE make_list

  FUNCTION c_index ( ri ) RESULT ( ci )
    IMPLICIT NONE
    INTEGER, DIMENSION(3)             :: ci ! Returns 3D cell index in range 0..sc-1, calculated from
    REAL,    DIMENSION(3), INTENT(in) :: ri ! position in box = 1 units

    ! We do not want to do any periodic imaging here, so as to cope with 
    ! Lees-Edwards boundaries as well as normal ones
    ! But we must check that ri is within bounds
    IF ( ANY ( ABS(ri) > 0.5 ) ) THEN ! Should never happen
       WRITE ( unit=error_unit, fmt='(a,3f15.6)') 'Atom not in main box', ri
       STOP 'Error in c_index'
    END IF

    ci(:) = FLOOR ( ( ri(:) + 0.5 ) * REAL(sc) ) ! The index formula

    ! Guard against small chance of roundoff error
    WHERE ( ci(:) < 0    ) ci(:) = 0
    WHERE ( ci(:) > sc-1 ) ci(:) = sc-1

  END FUNCTION c_index

  SUBROUTINE create_in_list ( i, ci ) ! Routine to create atom i in cell ci
    IMPLICIT NONE
    INTEGER,               INTENT(in) :: i  ! Index of atom
    INTEGER, DIMENSION(3), INTENT(in) :: ci ! 3D index of cell in which i lies
    
    list(i)                 = head(ci(1),ci(2),ci(3)) ! Transfer old head to list
    head(ci(1),ci(2),ci(3)) = i                       ! Atom i becomes new head for this list
    c(:,i)                  = ci(:)                   ! Store 3D index in array

  END SUBROUTINE create_in_list

  SUBROUTINE destroy_in_list ( i, ci ) ! Routine to destroy atom i in cell ci
    IMPLICIT NONE
    INTEGER,               INTENT(in) :: i  ! Index of atom
    INTEGER, DIMENSION(3), INTENT(in) :: ci ! 3D index of cell in which i lies

    INTEGER :: this, next
    
    this = head(ci(1),ci(2),ci(3)) ! Locate head of list corresponding to cell

    IF ( this == i ) THEN ! Atom i is the head atom in that cell

       head(ci(1),ci(2),ci(3)) = list(i) ! Simply point head at next atom, we're done

    ELSE ! Atom i lies further down the list

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

  SUBROUTINE move_in_list ( i, ci ) ! Routine to move atom i from current cell to ci
    IMPLICIT NONE
    INTEGER,               INTENT(in) :: i
    INTEGER, DIMENSION(3), INTENT(in) :: ci
    
    IF ( ALL ( ci(:) == c(:,i) ) ) RETURN ! No need to do anything

    CALL destroy_in_list ( i, c(:,i) ) ! Remove atom i from old cell
    CALL create_in_list  ( i, ci(:)  ) ! Add atom i to new cell

  END SUBROUTINE move_in_list

  SUBROUTINE check_list ( n, r ) ! Routine to check consistency of cell lists
    IMPLICIT NONE
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

  FUNCTION neighbours ( n, i, ci, half ) RESULT ( j_list )
    IMPLICIT NONE
    INTEGER,               INTENT(in) :: n      ! Number of atoms
    INTEGER,               INTENT(in) :: i      ! Atom whose neighbours are required
    INTEGER, DIMENSION(3), INTENT(in) :: ci     ! Cell of atom of interest
    LOGICAL,               INTENT(in) :: half   ! Determining the range of neighbours searched
    INTEGER, DIMENSION(n)             :: j_list ! Resulting list of indices

    ! This routine uses the link-list cell structure to fill out the array j_list
    ! with possible neighbours of atom i, padding with zeroes
    ! If half==.false., cell ci and all 26 surrounding cells are searched.
    ! If half==.true., cell ci, and just 13 of the neighbour cells, are searched
    ! and moreover, in ci, we only look down-list making use of list(i)
    ! There is a subtlety: using list(i) assumes that our interest is in the cells that
    ! are neighbours of c(:,i), i.e. that ci(:)==c(:,i), and we check for this explicitly.
    ! In other words, we assume that atom i has not moved since list was constructed.
    ! If half==.false., particle i might be in a very different position, and ci might be
    ! very different to c(:,i) but in that case we make no use of list(i), in normal use

    ! We have a cubic cell lattice
    ! Set up vectors to each cell in the 3x3x3 neighbourhood of a given cell
    ! To work properly, these are listed with inversion symmetry about (0,0,0)
    INTEGER,                      PARAMETER :: nk = 13 
    INTEGER, DIMENSION(3,-nk:nk), PARAMETER :: d = RESHAPE( [ &
         &   -1,-1,-1,    0,-1,-1,    1,-1,-1, &
         &   -1, 1,-1,    0, 1,-1,    1, 1,-1, &
         &   -1, 0,-1,    1, 0,-1,    0, 0,-1, &
         &    0,-1, 0,    1,-1, 0,   -1,-1, 0, &
         &   -1, 0, 0,    0, 0, 0,    1, 0, 0, &
         &    1, 1, 0,   -1, 1, 0,    0, 1, 0, &
         &    0, 0, 1,   -1, 0, 1,    1, 0, 1, &
         &   -1,-1, 1,    0,-1, 1,    1,-1, 1, &
         &   -1, 1, 1,    0, 1, 1,    1, 1, 1    ], [ 3, 2*nk+1 ] )

    INTEGER               :: k1, k2, k, j, nj
    INTEGER, DIMENSION(3) :: cj

    IF ( half ) THEN ! Check half neighbour cells and j downlist from i in current cell
       k1 = 0
       k2 = nk
       IF ( ANY ( ci(:) /= c(:,i) ) ) THEN ! should never happen
          WRITE ( unit=error_unit, fmt='(a,6i15)' ) 'Cell mismatch ', ci(:), c(:,i)
          STOP 'Error in get_neighbours'
       END IF
    ELSE ! Check every atom other than i in all cells
       k1 = -nk
       k2 =  nk
    END IF

    j_list = 0 ! Initialize with zero values everywhere
    nj     = 0 ! Next position in list to be filled

    DO k = k1, k2 ! Begin loop over neighbouring cells

       cj(:) = ci(:) + d(:,k)       ! Neighbour cell index
       cj(:) = MODULO ( cj(:), sc ) ! Periodic boundary correction

       IF ( k == 0 .AND. half ) THEN
          j = list(i) ! Check down-list from i in i-cell
       ELSE
          j = head(cj(1),cj(2),cj(3)) ! Check entire j-cell
       END IF

       DO ! Begin loop over j atoms in list

          IF ( j == 0 ) EXIT  ! Exhausted list

          IF ( j /= i ) THEN ! Skip self

             nj = nj + 1 ! Increment count of j atoms

             IF ( nj >= n ) THEN ! Check more than n-1 neighbours (should never happen)
                WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Neighbour error for j_list', nj, n
                STOP 'Impossible error in get_neighbours'
             END IF ! End check more than n-1 neighbours

             j_list(nj) = j       ! Store new j atom
          END IF

          j = list(j) ! Next atom in j cell

       END DO ! End loop over j atoms in list

    END DO ! End loop over neighbouring cells 

  END FUNCTION neighbours

END MODULE link_list_module
