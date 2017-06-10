! md_lj_ll_module.f90
! Force routine for MD simulation, LJ atoms, using link-lists
MODULE md_module

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
  PUBLIC :: force, hessian

  ! Public data
  INTEGER,                              PUBLIC :: n ! Number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v ! Velocities (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: f ! Forces (3,n)

  ! Public derived type
  TYPE, PUBLIC :: potential_type ! A composite variable for interactions comprising
     REAL    :: cut ! the potential energy cut (but not shifted) at r_cut and
     REAL    :: pot ! the potential energy cut-and-shifted at r_cut and
     REAL    :: vir ! the virial and
     REAL    :: lap ! the Laplacian and
     LOGICAL :: ovr ! a flag indicating overlap (i.e. pot too high to use)
   CONTAINS
     PROCEDURE :: add_potential_type
     GENERIC   :: OPERATOR(+) => add_potential_type
  END TYPE potential_type

CONTAINS

  FUNCTION add_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the sum of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%cut = a%cut  +   b%cut
    c%pot = a%pot  +   b%pot
    c%vir = a%vir  +   b%vir
    c%lap = a%lap  +   b%lap
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut-and-shifted version for dynamics'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut (but not shifted) version also calculated'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1' 
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
    USE link_list_module, ONLY : initialize_list
    IMPLICIT NONE
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

    REAL :: r_cut_box

    ALLOCATE ( r(3,n), v(3,n), f(3,n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)' ) 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

    CALL initialize_list ( n, r_cut_box )

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    USE link_list_module, ONLY : finalize_list
    IMPLICIT NONE

    DEALLOCATE ( r, v, f )
    CALL finalize_list

  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( box, r_cut, total )
    USE link_list_module, ONLY : make_list, sc, head, list
    IMPLICIT NONE
    REAL,                 INTENT(in)  :: box   ! Simulation box length
    REAL,                 INTENT(in)  :: r_cut ! Potential cutoff distance
    TYPE(potential_type), INTENT(out) :: total ! Composite of pot, vir, lap etc

    ! total%pot is the nonbonded cut-and-shifted potential energy for whole system
    ! total%cut is the nonbonded cut (but not shifted) potential energy for whole system
    ! total%vir is the corresponding virial
    ! total%lap is the corresponding Laplacian
    ! total%ovr is a warning flag that there is an overlap
    ! This routine also calculates forces and stores them in the array f
    ! Forces are derived from pot, not cut (which has a discontinuity)
    ! If total%ovr is set to .true., the forces etc should not be used
    ! Uses link lists

    INTEGER               :: i, j, ci1, ci2, ci3, k
    INTEGER, DIMENSION(3) :: ci, cj
    REAL                  :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL                  :: sr2, sr6, sr12, pot_cut
    REAL,    DIMENSION(3) :: rij, fij
    REAL,    PARAMETER    :: sr2_ovr = 1.77 ! Overlap threshold (pot > 100)
    TYPE(potential_type)  :: pair

    ! Set up vectors to half the cells in neighbourhood of 3x3x3 cells in cubic lattice
    ! The cells are chosen so that if (d1,d2,d3) appears, then (-d1,-d2,-d3) does not.
    INTEGER, PARAMETER :: nk = 13 
    INTEGER, DIMENSION(3,0:nk), PARAMETER :: d = RESHAPE( [ &
         &                0, 0, 0,    1, 0, 0, &
         &    1, 1, 0,   -1, 1, 0,    0, 1, 0, &
         &    0, 0, 1,   -1, 0, 1,    1, 0, 1, &
         &   -1,-1, 1,    0,-1, 1,    1,-1, 1, &
         &   -1, 1, 1,    0, 1, 1,    1, 1, 1    ], [ 3, nk+1 ] )

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    ! Calculate potential at cutoff
    sr2     = 1.0 / r_cut**2 ! in sigma=1 units
    sr6     = sr2 ** 3
    sr12    = sr6 **2
    pot_cut = sr12 - sr6 ! Without numerical factor 4

    r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundary conditions in box=1 units
    CALL make_list ( n, r )

    ! Initialize
    f     = 0.0
    total = potential_type (  pot=0.0, cut=0.0, vir=0.0, lap=0.0, ovr=.FALSE. )

    ! Triple loop over cells
    DO ci1 = 0, sc-1
       DO ci2 = 0, sc-1
          DO ci3 = 0, sc-1

             ci(:) = [ ci1, ci2, ci3 ] ! 3D index of i-cell
             i = head(ci1,ci2,ci3)     ! First i-atom in cell

             DO ! Begin loop over i-atoms in list
                IF ( i == 0 ) EXIT ! End of link list

                DO k = 0, nk ! Loop over neighbouring cells

                   IF ( k == 0 ) THEN
                      j = list(i) ! First j-atom is downlist from i in current cell
                   ELSE
                      cj(:) = ci(:) + d(:,k)          ! Neighbour j-cell 3D index
                      cj(:) = MODULO ( cj(:), sc )    ! Periodic boundary correction
                      j     = head(cj(1),cj(2),cj(3)) ! First j-atom in neighbour cell
                   END IF

                   DO ! Begin loop over j-atoms in list
                      IF ( j == 0 ) EXIT ! End of link list

                      IF ( j == i ) THEN ! This should never happen
                         WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Index error', i, j
                         STOP 'Impossible error in force' 
                      END IF

                      rij(:) = r(:,i) - r(:,j)           ! Separation vector
                      rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundary conditions in box=1 units
                      rij_sq = SUM ( rij**2 )            ! Squared separation

                      IF ( rij_sq < r_cut_box_sq ) THEN ! Check within cutoff

                         rij_sq   = rij_sq * box_sq ! Now in sigma=1 units
                         rij(:)   = rij(:) * box    ! Now in sigma=1 units
                         sr2      = 1.0 / rij_sq    ! (sigma/rij)**2
                         pair%ovr = sr2 > sr2_ovr   ! Overlap if too close

                         sr6      = sr2 ** 3
                         sr12     = sr6 ** 2
                         pair%cut = sr12 - sr6                    ! LJ pair potential (cut but not shifted)
                         pair%vir = pair%cut + sr12               ! LJ pair virial
                         pair%pot = pair%cut - pot_cut            ! LJ pair potential (cut-and-shifted)
                         pair%lap = ( 22.0*sr12 - 5.0*sr6 ) * sr2 ! LJ pair Laplacian
                         fij      = rij * pair%vir * sr2          ! LJ pair forces

                         total  = total  + pair
                         f(:,i) = f(:,i) + fij
                         f(:,j) = f(:,j) - fij

                      END IF ! End check within cutoff

                      j = list(j) ! Next j-atom
                   END DO ! End loop over j-atoms in list

                END DO ! End loop over neighbouring cells

                i = list(i) ! Next i-atom
             END DO ! End loop over i-atoms in list

          END DO
       END DO
    END DO
    ! End triple loop over cells

    ! Multiply results by numerical factors
    f         = f         * 24.0       ! 24*epsilon
    total%cut = total%cut * 4.0        ! 4*epsilon
    total%pot = total%pot * 4.0        ! 4*epsilon
    total%vir = total%vir * 24.0 / 3.0 ! 24*epsilon and divide virial by 3
    total%lap = total%lap * 24.0 * 2.0 ! 24*epsilon and factor 2 for ij and ji

  END SUBROUTINE force

  FUNCTION hessian ( box, r_cut ) RESULT ( hes )
    USE link_list_module, ONLY : make_list, sc, head, list
    IMPLICIT NONE

    REAL, INTENT(in)  :: box    ! Simulation box length
    REAL, INTENT(in)  :: r_cut  ! Potential cutoff distance
    REAL              :: hes    ! Result

    ! Calculates Hessian function (for 1/N correction to config temp)
    ! This routine is only needed in a constant-energy ensemble
    ! It is assumed that positions are in units where box = 1
    ! but the result is given in units where sigma = 1 and epsilon = 1
    ! It is assumed that forces have already been calculated 
    ! Uses link lists

    INTEGER               :: i, j, ci1, ci2, ci3, k
    INTEGER, DIMENSION(3) :: ci, cj
    REAL                  :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL                  :: sr2, sr6, sr8, sr10, rf, ff, v1, v2
    REAL,    DIMENSION(3) :: rij, fij

    ! Set up vectors to half the cells in neighbourhood of 3x3x3 cells in cubic lattice
    ! The cells are chosen so that if (d1,d2,d3) appears, then (-d1,-d2,-d3) does not.
    INTEGER, PARAMETER :: nk = 13 
    INTEGER, DIMENSION(3,0:nk), PARAMETER :: d = RESHAPE( [ &
         &                0, 0, 0,    1, 0, 0, &
         &    1, 1, 0,   -1, 1, 0,    0, 1, 0, &
         &    0, 0, 1,   -1, 0, 1,    1, 0, 1, &
         &   -1,-1, 1,    0,-1, 1,    1,-1, 1, &
         &   -1, 1, 1,    0, 1, 1,    1, 1, 1    ], [ 3, nk+1 ] )

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    hes = 0.0

    r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundary conditions in box=1 units
    CALL make_list ( n, r )

    ! Triple loop over cells
    DO ci1 = 0, sc-1
       DO ci2 = 0, sc-1
          DO ci3 = 0, sc-1

             ci(:) = [ ci1, ci2, ci3 ] ! 3D index of i-cell
             i = head(ci1,ci2,ci3)     ! First i-atom in cell

             DO ! Begin loop over i-atoms in list
                IF ( i == 0 ) EXIT ! End of link list

                DO k = 0, nk ! Loop over neighbouring cells

                   IF ( k == 0 ) THEN
                      j = list(i) ! First j-atom is downlist from i in current cell
                   ELSE
                      cj(:) = ci(:) + d(:,k)          ! Neighbour j-cell 3D index
                      cj(:) = MODULO ( cj(:), sc )    ! Periodic boundary correction
                      j     = head(cj(1),cj(2),cj(3)) ! First j-atom in neighbour cell
                   END IF

                   DO ! Begin loop over j-atoms in list
                      IF ( j == 0 ) EXIT ! End of link list

                      IF ( j == i ) THEN ! This should never happen
                         WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Index error', i, j
                         STOP 'Impossible error in hessian' 
                      END IF

                      rij(:) = r(:,i) - r(:,j)           ! Separation vector
                      rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundary conditions in box=1 units
                      rij_sq = SUM ( rij**2 )            ! Squared separation

                      IF ( rij_sq < r_cut_box_sq ) THEN ! Check within cutoff

                         rij_sq = rij_sq * box_sq ! Now in sigma=1 units
                         rij(:) = rij(:) * box    ! Now in sigma=1 units
                         fij(:) = f(:,i) - f(:,j) ! Difference in forces

                         ff   = DOT_PRODUCT(fij,fij)
                         rf   = DOT_PRODUCT(rij,fij)
                         sr2  = 1.0 / rij_sq
                         sr6  = sr2 ** 3
                         sr8  = sr6 * sr2
                         sr10 = sr8 * sr2
                         v1   = 24.0 * ( 1.0 - 2.0 * sr6 ) * sr8
                         v2   = 96.0 * ( 7.0 * sr6 - 2.0 ) * sr10
                         hes  = hes + v1 * ff + v2 * rf**2

                      END IF ! End check within cutoff

                      j = list(j) ! Next j-atom
                   END DO ! End loop over j-atoms in list

                END DO ! End loop over neighbouring cells

                i = list(i) ! Next i-atom
             END DO ! End loop over i-atoms in list

          END DO
       END DO
    END DO
    ! End triple loop over cells

  END FUNCTION hessian

END MODULE md_module

