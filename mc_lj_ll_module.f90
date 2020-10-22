! mc_lj_ll_module.f90
! Energy and move routines for MC, LJ potential, link-lists
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

  ! This (together with link_list_module.f90) functions as a drop-in replacement
  ! for mc_lj_module.f90 which should work in our constant-NVT NPT and zVT programs.

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, extend_arrays, deallocate_arrays
  PUBLIC :: box_check, potential_1, potential, force_sq
  PUBLIC :: move, create, destroy

  ! Public data
  INTEGER,                              PUBLIC :: n ! Number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n)

  ! Private data
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f      ! Forces for force_sq calculation (3,n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: j_list ! List of j-neighbours (n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: tmp_r  ! Temporary array for dynamic reallocation of r

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range

  ! Public derived type
  TYPE, PUBLIC :: potential_type ! A composite variable for interactions comprising
     REAL    :: pot ! the potential energy cut (but not shifted) at r_cut, and
     REAL    :: vir ! the virial and
     REAL    :: lap ! the Laplacian and
     LOGICAL :: ovr ! a flag indicating overlap (i.e. pot too high to use)
   CONTAINS
     PROCEDURE :: add_potential_type
     PROCEDURE :: subtract_potential_type
     GENERIC   :: OPERATOR(+) => add_potential_type
     GENERIC   :: OPERATOR(-) => subtract_potential_type
  END TYPE potential_type

CONTAINS

  FUNCTION add_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the sum of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot  +   b%pot
    c%vir = a%vir  +   b%vir
    c%lap = a%lap  +   b%lap
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  FUNCTION subtract_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the difference of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot  -   b%pot
    c%vir = a%vir  -   b%vir
    c%lap = a%lap  -   b%lap
    c%ovr = a%ovr .OR. b%ovr ! This is meaningless, but inconsequential
  END FUNCTION subtract_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut (but not shifted)'
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

    ALLOCATE ( r(3,n), f(3,n), j_list(n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)') 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

    CALL initialize_list ( n, r_cut_box )

  END SUBROUTINE allocate_arrays

  SUBROUTINE extend_arrays ( n_new )
    USE link_list_module, ONLY : extend_list
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n_new ! New size of arrays

    INTEGER :: n_old

    n_old = SIZE ( r, dim=2 )

    IF ( n_new <= n_old ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)') 'Unexpected array size ', n_old, n_new
       STOP 'Error in extend_arrays'
    END IF

    ! Reallocate r array, preserving existing data
    ALLOCATE ( tmp_r(3,n_new) )
    tmp_r(:,1:n_old)  = r
    tmp_r(:,n_old+1:) = 0.0
    CALL MOVE_ALLOC ( tmp_r, r )

    ! Reallocate f array, no need to preserve data
    DEALLOCATE ( f )
    ALLOCATE ( f(3,n_new) )

    ! Reallocate j_list array, no need to preserve data
    DEALLOCATE ( j_list )
    ALLOCATE ( j_list(n_new) )

    CALL extend_list ( n_new )

  END SUBROUTINE extend_arrays

  SUBROUTINE deallocate_arrays
    USE link_list_module, ONLY : finalize_list
    IMPLICIT NONE

    DEALLOCATE ( r, f, j_list )
    CALL finalize_list

  END SUBROUTINE deallocate_arrays

  SUBROUTINE box_check ( box, r_cut )
    USE link_list_module, ONLY : make_list
    IMPLICIT NONE
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

    REAL :: r_cut_box

    r_cut_box = r_cut / box
    CALL make_list ( n, r_cut_box, r )

  END SUBROUTINE box_check

  FUNCTION potential ( box, r_cut ) RESULT ( total )
    USE link_list_module, ONLY : make_list
    IMPLICIT NONE
    TYPE(potential_type) :: total ! Returns a composite of pot, vir etc
    REAL, INTENT(in)     :: box   ! Simulation box length
    REAL, INTENT(in)     :: r_cut ! Potential cutoff

    ! total%pot is the nonbonded cut (not shifted) potential energy for whole system
    ! total%vir is the corresponding virial for whole system
    ! total%lap is the corresponding Laplacian for whole system
    ! total%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of total%pot etc should not be used
    ! Actual calculation is performed by function potential_1

    ! This routine uses linked lists.
    ! We assume that the box length may vary, so we call the make_list function at the start.
    ! The list for fixed box length is maintained by create_in_list, destroy_in_list,
    ! and move_in_list which are called in other routines below.

    TYPE(potential_type) :: partial
    INTEGER              :: i
    REAL                 :: r_cut_box

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Impossible error in potential'
    END IF

    r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
    r_cut_box = r_cut / box
    CALL make_list ( n, r_cut_box, r )

    total = potential_type ( pot=0.0, vir=0.0, lap=0.0, ovr=.FALSE. ) ! Initialize

    DO i = 1, n

       partial = potential_1 ( r(:,i), i, box, r_cut, gt )

       IF ( partial%ovr ) THEN
          total%ovr = .TRUE. ! Overlap detected
          RETURN             ! Return immediately
       END IF

       total = total + partial

    END DO

    total%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential

  FUNCTION potential_1 ( ri, i, box, r_cut, j_range ) RESULT ( partial )
    USE link_list_module, ONLY : c, c_index, neighbours
    IMPLICIT NONE
    TYPE(potential_type)           :: partial ! Returns a composite of pot, vir etc for given atom
    REAL, DIMENSION(3), INTENT(in) :: ri      ! Coordinates of atom of interest
    INTEGER,            INTENT(in) :: i       ! Index of atom of interest
    REAL,               INTENT(in) :: box     ! Simulation box length
    REAL,               INTENT(in) :: r_cut   ! Potential cutoff distance
    INTEGER, OPTIONAL,  INTENT(in) :: j_range ! Optional partner index range

    ! partial%pot is the nonbonded cut (not shifted) potential energy of atom ri with a set of other atoms
    ! partial%vir is the corresponding virial of atom ri
    ! partial%lap is the corresponding Laplacian of atom ri
    ! partial%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the value of partial%pot should not be used
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to "half" which is
    ! actually implemented in the neighbours routine

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    ! This routine uses linked lists.
    ! We assume that the box length may vary, but that the make_list function
    ! has been called since the last change in box length.

    INTEGER               :: j, jj
    LOGICAL               :: half
    REAL                  :: r_cut_box, r_cut_box_sq, box_sq
    REAL                  :: sr2, sr6, sr12, rij_sq
    REAL,    DIMENSION(3) :: rij
    INTEGER, DIMENSION(3) :: ci
    REAL, PARAMETER       :: sr2_ovr = 1.77 ! overlap threshold (pot > 100)
    TYPE(potential_type)  :: pair
    
    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in potential_1'
    END IF

    ci = c_index ( ri ) ! Cell in which ri lies (not necessarily same as c(:,i))

    half = PRESENT ( j_range)
    IF ( half .AND. ANY ( ci(:) /= c(:,i) ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,6i15)' ) 'Cell mismatch ', ci(:), c(:,i)
       STOP 'Error in potential_1'
    END IF

    j_list = neighbours ( n, i, ci, half ) ! Put neighbours in j_list

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box**2
    box_sq       = box**2

    partial = potential_type ( pot=0.0, vir=0.0, lap=0.0, ovr=.FALSE. ) ! Initialize

    jj = 0

    DO ! Loop until no more entries in j_list
       jj = jj + 1        ! Next entry
       j  = j_list(jj)    ! Get neighbour index
       IF ( j == 0 ) EXIT ! List exhausted

       IF ( i == j ) CYCLE ! Skip self ( should never happen)

       rij(:) = ri(:) - r(:,j)            ! Separation vector
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )            ! Squared separation in box=1 units

       IF ( rij_sq < r_cut_box_sq ) THEN ! Check within range

          rij_sq   = rij_sq * box_sq ! Now in sigma=1 units
          sr2      = 1.0 / rij_sq    ! (sigma/rij)**2
          pair%ovr = sr2 > sr2_ovr   ! Overlap if too close

          IF ( pair%ovr ) THEN
             partial%ovr = .TRUE. ! Overlap detected
             RETURN               ! Return immediately
          END IF

          sr6      = sr2**3
          sr12     = sr6**2
          pair%pot = sr12 - sr6                    ! LJ pair potential (cut but not shifted)
          pair%vir = pair%pot + sr12               ! LJ pair virial
          pair%lap = ( 22.0*sr12 - 5.0*sr6 ) * sr2 ! LJ pair Laplacian

          partial = partial + pair

       END IF ! End check within range

    END DO ! End loop until no more entries in j_list

    ! Numerical factors
    partial%pot = partial%pot * 4.0        ! 4*epsilon
    partial%vir = partial%vir * 24.0 / 3.0 ! 24*epsilon and divide virial by 3
    partial%lap = partial%lap * 24.0 * 2.0 ! 24*epsilon and factor 2 for ij and ji
    partial%ovr = .FALSE.                  ! No overlaps detected (redundant but for clarity)

  END FUNCTION potential_1

  FUNCTION force_sq ( box, r_cut ) RESULT ( fsq )
    USE link_list_module, ONLY : make_list, c, neighbours
    IMPLICIT NONE
    REAL             :: fsq   ! Returns total squared force
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

    ! Calculates total squared force (using array f)
    ! Uses link lists
    ! We assume that the box length may vary, so we call the make_list function at the start.

    INTEGER               :: i, j, jj
    REAL                  :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL                  :: sr2, sr6, sr12
    REAL,    DIMENSION(3) :: rij, fij
    
    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2
    CALL make_list ( n, r_cut_box, r )

    ! Initialize
    f = 0.0

    DO i = 1, n ! Outer loop over i

       j_list = neighbours ( n, i, c(:,i), half=.TRUE. ) ! Put neighbours in j_list

       jj = 0

       DO ! Loop until no more entries in j_list
          jj = jj + 1        ! Next entry
          j  = j_list(jj)    ! Get neighbour index
          IF ( j == 0 ) EXIT ! List exhausted

          IF ( i == j ) CYCLE ! Skip self ( should never happen)

          rij(:) = r(:,i) - r(:,j)           ! Separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )            ! Squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN ! Check within cutoff

             rij_sq = rij_sq * box_sq ! Now in sigma=1 units
             rij(:) = rij(:) * box    ! Now in sigma=1 units
             sr2    = 1.0 / rij_sq
             sr6   = sr2 ** 3
             sr12  = sr6 ** 2
             fij   = rij * (2.0*sr12 - sr6) / rij_sq ! LJ pair forces
             f(:,i) = f(:,i) + fij
             f(:,j) = f(:,j) - fij

          END IF ! End check within cutoff

       END DO ! End loop until no more entries in j_list

    END DO ! End outer loop over i

    f   = f * 24.0 ! Numerical factor
    fsq = SUM ( f**2 )

  END FUNCTION force_sq

  SUBROUTINE move ( i, ri )
    USE link_list_module, ONLY : c_index, move_in_list
    IMPLICIT NONE
    INTEGER,               INTENT(in) :: i
    REAL,    DIMENSION(3), INTENT(in) :: ri

    INTEGER, DIMENSION(3) :: ci
    
    r(:,i) = ri                ! New position
    ci(:)  = c_index ( ri(:) ) ! New cell index
    CALL move_in_list ( i, ci(:) )

  END SUBROUTINE move

  SUBROUTINE create ( ri )
    USE link_list_module, ONLY : c_index, create_in_list
    IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(in) :: ri

    INTEGER, DIMENSION(3) :: ci

    n      = n+1               ! Increase number of atoms
    r(:,n) = ri(:)             ! Add new atom at the end
    ci(:)  = c_index ( ri(:) ) ! New cell index
    CALL create_in_list ( n, ci )

  END SUBROUTINE create

  SUBROUTINE destroy ( i )
    USE link_list_module, ONLY : destroy_in_list, move_in_list, c
    IMPLICIT NONE
    INTEGER, INTENT(in) :: i

    r(:,i) = r(:,n) ! Replace atom i coordinates with atom n
    CALL destroy_in_list ( n, c(:,n) )
    CALL move_in_list ( i, c(:,n) )
    n = n - 1  ! Reduce number of atoms

  END SUBROUTINE destroy

END MODULE mc_module
