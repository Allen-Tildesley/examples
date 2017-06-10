! mc_gibbs_lj_module.f90
! Energy and move routines for Gibbs MC, LJ potential
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

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: potential_1, potential
  PUBLIC :: move, swap

  ! Public data
  INTEGER, DIMENSION(2),                PUBLIC :: n ! Number of atoms in each box
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n(1)+n(2))

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range

  ! Public derived type
  TYPE, PUBLIC :: potential_type   ! A composite variable for interactions comprising
     REAL    :: pot ! the potential energy cut at r_cut and
     REAL    :: vir ! the virial and
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
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  FUNCTION subtract_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the difference of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot  -   b%pot
    c%vir = a%vir  -   b%vir
    c%ovr = a%ovr .OR. b%ovr ! This is meaningless, but inconsequential
  END FUNCTION subtract_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut (but not shifted)'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'    
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Two simulation boxes'

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
    IMPLICIT NONE
    REAL, DIMENSION(2), INTENT(in) :: box   ! Simulation box lengths
    REAL,               INTENT(in) :: r_cut ! Potential cutoff distance

    REAL :: r_cut_box

    ALLOCATE ( r(3,SUM(n)) )

    r_cut_box = r_cut / MINVAL ( box )
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)') 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r )

  END SUBROUTINE deallocate_arrays

  FUNCTION potential ( i1, i2, box, r_cut ) RESULT ( total )
    IMPLICIT NONE
    TYPE(potential_type) :: total  ! Returns a composite of pot, vir etc
    INTEGER, INTENT(in)  :: i1, i2 ! Index range defining system: (1,n(1)) or (n(1)+1,n(1)+n(2))
    REAL,    INTENT(in)  :: box    ! Simulation box length
    REAL,    INTENT(in)  :: r_cut  ! Potential cutoff distance

    ! total%pot is the nonbonded cut (not shifted) potential energy for whole system
    ! total%vir is the corresponding virial for whole system
    ! total%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of total%pot etc should not be used
    ! Actual calculation is performed by function potential_1

    TYPE(potential_type) :: partial ! Atomic contribution to total
    INTEGER              :: i

    IF ( SUM(n) > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', SUM(n), SIZE(r,dim=2)
       STOP 'Error in potential'
    END IF

    total = potential_type ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO i = i1, i2

       partial = potential_1 ( i1, i2, r(:,i), i, box, r_cut, gt )

       IF ( partial%ovr ) THEN
          total%ovr = .TRUE. ! Overlap detected
          RETURN             ! Return immediately
       END IF

       total = total + partial

    END DO

    total%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential

  FUNCTION potential_1 ( i1, i2, ri, i, box, r_cut, j_range ) RESULT ( partial )
    IMPLICIT NONE
    TYPE(potential_type)            :: partial ! Returns a composite of pot, vir etc for given atom
    INTEGER,            INTENT(in)  :: i1, i2  ! Index range defining system: (1,n(1)) or (n(1)+1,n(1)+n(2))
    REAL, DIMENSION(3), INTENT(in)  :: ri      ! Coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i       ! Index of atom of interest
    REAL,               INTENT(in)  :: box     ! Simulation box length
    REAL,               INTENT(in)  :: r_cut   ! Potential cutoff distance
    INTEGER, OPTIONAL,  INTENT(in)  :: j_range ! Optional partner index range

    ! partial%pot is the nonbonded cut (not shifted) potential energy of atom ri with a set of other atoms
    ! partial%vir is the corresponding virial of atom ri
    ! partial%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the values of partial%pot etc should not be used
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER              :: j, j1, j2
    REAL                 :: r_cut_box, r_cut_box_sq, box_sq
    REAL                 :: sr2, sr6, sr12, rij_sq
    REAL, DIMENSION(3)   :: rij
    REAL, PARAMETER      :: sr2_ovr = 1.77 ! overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    IF ( SUM(n) > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', SUM(n), SIZE(r,dim=2)
       STOP 'Error in potential_1'
    END IF

    IF ( PRESENT ( j_range ) ) THEN
       SELECT CASE ( j_range )
       CASE ( lt ) ! j < i
          j1 = i1
          j2 = i-1
       CASE ( gt ) ! j > i
          j1 = i+1
          j2 = i2
       CASE default ! should never happen
          WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error ', j_range
          STOP 'Impossible error in potential_1'
       END SELECT
    ELSE
       j1 = i1
       j2 = i2
    END IF

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box**2
    box_sq       = box**2

    partial = potential_type ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( i == j ) CYCLE ! Skip self

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
          pair%pot = sr12 - sr6      ! LJ pair potential (cut but not shifted)
          pair%vir = pair%pot + sr12 ! LJ pair virial

          partial = partial + pair 

       END IF ! End check within range

    END DO ! End loop over selected range of partners

    ! Include numerical factors
    partial%pot = partial%pot * 4.0        ! 4*epsilon
    partial%vir = partial%vir * 24.0 / 3.0 ! 24*epsilon and divide virial by 3
    partial%ovr = .FALSE.                  ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential_1

  SUBROUTINE move ( i, ri )
    IMPLICIT NONE
    INTEGER,               INTENT(in) :: i
    REAL,    DIMENSION(3), INTENT(in) :: ri

    r(:,i) = ri ! New position

  END SUBROUTINE move

  SUBROUTINE swap ( i, ri )
    IMPLICIT NONE
    INTEGER,            INTENT(in) :: i  ! Index of particle to be destroyed
    REAL, DIMENSION(3), INTENT(in) :: ri ! Position of particle to be created

    INTEGER :: t

    IF ( i <= n(1) ) THEN ! Destroy in system 1, create in system 2
       t      = n(1)          ! Last atom in system 1
       r(:,i) = r(:,t)        ! Replace i coordinates with t coordinates
       r(:,t) = ri            ! New particle coordinates
       n(:)   = n(:) + [-1,1] ! Move boundary down to include t
    ELSE ! Destroy in system 2, create in system 1
       t      = n(1)+1        ! First atom in system 2
       r(:,i) = r(:,t)        ! Replace i coordinates with t coordinates
       r(:,t) = ri            ! New particle coordinates
       n(:)   = n(:) + [1,-1] ! Move boundary up to include t
    END IF

    IF ( MINVAL(n) <= 0 ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)') 'Number of particles is zero', n
       STOP 'Error in mc_gibbs_lj'
    END IF

  END SUBROUTINE swap

END MODULE mc_module
