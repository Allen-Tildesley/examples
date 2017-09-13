! qmc_pi_lj_module.f90
! Energy and move routines for PIMC simulation, LJ potential
MODULE qmc_module

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
  PUBLIC :: potential_1, potential, spring_1, spring

  ! Public data
  INTEGER,                                PUBLIC :: n ! Number of atoms
  INTEGER,                                PUBLIC :: p ! Number of beads
  REAL,    DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n,p)

  ! Private data
  INTEGER, PARAMETER :: lt = -1, gt = 1 ! j-range and l-range options

  ! Public derived type
  TYPE, PUBLIC :: potential_type ! A composite variable for interaction energies comprising
     REAL    :: pot ! the cut (but not shifted) potential energy and
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

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
    IMPLICIT NONE
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

    REAL :: r_cut_box

    ALLOCATE ( r(3,n,p) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)') 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r )

  END SUBROUTINE deallocate_arrays

  FUNCTION potential ( box, r_cut ) RESULT ( total )
    IMPLICIT NONE
    TYPE(potential_type) :: total ! Returns a composite of pot, vir and ovr
    REAL, INTENT(in)     :: box   ! Simulation box length
    REAL, INTENT(in)     :: r_cut ! Potential cutoff distance

    ! total%pot is the nonbonded cut (not shifted) classical potential energy for whole system
    ! total%vir is the corresponding virial for whole system
    ! total%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of total%pot etc should not be used
    ! Actual calculation is performed by function potential_1

    TYPE(potential_type) :: partial ! Atomic contribution to total
    INTEGER              :: i, k

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in potential'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in potential'
    END IF

    total = potential_type ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO k = 1, p ! Loop over ring polymers
       DO i = 1, n - 1 ! Loop over atoms within polymer

          partial = potential_1 ( r(:,i,k), i, k, box, r_cut, gt )

          IF ( partial%ovr ) THEN
             total%ovr = .TRUE. ! Overlap detected
             RETURN             ! Return immediately
          END IF

          total = total + partial

       END DO ! End loop over atoms within polymer
    END DO ! End loop over ring polymers

    total%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential

  FUNCTION potential_1 ( rik, i, k, box, r_cut, j_range ) RESULT ( partial )
    IMPLICIT NONE
    TYPE(potential_type)           :: partial ! Returns a composite of pot, vir and ovr for given atom
    REAL, DIMENSION(3), INTENT(in) :: rik     ! Coordinates of atom of interest
    INTEGER,            INTENT(in) :: i, k    ! Index and polymer id of atom of interest
    REAL,               INTENT(in) :: box     ! Simulation box length
    REAL,               INTENT(in) :: r_cut   ! Potential cutoff distance
    INTEGER, OPTIONAL,  INTENT(in) :: j_range ! Optional partner index range

    ! partial%pot is the nonbonded cut (not shifted) classical potential energy of atom rik with a set of other atoms
    ! partial%vir is the corresponding virial of atom rik
    ! partial%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the value of partial%pot should not be used
    ! The coordinates in rik are not necessarily identical with those in r(:,i,k)
    ! The partner atoms always have the same polymer index k
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER              :: j, j1, j2
    REAL                 :: r_cut_box, r_cut_box_sq, box_sq
    REAL                 :: sr2, sr6, sr12, r_ik_jk_sq
    REAL, DIMENSION(3)   :: r_ik_jk
    REAL, PARAMETER      :: sr2_ovr = 1.77 ! overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in potential_1'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in potential_1'
    END IF

    IF ( PRESENT ( j_range ) ) THEN
       SELECT CASE ( j_range )
       CASE ( lt ) ! j < i
          j1 = 1
          j2 = i-1
       CASE ( gt ) ! j > i
          j1 = i+1
          j2 = n
       CASE default ! should never happen
          WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error ', j_range
          STOP 'Impossible error in potential_1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box**2
    box_sq       = box**2

    partial = potential_type ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( i == j ) CYCLE ! Skip self

       r_ik_jk(:) = rik(:) - r(:,j,k)                 ! Separation vector
       r_ik_jk(:) = r_ik_jk(:) - ANINT ( r_ik_jk(:) ) ! Periodic boundaries in box=1 units
       r_ik_jk_sq = SUM ( r_ik_jk**2 )                ! Squared separation in box=1 units

       IF ( r_ik_jk_sq < r_cut_box_sq ) THEN ! Check within range

          r_ik_jk_sq = r_ik_jk_sq * box_sq ! now in sigma=1 units
          sr2        = 1.0 / r_ik_jk_sq    ! (sigma/rikjk)**2
          pair%ovr   = sr2 > sr2_ovr       ! Overlap if too close

          IF ( pair%ovr ) THEN
             partial%ovr = .TRUE. ! Overlap detected
             RETURN               ! Return immediately
          END IF

          sr6      = sr2**3
          sr12     = sr6**2
          pair%pot = sr12 - sr6      ! LJ potential (cut but not shifted)
          pair%vir = pair%pot + sr12 ! LJ pair virial

          partial = partial + pair

       END IF ! End check with range

    END DO ! End loop over selected range of partners

    ! Include numerical factors
    partial%pot = partial%pot * 4.0 / REAL ( p ) ! Classical potentials are weaker by a factor p
    partial%vir = partial%vir * 24.0 / REAL(3*p) ! 24*epsilon and divide virial by 3 & by p
    partial%ovr = .FALSE.                        ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential_1

  FUNCTION spring ( box, k_spring ) RESULT ( total )
    IMPLICIT NONE
    REAL             :: total    ! Returns quantum spring potential 
    REAL, INTENT(in) :: box      ! Simulation box length
    REAL, INTENT(in) :: k_spring ! Spring potential constant

    ! total is the quantum spring potential for whole system
    ! Actual calculation is performed by subroutine spring_1

    INTEGER :: i, k, kp
    REAL    :: partial_pot

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in spring'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in spring'
    END IF

    total = 0.0

    DO k = 1, p ! Loop over atoms within polymer

       kp = k+1
       IF ( kp > p ) kp = 1
       
       DO i = 1, n ! Loop over ring polymers

          partial_pot = spring_1 ( r(:,i,k), r(:,i,kp), box, k_spring )
          total       = total + partial_pot

       END DO ! End loop over ring polymers

    END DO ! End loop over atoms within polymer

  END FUNCTION spring

  FUNCTION spring_1 ( rik, ril, box, k_spring ) RESULT ( partial_pot )
    IMPLICIT NONE
    REAL                           :: partial_pot ! Returns quantum potential for given atom
    REAL, DIMENSION(3), INTENT(in) :: rik         ! Coordinates of atom of interest
    REAL, DIMENSION(3), INTENT(in) :: ril         ! Coordinates of other atom of interest
    REAL,               INTENT(in) :: box         ! Simulation box length
    REAL,               INTENT(in) :: k_spring    ! Spring potential constant

    ! partial_pot is the quantum spring energy of atom (i,k) in rik for polymer bead k
    ! with its neighbour (i,l) in the same polymer ring
    ! The coordinates in rik are not necessarily identical with those in r(:,i,k)
    ! The partner atom index is always the same as i

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    REAL               :: r_ik_il_sq
    REAL, DIMENSION(3) :: r_ik_il

    r_ik_il(:)  = rik(:) - ril(:)
    r_ik_il(:)  = r_ik_il(:) - ANINT ( r_ik_il(:) ) ! Periodic boundaries in box=1 units
    r_ik_il_sq  = SUM ( r_ik_il**2 ) *  box**2      ! Squared distance in sigma=1 units
    partial_pot = 0.5 * k_spring * r_ik_il_sq       ! Spring potential

  END FUNCTION spring_1

END MODULE qmc_module
