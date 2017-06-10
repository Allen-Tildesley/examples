! md_lj_le_module.f90
! Force routine for MD, LJ atoms, Lees-Edwards boundaries
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
  PUBLIC :: force

  ! Public data
  INTEGER,                              PUBLIC :: n ! Number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v ! Velocities (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: f ! Forces (3,n)

  ! Private data
  REAL, PARAMETER :: r_cut_sq = 2.0**(1.0/3.0), r_cut = SQRT(r_cut_sq) ! Minimum and cutoff of WCA LJ potential

  ! Public derived type
  TYPE, PUBLIC :: potential_type ! A composite variable for interactions comprising
     REAL    :: pot ! the potential energy
     REAL    :: vir ! the virial and
     REAL    :: lap ! the Laplacian and
     REAL    :: pyx ! the off-diagonal virial part of the pressure tensor
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
    c%pot = a%pot  +   b%pot
    c%vir = a%vir  +   b%vir
    c%pyx = a%pyx  +   b%pyx
    c%lap = a%lap  +   b%lap
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)' ) 'WCA shifted Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'  
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Lees-Edwards boundaries'

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box )
    IMPLICIT NONE
    REAL, INTENT(in) :: box   ! Simulation box length

    REAL :: r_cut_box

    ALLOCATE ( r(3,n), v(3,n), f(3,n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)' ) 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, v, f )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( box, strain, total )
    IMPLICIT NONE
    REAL,                 INTENT(in)  :: box    ! Simulation box length
    REAL,                 INTENT(in)  :: strain ! Shear strain
    TYPE(potential_type), INTENT(out) :: total  ! Composite of pot, vir, lap etc

    ! total%pot is the WCA LJ potential energy for whole system
    ! total%vir is the corresponding virial
    ! total%pyx is the virial part of the yx-element of the pressure tensor
    ! total%lap is the corresponding Laplacian
    ! total%ovr is a warning flag that there is an overlap
    ! This routine also calculates forces and stores them in the array f
    ! If total%ovr is set to .true., the forces etc should not be used

    ! It is assumed that positions are in units where box = 1
    ! Forces are calculated in units where sigma = 1 and epsilon = 1
    ! Lees-Edwards boundaries, in sliding brick arrangement
    ! Flow/gradient/vorticity directions are x/y/z == 1/2/3

    INTEGER              :: i, j
    REAL                 :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL                 :: sr2, sr6, sr12
    REAL, DIMENSION(3)   :: rij, fij
    REAL, PARAMETER      :: sr2_ovr = 1.77 ! overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    ! Initialize
    f     = 0.0
    total = potential_type (  pot=0.0, vir=0.0, pyx=0.0, lap=0.0, ovr=.FALSE. )

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms

          rij(:) = r(:,i) - r(:,j)                    ! Separation vector
          rij(1) = rij(1) - ANINT ( rij(2) ) * strain ! Extra correction in box=1 units
          rij(:) = rij(:) - ANINT ( rij(:) )          ! Periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )                     ! Squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN ! Check within cutoff

             rij_sq   = rij_sq * box_sq ! Now in sigma=1 units
             rij(:)   = rij(:) * box    ! Now in sigma=1 units
             sr2      = 1.0 / rij_sq    ! (sigma/rij)**2
             pair%ovr = sr2 > sr2_ovr   ! Overlap if too close

             sr6      = sr2 ** 3
             sr12     = sr6 ** 2
             pair%pot = sr12 - sr6                    ! LJ pair potential (cut but not shifted)
             pair%vir = pair%pot + sr12               ! LJ pair virial
             pair%pot = pair%pot + 0.25               ! WCA LJ pair potential (cut-and-shifted)
             pair%lap = ( 22.0*sr12 - 5.0*sr6 ) * sr2 ! LJ pair Laplacian
             fij      = rij * pair%vir * sr2          ! LJ pair forces
             pair%pyx = rij(2)*fij(1)                 ! Off-diagonal element of pressure tensor

             total  = total  + pair
             f(:,i) = f(:,i) + fij
             f(:,j) = f(:,j) - fij

          END IF ! End check within cutoff

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    ! Multiply results by numerical factors
    f         = f         * 24.0       ! 24*epsilon
    total%pot = total%pot * 4.0        ! 4*epsilon
    total%vir = total%vir * 24.0 / 3.0 ! 24*epsilon and divide virial by 3
    total%pyx = total%pyx * 24.0       ! 24*epsilon
    total%lap = total%lap * 24.0 * 2.0 ! 24*epsilon and factor 2 for ij and ji

  END SUBROUTINE force

END MODULE md_module

