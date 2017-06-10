! md_lj_omp_module.f90
! Force routine for MD simulation, Lennard-Jones atoms, OpenMP
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
  USE, INTRINSIC :: omp_lib

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

  ! Private data
  REAL, SAVE :: wall_time

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

    WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Max # of OpenMP threads = ', omp_get_max_threads()
    wall_time = omp_get_wtime()

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'OpenMP walltime = ', omp_get_wtime()-wall_time
    WRITE ( unit=output_unit, fmt='(a)'           ) 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
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

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, v, f )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( box, r_cut, total )
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

    ! It is assumed that positions are in units where box = 1
    ! Forces are calculated in units where sigma = 1 and epsilon = 1

    INTEGER            :: i, j
    REAL               :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL               :: sr2, sr6, sr12, pot_cut
    REAL               :: pair_cut, total_cut
    REAL               :: pair_pot, total_pot
    REAL               :: pair_vir, total_vir
    REAL               :: pair_lap, total_lap
    LOGICAL            :: pair_ovr, total_ovr
    REAL, DIMENSION(3) :: rij, fij
    REAL, PARAMETER    :: sr2_ovr = 1.77 ! Overlap threshold (pot > 100)

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    ! Calculate potential at cutoff
    sr2     = 1.0 / r_cut**2 ! in sigma=1 units
    sr6     = sr2 ** 3
    sr12    = sr6 **2
    pot_cut = sr12 - sr6 ! Without numerical factor 4

    ! Initialize (redundant, since these variables should be initialized in the reduction clause)
    total_cut = 0.0
    total_pot = 0.0
    total_vir = 0.0
    total_lap = 0.0
    total_ovr = .false.
    f         = 0.0
    
    ! Each thread must use its own private copies of most of the variables
    ! The shared variables, n, r, r_cut_box_sq, box, box_sq, pot_cut are never changed
    ! The reduction variables are private, separately accumulated in each thread,
    ! and then added together to give global totals at the end
    ! They are automatically initialized to zero at the start of the loop
    ! The schedule(static,1) attempts to share the i-loop iterations between threads
    ! fairly evenly, given the i-dependence of the work in the inner j-loop
    ! Feel free to experiment with this

    !$omp parallel do &
    !$omp& default(none), private(i,j,rij,rij_sq,sr2,sr6,sr12,pair_cut,pair_pot,pair_vir,pair_lap,pair_ovr,fij), &
    !$omp& shared(n,r,r_cut_box_sq,box,box_sq,pot_cut), &
    !$omp& reduction(+:total_cut,total_pot,total_vir,total_lap,f), reduction(.or.:total_ovr), &
    !$omp& schedule(static,1)

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms

          rij(:) = r(:,i) - r(:,j)           ! Separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )            ! Squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN

             rij_sq   = rij_sq * box_sq ! Now in sigma=1 units
             rij(:)   = rij(:) * box    ! Now in sigma=1 units
             sr2      = 1.0 / rij_sq    ! (sigma/rij)**2
             pair_ovr = sr2 > sr2_ovr   ! Overlap if too close

             sr6      = sr2 ** 3
             sr12     = sr6 ** 2
             pair_cut = sr12 - sr6                    ! LJ pair potential (cut but not shifted)
             pair_vir = pair_cut + sr12               ! LJ pair virial
             pair_pot = pair_cut - pot_cut            ! LJ pair potential (cut-and-shifted)
             pair_lap = ( 22.0*sr12 - 5.0*sr6 ) * sr2 ! LJ pair Laplacian
             fij      = rij * pair_vir * sr2          ! LJ pair forces

             total_cut = total_cut  +   pair_cut
             total_pot = total_pot  +   pair_pot
             total_vir = total_vir  +   pair_vir
             total_lap = total_lap  +   pair_lap
             total_ovr = total_ovr .or. pair_ovr

             f(:,i) = f(:,i) + fij
             f(:,j) = f(:,j) - fij

          END IF

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    ! Multiply results by numerical factors
    f         = f         * 24.0       ! 24*epsilon
    total_cut = total_cut * 4.0        ! 4*epsilon
    total_pot = total_pot * 4.0        ! 4*epsilon
    total_vir = total_vir * 24.0 / 3.0 ! 24*epsilon and divide virial by 3
    total_lap = total_lap * 24.0 * 2.0 ! 24*epsilon and factor 2 for ij and ji

    total = potential_type ( cut=total_cut, pot=total_pot, vir=total_vir, lap=total_lap, ovr=total_ovr )

  END SUBROUTINE force

  FUNCTION hessian ( box, r_cut ) RESULT ( total_hes )
    IMPLICIT NONE
    REAL             :: total_hes ! Returns the total Hessian
    REAL, INTENT(in) :: box       ! Simulation box length
    REAL, INTENT(in) :: r_cut     ! Potential cutoff distance

    ! Calculates Hessian function (for 1/N correction to config temp)
    ! This routine is only needed in a constant-energy ensemble
    ! It is assumed that positions are in units where box = 1
    ! but the result is given in units where sigma = 1 and epsilon = 1
    ! It is assumed that forces have already been calculated 

    INTEGER            :: i, j
    REAL               :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL               :: sr2, sr6, sr8, sr10, rf, ff, v1, v2
    REAL, DIMENSION(3) :: rij, fij

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    total_hes = 0.0

    ! Similar comments apply here as for the double loop in the force routine

    !$omp parallel do &
    !$omp& default(none), private(i,j,rij,fij,rij_sq,ff,rf,sr2,sr6,sr8,sr10,v1,v2), &
    !$omp& shared(n,r,f,r_cut_box_sq,box_sq,box), &
    !$omp& reduction(+:total_hes), &
    !$omp& schedule(static,1)
    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms

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

             total_hes = total_hes + v1 * ff + v2 * rf**2

          END IF ! End check within cutoff

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

  END FUNCTION hessian

END MODULE md_module
