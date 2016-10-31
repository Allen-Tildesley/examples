! md_lj_module.f90
! Force routine for MD simulation, Lennard-Jones atoms
MODULE md_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: force, hessian, potential_lrc, pressure_lrc

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

  SUBROUTINE introduction ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut-and-shifted version for dynamics'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut (but not shifted) version also calculated'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'    
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'  

  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
    IMPLICIT NONE
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

    REAL :: r_cut_box

    ALLOCATE ( r(3,n), v(3,n), f(3,n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'r_cut/box too large ', r_cut_box
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

    INTEGER              :: i, j
    REAL                 :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL                 :: sr2, sr6, sr12, pot_cut
    REAL, DIMENSION(3)   :: rij, fij
    REAL, PARAMETER      :: sr2_ovr = 1.77 ! Overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    ! Calculate potential at cutoff
    sr2     = 1.0 / r_cut**2 ! in sigma=1 units
    sr6     = sr2 ** 3
    sr12    = sr6 **2
    pot_cut = sr12 - sr6 ! Without numerical factor 4

    ! Initialize
    f     = 0.0
    total = potential_type ( pot=0.0, cut=0.0, vir=0.0, lap=0.0, ovr=.FALSE. )

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms

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

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    ! Multiply results by numerical factors
    f         = f         * 24.0       ! 24*epsilon
    total%cut = total%cut * 4.0        ! 4*epsilon
    total%pot = total%pot * 4.0        ! 4*epsilon
    total%vir = total%vir * 24.0 / 3.0 ! 24*epsilon and divide virial by 3
    total%lap = total%lap * 24.0 * 2.0 ! 24*epsilon and factor 2 for ij and ji

  END SUBROUTINE force

  FUNCTION potential_lrc ( density, r_cut )
    IMPLICIT NONE
    REAL                :: potential_lrc ! Returns long-range energy/atom correction
    REAL,    INTENT(in) :: density       ! Number density N/V
    REAL,    INTENT(in) :: r_cut         ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones energy per atom
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3           = 1.0 / r_cut**3
    potential_lrc = pi * ( (8.0/9.0)  * sr3**3  - (8.0/3.0)  * sr3 ) * density

  END FUNCTION potential_lrc

  FUNCTION pressure_lrc ( density, r_cut )
    IMPLICIT NONE
    REAL                :: pressure_lrc ! Returns long-range pressure correction
    REAL,    INTENT(in) :: density      ! Number density N/V
    REAL,    INTENT(in) :: r_cut        ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones pressure
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3          = 1.0 / r_cut**3
    pressure_lrc = pi * ( (32.0/9.0) * sr3**3  - (16.0/3.0) * sr3 ) * density**2

  END FUNCTION pressure_lrc

  FUNCTION hessian ( box, r_cut ) RESULT ( hes )
    IMPLICIT NONE
    REAL             :: hes   ! Returns the total Hessian
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

    ! Calculates Hessian function (for 1/N correction to config temp)
    ! This routine is only needed in a constant-energy ensemble
    ! It is assumed that positions are in units where box = 1
    ! but the result is given in units where sigma = 1 and epsilon = 1
    ! It is assumed that forces have already been calculated in array f

    INTEGER            :: i, j
    REAL               :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL               :: sr2, sr6, sr8, sr10, rf, ff, v1, v2
    REAL, DIMENSION(3) :: rij, fij

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    hes = 0.0

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
             hes  = hes + v1 * ff + v2 * rf**2

          END IF ! End check within cutoff

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

  END FUNCTION hessian

END MODULE md_module
