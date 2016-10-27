! md_lj_le_module.f90
! Force routine for MD, LJ atoms, Lees-Edwards boundaries
MODULE md_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: force, potential_lrc, pressure_lrc

  ! Public data
  INTEGER,                              PUBLIC :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v ! velocities (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: f ! forces (3,n)

CONTAINS

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

  SUBROUTINE force ( box, r_cut, strain, pot, cut, vir, lap, overlap )
    IMPLICIT NONE
    REAL,    INTENT(in)  :: box     ! Simulation box length
    REAL,    INTENT(in)  :: r_cut   ! Potential cutoff distance
    REAL,    INTENT(in)  :: strain  ! Shear strain
    REAL,    INTENT(out) :: pot     ! Cut-and-shifted total potential energy
    REAL,    INTENT(out) :: cut     ! Cut (but not shifted) total potential energy
    REAL,    INTENT(out) :: vir     ! Total virial
    REAL,    INTENT(out) :: lap     ! Total Laplacian
    LOGICAL, INTENT(out) :: overlap ! Warning flag that there is an overlap

    ! Calculates forces in array f, and also pot, cut etc
    ! If overlap is set to .true., the forces etc should not be used
    ! It is assumed that positions are in units where box = 1
    ! Forces are calculated in units where sigma = 1 and epsilon = 1
    ! Lees-Edwards boundaries, in sliding brick arrangement
    ! Flow/gradient/vorticity directions are x/y/z == 1/2/3

    INTEGER            :: i, j, ncut
    REAL               :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL               :: sr2, sr6, sr12, cutij, virij, lapij
    REAL, DIMENSION(3) :: rij, fij
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    ! Initialize
    f       = 0.0
    pot     = 0.0
    cut     = 0.0
    vir     = 0.0
    lap     = 0.0
    overlap = .FALSE.
    ncut    = 0

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms

          rij(:) = r(:,i) - r(:,j)                    ! Separation vector
          rij(1) = rij(1) - ANINT ( rij(2) ) * strain ! Extra correction in box=1 units
          rij(:) = rij(:) - ANINT ( rij(:) )          ! Periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )                     ! Squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN ! Check within cutoff

             rij_sq = rij_sq * box_sq ! Now in sigma=1 units
             rij(:) = rij(:) * box    ! Now in sigma=1 units
             sr2    = 1.0 / rij_sq

             IF ( sr2 > sr2_overlap ) overlap = .TRUE. ! Overlap detected

             sr6   = sr2 ** 3
             sr12  = sr6 ** 2
             cutij = sr12 - sr6                    ! LJ pair potential (cut but not shifted)
             virij = cutij + sr12                  ! LJ pair virial
             lapij = ( 22.0*sr12 - 5.0*sr6 ) * sr2 ! LJ pair Laplacian
             fij   = rij * virij / rij_sq          ! LJ pair forces

             cut    = cut    + cutij
             vir    = vir    + virij
             lap    = lap    + lapij
             f(:,i) = f(:,i) + fij
             f(:,j) = f(:,j) - fij
             ncut   = ncut   + 1

          END IF ! end check within cutoff

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    ! Calculate cut&shifted potential
    sr2   = 1.0 / r_cut**2 ! in sigma=1 units
    sr6   = sr2 ** 3
    sr12  = sr6 **2
    cutij = sr12 - sr6
    pot   = cut - REAL ( ncut ) * cutij

    ! Multiply results by numerical factors
    f   = f   * 24.0
    cut = cut * 4.0
    pot = pot * 4.0
    vir = vir * 24.0 / 3.0
    lap = lap * 24.0 * 2.0

  END SUBROUTINE force

  FUNCTION potential_lrc ( density, r_cut )
    IMPLICIT NONE
    REAL                :: potential_lrc ! Returns long-range energy/atom
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
    REAL                :: pressure_lrc ! Returns long-range pressure
    REAL,    INTENT(in) :: density      ! Number density N/V
    REAL,    INTENT(in) :: r_cut        ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones pressure
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3          = 1.0 / r_cut**3
    pressure_lrc = pi * ( (32.0/9.0) * sr3**3  - (16.0/3.0) * sr3 ) * density**2

  END FUNCTION pressure_lrc

END MODULE md_module

