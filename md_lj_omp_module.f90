! md_lj_omp_module.f90
! Force routine for MD simulation, Lennard-Jones atoms, OpenMP
MODULE md_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
  USE, INTRINSIC :: omp_lib

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

  REAL, SAVE :: wall_time

CONTAINS

  SUBROUTINE introduction ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut-and-shifted version for dynamics'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut (but not shifted) version also calculated'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'    
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'

    WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Max # of OpenMP threads = ', omp_get_max_threads()
    wall_time = omp_get_wtime()

  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'OpenMP walltime = ', omp_get_wtime()-wall_time
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
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, v, f )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( box, r_cut, pot, cut, vir, lap, overlap )
    IMPLICIT NONE
    REAL,    INTENT(in)  :: box     ! Simulation box length
    REAL,    INTENT(in)  :: r_cut   ! Potential cutoff distance
    REAL,    INTENT(out) :: pot     ! Cut-and-shifted total potential energy
    REAL,    INTENT(out) :: cut     ! Cut (but not shifted) total potential energy
    REAL,    INTENT(out) :: vir     ! Total virial
    REAL,    INTENT(out) :: lap     ! Total Laplacian
    LOGICAL, INTENT(out) :: overlap ! Warning flag that there is an overlap

    ! Calculates forces in array f, and also pot, cut etc  
    ! If overlap is set to .true., the forces etc should not be used
    ! It is assumed that positions are in units where box = 1
    ! Forces are calculated in units where sigma = 1 and epsilon = 1

    INTEGER            :: i, j, ncut
    REAL               :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL               :: sr2, sr6, sr12, cutij, virij, lapij
    REAL, DIMENSION(3) :: rij, fij
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! Overlap threshold

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    ! Initialize (redundant: all these, except pot, are initialized in reduction clause)
    f       = 0.0
    pot     = 0.0
    cut     = 0.0
    vir     = 0.0
    lap     = 0.0
    overlap = .FALSE.
    ncut    = 0

    ! Each thread must use its own private copies of most of the variables
    ! The shared variables, n, r_cut_box_sq, box, box_sq, are never changed
    ! The reduction variables are private, separately accumulated in each thread,
    ! and then added together to give global totals at the end
    ! The schedule(static,1) attempts to share the i-loop iterations between threads
    ! fairly evenly, given the i-dependence of the work in the inner j-loop
    ! Feel free to experiment with this

    !$omp parallel do &
    !$omp& default(shared), private(i,j,rij,rij_sq,sr2,sr6,sr12,cutij,virij,lapij,fij) &
    !$omp& reduction(+:cut,vir,lap,ncut,f), reduction(.or.:overlap) &
    !$omp& schedule(static,1)
    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms

          rij(:) = r(:,i) - r(:,j)           ! Separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )            ! Squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN

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

          END IF

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    ! Calculate shifted potential
    sr2   = 1.0 / r_cut**2 ! in sigma=1 units
    sr6   = sr2**3
    sr12  = sr6**2
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
    REAL,    INTENT(in) :: density    ! Number density N/V
    REAL,    INTENT(in) :: r_cut      ! Cutoff distance

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

  FUNCTION hessian ( box, r_cut ) RESULT ( hes )
    IMPLICIT NONE
    REAL             :: hes   ! Returns the total Hessian
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

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

    hes = 0.0

    ! Similar comments apply here as for the double loop in the force routine

    !$omp parallel do &
    !$omp& default(shared), private(i,j,rij,fij,rij_sq,ff,rf,sr2,sr6,sr8,sr10,v1,v2) &
    !$omp& reduction(+:hes) &
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
             hes  = hes + v1 * ff + v2 * rf**2

          END IF ! End check within cutoff

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

  END FUNCTION hessian

END MODULE md_module
