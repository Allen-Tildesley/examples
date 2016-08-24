! md_lj_omp_module.f90
! Force routine for MD simulation, Lennard-Jones atoms, OpenMP
MODULE md_module

  ! TODO MPA complete this code
  
    USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, v, f
  PUBLIC :: model_description, allocate_arrays, deallocate_arrays
  public :: force, energy_lrc

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,:)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: v ! velocities (3,:)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f ! forces (3,:)

  REAL, PARAMETER :: sigma = 1.0 ! LJ diameter (unit of length)
  REAL, PARAMETER :: epslj = 1.0 ! LJ well depth (unit of energy)

CONTAINS

  SUBROUTINE model_description ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ',     sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epsilon = ', epslj    
  END SUBROUTINE model_description

  SUBROUTINE allocate_arrays ( box, r_cut )
    REAL, INTENT(in) :: box   ! simulation box length
    REAL, INTENT(in) :: r_cut ! potential cutoff distance

    REAL :: r_cut_box

    ALLOCATE ( r(3,n), v(3,n), f(3,n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, v, f )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( box, r_cut, pot, pot_sh, vir )
    REAL, INTENT(in)  :: box    ! simulation box length
    REAL, INTENT(in)  :: r_cut  ! potential cutoff distance
    REAL, INTENT(out) :: pot    ! total potential energy
    REAL, INTENT(out) :: pot_sh ! potential shifted to zero at cutoff
    REAL, INTENT(out) :: vir    ! virial

    ! Calculates potential (unshifted and shifted), virial and forces
    ! It is assumed that positions are in units where box = 1

    INTEGER            :: i, j, n_cut
    REAL               :: r_cut_box, r_cut_box_sq, box_sq, rij_sq, sr2, sr6, sr12, potij, virij
    REAL, DIMENSION(3) :: rij, fij

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    f     = 0.0
    pot   = 0.0
    vir   = 0.0
    n_cut = 0

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms

          rij(:) = r(:,i) - r(:,j)           ! separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )            ! squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN

             rij_sq = rij_sq * box_sq ! now in sigma=1 units
             rij(:) = rij(:) * box    ! now in sigma=1 units
             sr2    = 1.0 / rij_sq
             sr6    = sr2 ** 3
             sr12   = sr6 ** 2
             potij  = sr12 - sr6
             virij  = potij + sr12
             pot    = pot + potij
             vir    = vir + virij
             fij    = rij * virij / rij_sq
             f(:,i) = f(:,i) + fij
             f(:,j) = f(:,j) - fij
             n_cut  = n_cut + 1

          END IF

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    ! calculate shifted potential
    sr2    = 1.0 / r_cut**2 ! in sigma=1 units
    sr6    = sr2**3
    sr12   = sr6**2
    potij  = sr12 - sr6
    pot_sh = pot - REAL ( n_cut ) * potij

    ! multiply results by numerical factors

    f      = f      * 24.0
    pot    = pot    * 4.0
    pot_sh = pot_sh * 4.0
    vir    = vir    * 24.0 / 3.0

  END SUBROUTINE force

  SUBROUTINE energy_lrc ( n, box, r_cut, pot, vir )
    INTEGER, INTENT(in)  :: n        ! number of atoms
    REAL,    INTENT(in)  :: box      ! simulation box length
    REAL,    INTENT(in)  :: r_cut    ! potential cutoff distance
    REAL,    INTENT(out) :: pot, vir ! potential and virial

    ! Calculates long-range corrections for Lennard-Jones potential and virial
    ! These are the corrections to the total values
    ! Results are in LJ units where sigma = 1, epsilon = 1

    REAL               :: sr3, density
    REAL, PARAMETER    :: pi = 4.0 * ATAN(1.0)

    sr3     = 1.0 / r_cut**3
    pot     = (8.0/9.0)  * sr3**3 - (8.0/3.0)  * sr3
    vir     = (32.0/9.0) * sr3**3 - (32.0/6.0) * sr3
    density =  REAL(n)/box**3
    pot     = pot * pi * density * REAL(n)
    vir     = vir * pi * density * REAL(n)

  END SUBROUTINE energy_lrc

END MODULE md_module
