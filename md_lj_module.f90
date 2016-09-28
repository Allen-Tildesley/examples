! md_lj_module.f90
! Force routine for MD simulation, Lennard-Jones atoms
MODULE md_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, r, v, f
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: force, hessian, energy_lrc

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: v ! velocities (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f ! forces (3,n)

  REAL, PARAMETER :: sigma = 1.0 ! LJ diameter (unit of length)
  REAL, PARAMETER :: epslj = 1.0 ! LJ well depth (unit of energy)

CONTAINS

  SUBROUTINE introduction ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ',     sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epsilon = ', epslj    
  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output
    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'
  END SUBROUTINE conclusion

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

  SUBROUTINE force ( box, r_cut, pot, pot_sh, vir, lap )
    REAL, INTENT(in)  :: box    ! simulation box length
    REAL, INTENT(in)  :: r_cut  ! potential cutoff distance
    REAL, INTENT(out) :: pot    ! total potential energy
    REAL, INTENT(out) :: pot_sh ! potential shifted to zero at cutoff
    REAL, INTENT(out) :: vir    ! virial
    REAL, INTENT(out) :: lap    ! Laplacian

    ! Calculates potential (unshifted and shifted), virial, Laplacian, and forces
    ! It is assumed that positions are in units where box = 1
    ! Forces are calculated in units where sigma = 1 and epsilon = 1

    INTEGER            :: i, j, n_cut
    REAL               :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL               :: sr2, sr6, sr12, potij, virij, lapij
    REAL, DIMENSION(3) :: rij, fij

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    f     = 0.0
    pot   = 0.0
    vir   = 0.0
    lap   = 0.0
    n_cut = 0

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms

          rij(:) = r(:,i) - r(:,j)           ! separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )            ! squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN ! check within cutoff

             rij_sq = rij_sq * box_sq ! now in sigma=1 units
             rij(:) = rij(:) * box    ! now in sigma=1 units
             sr2    = 1.0 / rij_sq
             sr6    = sr2 ** 3
             sr12   = sr6 ** 2
             potij  = sr12 - sr6
             virij  = potij + sr12
             lapij  = ( 22.0*sr12 - 5.0*sr6 ) * sr2
             pot    = pot + potij
             vir    = vir + virij
             lap    = lap + lapij
             fij    = rij * virij / rij_sq
             f(:,i) = f(:,i) + fij
             f(:,j) = f(:,j) - fij
             n_cut  = n_cut + 1

          END IF ! End check within cutoff

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    ! calculate shifted potential
    sr2    = 1.0 / r_cut**2 ! in sigma=1 units
    sr6    = sr2 ** 3
    sr12   = sr6 **2
    potij  = sr12 - sr6
    pot_sh = pot - REAL ( n_cut ) * potij

    ! multiply results by numerical factors
    f      = f      * 24.0
    pot    = pot    * 4.0
    pot_sh = pot_sh * 4.0
    vir    = vir    * 24.0 / 3.0
    lap    = lap    * 24.0 * 2.0

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
    density = REAL(n)/box**3
    pot     = pot * pi * density * REAL(n)
    vir     = vir * pi * density * REAL(n)

  END SUBROUTINE energy_lrc

  FUNCTION hessian ( box, r_cut ) RESULT ( hes )
    REAL, INTENT(in)  :: box    ! simulation box length
    REAL, INTENT(in)  :: r_cut  ! potential cutoff distance
    REAL              :: hes    ! result

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

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms
          rij(:) = r(:,i) - r(:,j)           ! separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )            ! squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN ! check within cutoff

             rij_sq = rij_sq * box_sq ! now in sigma=1 units
             rij(:) = rij(:) * box    ! now in sigma=1 units
             fij(:) = f(:,i) - f(:,j) ! difference in forces

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
