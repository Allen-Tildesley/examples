! smc_lj_module.f90
! Energy, force, and move routines for SMC, LJ potential
MODULE smc_module
  
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, r, r_old, v, f, f1, move, lt, ne, gt
  PUBLIC :: model_description, allocate_arrays, deallocate_arrays
  public :: force, force1, energy_lrc

  INTEGER                              :: n     ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r     ! positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r_old ! old positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: v     ! velocities (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f     ! forces (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f1    ! forces (3,n) used in force1
  logical, dimension(:,:), allocatable :: move  ! mask for multi-atom moves

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options
  REAL,    PARAMETER :: sigma = 1.0 ! LJ diameter (unit of length)
  REAL,    PARAMETER :: epslj = 1.0 ! LJ well depth (unit of energy)

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

    ALLOCATE ( r(3,n), r_old(3,n), v(3,n), f(3,n), f1(3,n), move(3,n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, r_old, v, f, f1, move )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( box, r_cut, pot, pot_sh, vir )
    REAL, INTENT(in)  :: box    ! simulation box length
    REAL, INTENT(in)  :: r_cut  ! potential cutoff distance
    REAL, INTENT(out) :: pot    ! total potential energy
    REAL, INTENT(out) :: pot_sh ! potential shifted to zero at cutoff
    REAL, INTENT(out) :: vir    ! virial

    ! Calculates potential (unshifted and shifted), virial and forces
    ! Actual calculation is performed by subroutine energy_1
    ! Forces are accumulated in f; contributions calculated by force1 are in f1

    INTEGER            :: i
    REAL               :: pot1, pot1_sh, vir1

    f      = 0.0
    pot    = 0.0
    pot_sh = 0.0
    vir    = 0.0

    DO i = 1, n - 1 ! Begin loop over atoms

       call force1 ( box, r_cut, i, gt, pot1, pot1_sh, vir1 )
       f      = f + f1
       pot    = pot + pot1
       pot_sh = pot_sh + pot1_sh
       vir    = vir + vir1
      
    END DO ! End loop over atoms

  END SUBROUTINE force

  SUBROUTINE force1 ( box, r_cut, i, j_range, pot, pot_sh, vir )
    REAL,    INTENT(in)  :: box     ! simulation box length
    REAL,    INTENT(in)  :: r_cut   ! potential cutoff distance
    integer, intent(in)  :: i       ! index of atom of interest
    integer, intent(in)  :: j_range ! index range of partner atoms
    REAL,    INTENT(out) :: pot     ! total potential energy
    REAL,    INTENT(out) :: pot_sh  ! potential shifted to zero at cutoff
    REAL,    INTENT(out) :: vir     ! virial

    ! Calculates potential (unshifted and shifted), virial and forces f1
    ! It is assumed that positions are in units where box = 1
    ! Forces are calculated in units where sigma = 1 and epsilon = 1

    INTEGER            :: j, j1, j2, n_cut
    REAL               :: r_cut_box, r_cut_box_sq, box_sq, rij_sq, sr2, sr6, sr12, potij, virij
    REAL, DIMENSION(3) :: rij, fij

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    f1    = 0.0
    pot   = 0.0
    vir   = 0.0
    n_cut = 0

    SELECT CASE ( j_range )
    CASE ( lt ) ! j < i
       j1 = 1
       j2 = i-1
    CASE ( gt ) ! j > i
       j1 = i+1
       j2 = n
    CASE ( ne ) ! j /= i
       j1 = 1
       j2 = n
    END SELECT

       DO j = j1, j2 ! Begin loop over atoms

          if ( j == i ) cycle
          
          rij(:) = r(:,i) - r(:,j)           ! separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )            ! squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN ! check within cutoff

             rij_sq  = rij_sq * box_sq ! now in sigma=1 units
             rij(:)  = rij(:) * box    ! now in sigma=1 units
             sr2     = 1.0 / rij_sq
             sr6     = sr2 ** 3
             sr12    = sr6 ** 2
             potij   = sr12 - sr6
             virij   = potij + sr12
             pot     = pot + potij
             vir     = vir + virij
             fij     = rij * virij / rij_sq
             f1(:,i) = f1(:,i) + fij
             f1(:,j) = f1(:,j) - fij
             n_cut   = n_cut + 1

          END IF ! End check within cutoff

       END DO ! End inner loop over atoms

    ! calculate shifted potential
    sr2    = 1.0 / r_cut**2 ! in sigma=1 units
    sr6    = sr2 ** 3
    sr12   = sr6 **2
    potij  = sr12 - sr6
    pot_sh = pot - REAL ( n_cut ) * potij

    ! multiply results by numerical factors
    f1     = f1      * 24.0
    pot    = pot    * 4.0
    pot_sh = pot_sh * 4.0
    vir    = vir    * 24.0 / 3.0

  END SUBROUTINE force1


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

END MODULE smc_module
