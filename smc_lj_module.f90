! smc_lj_module.f90
! Energy, force, and move routines for SMC, LJ potential
MODULE smc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, r, r_old, v, f, move
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: force, force_1
  public :: potovr

  INTEGER                              :: n     ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r     ! positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r_old ! old positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: v     ! velocities (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f     ! forces (3,n)
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: move  ! mask for multi-atom moves

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range
  REAL,    PARAMETER :: sigma = 1.0     ! LJ diameter (unit of length)
  REAL,    PARAMETER :: epslj = 1.0     ! LJ well depth (unit of energy)

  ! The use of n_max in this module is clumsy. We do it this way because at the time of writing
  ! gfortran does not implement parameterized derived types (part of the Fortran 2003 standard)
  INTEGER, PARAMETER :: n_max = 256 
  TYPE potovr ! A composite variable for interactions comprising
     REAL                     :: pot ! the potential energy and
     REAL                     :: vir ! the virial and
     REAL, DIMENSION(3,n_max) :: f   ! the forces and
     LOGICAL                  :: ovr ! a flag indicating overlap (i.e. pot too high to use)
  END TYPE potovr

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_potovr
  END INTERFACE OPERATOR (+)

CONTAINS

  FUNCTION add_potovr ( a, b ) RESULT (c)
    TYPE(potovr)             :: c    ! Result is the sum of the two inputs
    TYPE(potovr), INTENT(in) :: a, b
    c%pot = a%pot +    b%pot
    c%vir = a%vir +    b%vir
    c%f   = a%f   +    b%f
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potovr

  SUBROUTINE introduction ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Shifted Lennard-Jones potential'
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

    IF ( n > n_max ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i5)' ) 'n too large ', n, n_max
       STOP 'Error in allocate_arrays'
    END IF

    ALLOCATE ( r(3,n), r_old(3,n), v(3,n), f(3,n), move(3,n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, r_old, v, f, move )
  END SUBROUTINE deallocate_arrays

  FUNCTION force ( box, r_cut )
    TYPE(potovr)     :: force  ! Returns composite of forces, pot, vir and ovr
    REAL, INTENT(in) :: box    ! simulation box length
    REAL, INTENT(in) :: r_cut  ! potential cutoff distance

    ! Calculates potential (unshifted and shifted), virial and forces
    ! Actual calculation is performed by function force_1

    INTEGER            :: i
    TYPE(potovr) :: force_i

    force = potovr ( f=0.0, pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO i = 1, n - 1 ! Begin loop over atoms

       force_i = force_1 ( box, r_cut, i, gt )
       IF ( force_i%ovr ) THEN
          force%ovr = .TRUE. ! Overlap detected
          RETURN             ! Return immediately
       END IF
       force = force + force_i

    END DO ! End loop over atoms

    force%ovr = .FALSE. ! No overlaps detected (redundant but for clarity)

  END FUNCTION force

  FUNCTION force_1 ( box, r_cut, i, j_range ) RESULT ( force )
    TYPE(potovr)                  :: force   ! Returns composite of forces, pot, vir and ovr
    REAL,              INTENT(in) :: box     ! simulation box length
    REAL,              INTENT(in) :: r_cut   ! potential cutoff distance
    INTEGER,           INTENT(in) :: i       ! index of atom of interest
    INTEGER, OPTIONAL, INTENT(in) :: j_range ! Optional partner index range

    ! Calculates potential, virial and forces f1
    ! It is assumed that positions are in units where box = 1
    ! Forces are calculated in units where sigma = 1 and epsilon = 1
    ! Note that we use a shifted LJ potential here

    INTEGER            :: j, j1, j2
    REAL               :: r_cut_box, r_cut_box_sq, box_sq
    REAL               :: rij_sq, sr2, sr6, sr12, potij, virij, pot_cut
    REAL, DIMENSION(3) :: rij, fij
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

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
          STOP 'Impossible error in force1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    ! calculate shifted potential
    sr2     = 1.0 / r_cut**2 ! in sigma=1 units
    sr6     = sr2 ** 3
    sr12    = sr6 **2
    pot_cut = sr12 - sr6 ! without numerical factors

    force = potovr ( f=0.0, pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2 ! Begin loop over atoms

       IF ( j == i ) CYCLE ! Skip self

       rij(:) = r(:,i) - r(:,j)           ! separation vector
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundary conditions in box=1 units
       rij_sq = SUM ( rij**2 )            ! squared separation

       IF ( rij_sq < r_cut_box_sq ) THEN ! check within cutoff

          rij_sq  = rij_sq * box_sq ! now in sigma=1 units
          rij(:)  = rij(:) * box    ! now in sigma=1 units
          sr2     = 1.0 / rij_sq

          IF ( sr2 > sr2_overlap ) THEN
             force%ovr = .TRUE. ! Overlap detected
             RETURN             ! Return immediately
          END IF

          sr6          = sr2 ** 3
          sr12         = sr6 ** 2
          potij        = sr12 - sr6
          force%pot    = force%pot + potij - pot_cut ! shifted potential
          virij        = potij + sr12
          force%vir    = force%vir + virij
          fij          = rij * virij / rij_sq
          force%f(:,i) = force%f(:,i) + fij
          force%f(:,j) = force%f(:,j) - fij

       END IF ! End check within cutoff

    END DO ! End inner loop over atoms

    ! multiply results by numerical factors
    force%f   = force%f   * 24.0
    force%pot = force%pot * 4.0
    force%vir = force%vir * 24.0 / 3.0

    force%ovr = .false. ! No overlaps detected (redundant but for clarity)

  END FUNCTION force_1

END MODULE smc_module
