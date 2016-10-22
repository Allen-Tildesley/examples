! smc_lj_module.f90
! Energy, force, and move routines for SMC, LJ potential
MODULE smc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: force, force_1, potential_lrc, pressure_lrc

  ! Public data
  INTEGER,                              PUBLIC :: n     ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r     ! positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r_old ! old positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v     ! velocities (3,n)
  REAL,    DIMENSION(:),   ALLOCATABLE, PUBLIC :: zeta  ! random numbers (n)
  LOGICAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: move  ! mask for multi-atom moves (3,n)

  ! Private data
  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range
  REAL,    PARAMETER :: sigma = 1.0     ! LJ diameter (unit of length)
  REAL,    PARAMETER :: epslj = 1.0     ! LJ well depth (unit of energy)

  ! Public derived type
  ! At the time of writing gfortran, and some other compilers, do not implement
  ! parameterized derived types (part of the Fortran 2003 standard);
  ! hence the rather clumsy use of the n_max parameter here.
  PUBLIC :: potential_type
  INTEGER, PARAMETER :: n_max = 256 
  TYPE potential_type ! A composite variable for interactions comprising
     REAL                     :: pot_s   ! the cut-and-shifted potential energy and
     REAL                     :: pot_c   ! the cut (but not shifted) potential energy and
     REAL                     :: vir     ! the virial and
     REAL, DIMENSION(3,n_max) :: f       ! the forces and
     LOGICAL                  :: overlap ! a flag indicating overlap (i.e. pot_s too high to use)
  END TYPE potential_type

  PUBLIC :: OPERATOR (+)
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_potential_type
  END INTERFACE OPERATOR (+)

  PUBLIC :: OPERATOR (-)
  INTERFACE OPERATOR (-)
     MODULE PROCEDURE subtract_potential_type
  END INTERFACE OPERATOR (-)

CONTAINS

  FUNCTION add_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)             :: c    ! Result is the sum of the two inputs
    TYPE(potential_type), INTENT(in) :: a, b
    c%pot_s   = a%pot_s    +   b%pot_s
    c%pot_c   = a%pot_c    +   b%pot_c
    c%vir     = a%vir      +   b%vir
    c%f       = a%f        +   b%f
    c%overlap = a%overlap .OR. b%overlap
  END FUNCTION add_potential_type

  FUNCTION subtract_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)             :: c    ! Result is the difference of the two inputs
    TYPE(potential_type), INTENT(in) :: a, b
    c%pot_s   = a%pot_s    -   b%pot_s
    c%pot_c   = a%pot_c    -   b%pot_c
    c%vir     = a%vir      -   b%vir
    c%f       = a%f        -   b%f
    c%overlap = a%overlap .OR. b%overlap ! meaningless but inconsequential
  END FUNCTION subtract_potential_type

  SUBROUTINE introduction ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)'           ) 'Cut and optionally shifted'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ',     sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epsilon = ', epslj    

  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
    IMPLICIT NONE
    REAL, INTENT(in) :: box   ! simulation box length
    REAL, INTENT(in) :: r_cut ! potential cutoff distance

    REAL :: r_cut_box

    IF ( n > n_max ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i5)' ) 'n too large ', n, n_max
       STOP 'Error in allocate_arrays'
    END IF

    ALLOCATE ( r(3,n), r_old(3,n), v(3,n), zeta(n), move(3,n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, r_old, v, zeta, move )

  END SUBROUTINE deallocate_arrays

  FUNCTION force ( box, r_cut ) RESULT ( system )
    IMPLICIT NONE
    TYPE(potential_type) :: system ! Returns composite of forces, pot, vir and overlap
    REAL, INTENT(in)     :: box    ! simulation box length
    REAL, INTENT(in)     :: r_cut  ! potential cutoff distance

    ! system%pot_s is the cut-and-shifted potential energy for whole system
    ! system%pot_c is the cut (but not shifted) version of the above
    ! system%vir is the corresponding virial for whole system
    ! system%f contains the forces on all the atoms
    ! system%overlap is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of system%pot_s etc should not be used
    ! Actual calculation is performed by function force_1

    INTEGER              :: i
    TYPE(potential_type) :: atom

    system = potential_type ( f=0.0, pot_s=0.0, pot_c=0.0, vir=0.0, overlap=.FALSE. ) ! Initialize

    DO i = 1, n - 1 ! Begin loop over atoms

       atom = force_1 ( i, box, r_cut, gt )

       IF ( atom%overlap ) THEN
          system%overlap = .TRUE. ! Overlap detected
          RETURN                  ! Return immediately
       END IF
       system = system + atom

    END DO ! End loop over atoms

    system%overlap = .FALSE. ! No overlaps detected (redundant but for clarity)

  END FUNCTION force

  FUNCTION force_1 ( i, box, r_cut, j_range ) RESULT ( atom )
    IMPLICIT NONE
    TYPE(potential_type)          :: atom    ! Returns composite of forces, pot, vir and overlap
    INTEGER,           INTENT(in) :: i       ! index of atom of interest
    REAL,              INTENT(in) :: box     ! simulation box length
    REAL,              INTENT(in) :: r_cut   ! potential cutoff distance
    INTEGER, OPTIONAL, INTENT(in) :: j_range ! Optional partner index range

    ! atom%pot_s is the cut-and-shifted potential energy of atom i with a set of other atoms
    ! atom%pot_c is the cut (but not shifted) version of the above
    ! atom%vir is the corresponding virial of atom i
    ! atom%f contains the force on i and the reaction forces on all other atoms due to i
    ! atom%overlap is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the values of atom%pot_s etc should not be used
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that positions are in units where box = 1
    ! Forces are calculated in units where sigma = 1 and epsilon = 1
    ! Note that we use a shifted LJ potential here

    INTEGER            :: j, j1, j2, ncut
    REAL               :: r_cut_box, r_cut_box_sq, box_sq
    REAL               :: rij_sq, sr2, sr6, sr12, potij, virij
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

    atom = potential_type ( f=0.0, pot_s=0.0, pot_c=0.0, vir=0.0, overlap=.FALSE. ) ! Initialize
    ncut = 0

    DO j = j1, j2 ! Begin loop over atoms

       IF ( j == i ) CYCLE ! Skip self

       rij(:) = r(:,i) - r(:,j)           ! Separation vector
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundary conditions in box=1 units
       rij_sq = SUM ( rij**2 )            ! Squared separation

       IF ( rij_sq < r_cut_box_sq ) THEN ! Check within cutoff

          rij_sq  = rij_sq * box_sq ! now in sigma=1 units
          rij(:)  = rij(:) * box    ! now in sigma=1 units
          sr2     = 1.0 / rij_sq

          IF ( sr2 > sr2_overlap ) THEN
             atom%overlap = .TRUE. ! Overlap detected
             RETURN                ! Return immediately
          END IF

          sr6          = sr2 ** 3
          sr12         = sr6 ** 2
          potij        = sr12 - sr6
          atom%pot_c   = atom%pot_c + potij ! LJ potential (cut but not shifted)
          virij        = potij + sr12
          atom%vir     = atom%vir + virij
          fij          = rij * virij / rij_sq
          atom%f(:,i)  = atom%f(:,i) + fij
          atom%f(:,j)  = atom%f(:,j) - fij
          ncut         = ncut + 1

       END IF ! End check within cutoff

    END DO ! End inner loop over atoms

    ! Calculate shifted potential
    sr2        = 1.0 / r_cut**2 ! in sigma=1 units
    sr6        = sr2 ** 3
    sr12       = sr6 **2
    atom%pot_s = atom%pot_c - real ( ncut ) * ( sr12 - sr6 )

    ! Multiply results by numerical factors
    atom%f        = atom%f * 24.0
    atom%pot_s    = atom%pot_s * 4.0
    atom%pot_c    = atom%pot_c * 4.0
    atom%vir      = atom%vir * 24.0 / 3.0
    atom%overlap  = .FALSE. ! No overlaps detected (redundant but for clarity)

  END FUNCTION force_1

  FUNCTION potential_lrc ( density, r_cut )
    IMPLICIT NONE
    REAL                :: potential_lrc ! Returns long-range energy/atom
    REAL,    INTENT(in) :: density       ! Number density N/V
    REAL,    INTENT(in) :: r_cut         ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones energy per atom
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3 = 1.0 / r_cut**3

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

    sr3 = 1.0 / r_cut**3

    pressure_lrc = pi * ( (32.0/9.0) * sr3**3  - (16.0/3.0) * sr3 ) * density**2

  END FUNCTION pressure_lrc

END MODULE smc_module
