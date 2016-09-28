! mc_lj_module.f90
! Energy and move routines for MC simulation, LJ potential
MODULE mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, r, lt, ne, gt
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  public :: resize, energy_1, energy, energy_lrc
  PUBLIC :: move, create, destroy

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,n)

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options
  REAL,    PARAMETER :: sigma = 1.0             ! Lennard-Jones diameter (unit of length)
  REAL,    PARAMETER :: epslj = 1.0             ! Lennard-Jones well depth (unit of energy)

CONTAINS

  SUBROUTINE introduction ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ',   sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epslj = ', epslj    
  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output
    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'
  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
    REAL, INTENT(in) :: box   ! simulation box length
    REAL, INTENT(in) :: r_cut ! potential cutoff distance

    REAL :: r_cut_box
    
    ALLOCATE ( r(3,n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE resize ! reallocates r array, twice as large

    ! This is employed by mc_zvt_lj, grand canonical ensemble
    
    REAL, DIMENSION(:,:), ALLOCATABLE :: tmp
    INTEGER                           :: n_old, n_new

    n_old = SIZE(r,dim=2)
    n_new = 2*n_old
    WRITE( unit=output_unit, fmt='(a,i10,a,i10)' ) 'Reallocating r from old ', n_old, ' to ', n_new

    ALLOCATE ( tmp(3,n_new) ) ! new size for r
    tmp(:,1:n_old) = r(:,:)   ! copy elements across

    CALL move_ALLOC ( tmp, r )

  END SUBROUTINE resize

  SUBROUTINE energy ( box, r_cut, overlap, pot, vir )
    REAL,    INTENT(in)  :: box        ! simulation box length
    REAL,    INTENT(in)  :: r_cut      ! potential cutoff distance
    LOGICAL, INTENT(out) :: overlap    ! shows if an overlap was detected
    REAL,    INTENT(out) :: pot, vir   ! potential and virial 

    ! Calculates potential and virial for whole system
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the values of pot and vir should not be used
    ! Actual calculation is performed by subroutine energy_1

    REAL               :: pot_i, vir_i
    INTEGER            :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy'
    END IF
    
    overlap  = .FALSE.
    pot      = 0.0
    vir      = 0.0

    DO i = 1, n - 1
       CALL energy_1 ( r(:,i), i, gt, box, r_cut, overlap, pot_i, vir_i )
       IF ( overlap ) EXIT ! jump out of loop
       pot  = pot  + pot_i
       vir  = vir  + vir_i
    END DO
    
  END SUBROUTINE energy

  SUBROUTINE energy_1 ( ri, i, j_range, box, r_cut, overlap, pot, vir )

    REAL, DIMENSION(3), INTENT(in)  :: ri         ! coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i, j_range ! index, and partner index range
    REAL,               INTENT(in)  :: box        ! simulation box length
    REAL,               INTENT(in)  :: r_cut      ! potential cutoff distance
    LOGICAL,            INTENT(out) :: overlap    ! shows if an overlap was detected
    REAL,               INTENT(out) :: pot, vir   ! potential and virial

    ! Calculates potential energy and virial of atom in ri
    ! with j/=i, j>i, or j<i depending on j_range
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the values of pot and vir should not be used
    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: j, j1, j2
    REAL               :: r_cut_box, r_cut_box_sq, box_sq
    REAL               :: sr2, sr6, rij_sq
    REAL, DIMENSION(3) :: rij
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy_1'
    END IF

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box**2
    box_sq       = box**2

    pot     = 0.0
    vir     = 0.0
    overlap = .FALSE.

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

    DO j = j1, j2

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )

       IF ( rij_sq < r_cut_box_sq ) THEN

          rij_sq = rij_sq * box_sq ! now in sigma=1 units
          sr2    = 1.0 / rij_sq    ! (sigma/rij)**2

          IF ( sr2 > sr2_overlap ) THEN
             overlap = .TRUE.
             EXIT ! jump out of loop
          END IF

          sr6 = sr2**3
          pot = pot + sr6**2 - sr6
          vir = vir + 2.0*sr6**2 - sr6

       END IF

    END DO
    pot = 4.0 * pot
    vir = 24.0 * vir
    vir = vir / 3.0

  END SUBROUTINE energy_1

  SUBROUTINE energy_lrc ( n, box, r_cut, pot, vir )
    INTEGER, INTENT(in)  :: n        ! number of atoms
    REAL,    INTENT(in)  :: box      ! simulation box length
    REAL,    INTENT(in)  :: r_cut    ! cutoff distance
    REAL,    INTENT(out) :: pot, vir ! potential and virial

    ! Calculates long-range corrections for Lennard-Jones potential and virial
    ! These are the corrections to the total values
    ! r_cut, box, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL               :: sr3, density
    REAL, PARAMETER    :: pi = 4.0 * ATAN(1.0)

    sr3 = ( 1.0 / r_cut ) ** 3
    pot = (8.0/9.0)  * sr3**3  - (8.0/3.0)  * sr3
    vir = (32.0/9.0) * sr3**3  - (32.0/6.0) * sr3

    density = REAL(n) / box**3
    pot     = pot * pi * density * REAL(n)
    vir     = vir * pi * density * REAL(n)

  END SUBROUTINE energy_lrc

  SUBROUTINE move ( i, ri )
    INTEGER,               INTENT(in) :: i
    REAL,    DIMENSION(3), INTENT(in) :: ri

    r(:,i) = ri
    
  END SUBROUTINE move

  SUBROUTINE create ( ri )
    REAL,    DIMENSION(3), INTENT(in) :: ri

    n        = n+1 ! increase number of atoms
    r(:,n)   = ri  ! add new atom at the end
    
  END SUBROUTINE create

  SUBROUTINE destroy ( i )
    INTEGER, INTENT(in) :: i

    r(:,i)    = r(:,n) ! replace atom i with atom n
    n         = n - 1  ! reduce number of atoms

  END SUBROUTINE destroy
  
END MODULE mc_module
