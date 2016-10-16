! mc_lj_module.f90
! Energy and move routines for MC simulation, LJ potential
MODULE mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, r
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: resize, energy_1, energy, energy_lrc, pressure_lrc, pressure_delta
  PUBLIC :: move, create, destroy
  PUBLIC :: pot_type

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,n)

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range
  REAL,    PARAMETER :: sigma = 1.0     ! Lennard-Jones diameter (unit of length)
  REAL,    PARAMETER :: epslj = 1.0     ! Lennard-Jones well depth (unit of energy)

  TYPE pot_type ! A composite variable for interaction energies comprising
     REAL    :: pot ! the potential energy and
     REAL    :: vir ! the virial and
     LOGICAL :: ovr ! a flag indicating overlap (i.e. pot too high to use)
  END TYPE pot_type

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_pot_type
  END INTERFACE OPERATOR (+)
  
CONTAINS

  FUNCTION add_pot_type ( a, b ) RESULT (c)
    TYPE(pot_type)             :: c    ! Result is the sum of the two inputs
    TYPE(pot_type), INTENT(in) :: a, b
    c%pot = a%pot +    b%pot
    c%vir = a%vir +    b%vir
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_pot_type
  
  SUBROUTINE introduction ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ',   sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epslj = ', epslj    
  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output
    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'
  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

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

  SUBROUTINE resize 

    ! Reallocates r array, twice as large
    ! This is employed by mc_zvt_lj, grand canonical ensemble

    REAL, DIMENSION(:,:), ALLOCATABLE :: tmp
    INTEGER                           :: n_old, n_new

    n_old = SIZE(r,dim=2)
    n_new = 2*n_old
    WRITE( unit=output_unit, fmt='(a,i10,a,i10)' ) 'Reallocating r from old ', n_old, ' to ', n_new

    ALLOCATE ( tmp(3,n_new) ) ! New size for r
    tmp(:,1:n_old) = r(:,:)   ! Copy elements across

    CALL move_ALLOC ( tmp, r )

  END SUBROUTINE resize

  FUNCTION energy ( box, r_cut )
    TYPE(pot_type)     :: energy ! Returns a composite of pot, vir and ovr
    REAL, INTENT(in) :: box    ! Simulation box length
    REAL, INTENT(in) :: r_cut  ! Potential cutoff distance

    ! energy%pot is the nonbonded potential energy for whole system
    ! energy%vir is the corresponding virial for whole system
    ! energy%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of energy%pot, energy%vir should not be used
    ! Actual calculation is performed by function energy_1

    TYPE(pot_type) :: energy_i
    INTEGER      :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy'
    END IF

    energy = pot_type ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize
    
    DO i = 1, n - 1
       energy_i = energy_1 ( r(:,i), i, box, r_cut, gt )
       IF ( energy_i%ovr ) THEN
          energy%ovr = .TRUE. ! Overlap detected
          RETURN              ! Return immediately
       END IF
       energy = energy + energy_i
    END DO

    energy%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)
    
  END FUNCTION energy

  FUNCTION energy_1 ( ri, i, box, r_cut, j_range ) RESULT ( energy )
    TYPE(pot_type)                    :: energy  ! Returns a composite of pot, vir and ovr
    REAL, DIMENSION(3), INTENT(in)  :: ri      ! Coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i       ! Index of atom of interest
    REAL,               INTENT(in)  :: box     ! Simulation box length
    REAL,               INTENT(in)  :: r_cut   ! Potential cutoff distance
    INTEGER, OPTIONAL,  INTENT(in)  :: j_range ! Optional partner index range

    ! energy%pot is the nonbonded potential energy of atom ri with a set of other atoms
    ! energy%vir is the corresponding virial of atom ri
    ! energy%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the value of energy%pot should not be used
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

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
          STOP 'Impossible error in energy_1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box**2
    box_sq       = box**2

    energy = pot_type ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )

       IF ( rij_sq < r_cut_box_sq ) THEN

          rij_sq = rij_sq * box_sq ! now in sigma=1 units
          sr2    = 1.0 / rij_sq    ! (sigma/rij)**2

          IF ( sr2 > sr2_overlap ) THEN
             energy%ovr = .TRUE. ! Overlap detected
             RETURN              ! Return immediately
          END IF

          sr6 = sr2**3
          energy%pot = energy%pot + sr6**2 - sr6
          energy%vir = energy%vir + 2.0*sr6**2 - sr6

       END IF

    END DO

    ! Include numerical factors
    energy%pot = 4.0 * energy%pot
    energy%vir = 24.0 * energy%vir
    energy%vir = energy%vir / 3.0

    energy%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION energy_1

  FUNCTION energy_lrc ( density, r_cut )
    REAL                :: energy_lrc ! Returns long-range energy/atom
    REAL,    INTENT(in) :: density    ! Number density N/V
    REAL,    INTENT(in) :: r_cut      ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones energy per atom
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3        = 1.0 / r_cut**3
    energy_lrc = pi * ( (8.0/9.0)  * sr3**3  - (8.0/3.0)  * sr3 ) * density

  END FUNCTION energy_lrc

  FUNCTION pressure_lrc ( density, r_cut )
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

  FUNCTION pressure_delta ( density, r_cut )
    REAL                :: pressure_delta ! Returns pressure delta correction
    REAL,    INTENT(in) :: density        ! Number density N/V
    REAL,    INTENT(in) :: r_cut          ! Cutoff distance

    ! Calculates correction for Lennard-Jones pressure
    ! due to discontinuity in the potential at r_cut
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3            = 1.0 / r_cut**3
    pressure_delta = pi * (8.0/3.0) * ( sr3**3  - sr3 ) * density**2

  END FUNCTION pressure_delta

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
