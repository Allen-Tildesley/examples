! mc_gibbs_lj_module.f90
! Energy and move routines for Gibbs MC, LJ potential
MODULE mc_module

  ! TODO (DJT) adapt for Gibbs simulation

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: potential_1, potential, potential_lrc, pressure_lrc, pressure_delta
  PUBLIC :: move, create, destroy

  ! Public data
  INTEGER,                              PUBLIC :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! positions (3,n)

  ! Private data
  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range

  ! Public derived type
  PUBLIC :: potential_type
  TYPE potential_type   ! A composite variable for interactions comprising
     REAL    :: pot_c   ! the potential energy cut at r_cut and
     REAL    :: pot_s   ! the potential energy cut and shifted to 0 at r_cut, and
     REAL    :: vir     ! the virial and
     LOGICAL :: overlap ! a flag indicating overlap (i.e. pot too high to use)
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
    c%pot_c   = a%pot_c    +   b%pot_c
    c%pot_s   = a%pot_s    +   b%pot_s
    c%vir     = a%vir      +   b%vir
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  FUNCTION subtract_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)             :: c    ! Result is the difference of the two inputs
    TYPE(potential_type), INTENT(in) :: a, b
    c%pot_c   = a%pot_c    -   b%pot_c
    c%pot_s   = a%pot_s    -   b%pot_s
    c%vir     = a%vir      -   b%vir
    c%ovr = a%ovr .OR. b%ovr ! this is meaningless, but inconsequential
  END FUNCTION subtract_potential_type

  SUBROUTINE introduction ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut and optionally shifted'
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

    ALLOCATE ( r(3,n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r )

  END SUBROUTINE deallocate_arrays

  FUNCTION potential ( box, r_cut ) RESULT ( total )
    IMPLICIT NONE
    TYPE(potential_type) :: total ! Returns a composite of pot, vir and overlap
    REAL, INTENT(in)     :: box   ! Simulation box length
    REAL, INTENT(in)     :: r_cut ! Potential cutoff distance

    ! total%pot_c is the nonbonded cut (not shifted) potential energy for whole system
    ! total%pot_s is the nonbonded cut-and-shifted potential energy for whole system
    ! total%vir is the corresponding virial for whole system
    ! total%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of total%pot_c etc should not be used
    ! Actual calculation is performed by function potential_1

    TYPE(potential_type) :: partial ! Atomic contribution to total
    INTEGER              :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in potential'
    END IF

    total = potential_type ( pot_c=0.0, pot_s=0.0, vir=0.0, overlap=.FALSE. ) ! Initialize

    DO i = 1, n - 1

       partial = potential_1 ( r(:,i), i, box, r_cut, gt )

       IF ( partial%ovr ) THEN
          total%ovr = .TRUE. ! Overlap detected
          RETURN                 ! Return immediately
       END IF

       total = total + partial

    END DO

    total%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential

  FUNCTION potential_1 ( ri, i, box, r_cut, j_range ) RESULT ( partial )
    IMPLICIT NONE
    TYPE(potential_type)            :: partial ! Returns a composite of pot, vir and overlap for given atom
    REAL, DIMENSION(3), INTENT(in)  :: ri      ! Coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i       ! Index of atom of interest
    REAL,               INTENT(in)  :: box     ! Simulation box length
    REAL,               INTENT(in)  :: r_cut   ! Potential cutoff distance
    INTEGER, OPTIONAL,  INTENT(in)  :: j_range ! Optional partner index range

    ! partial%pot_c is the nonbonded cut (not shifted) potential energy of atom ri with a set of other atoms
    ! partial%pot_s is the nonbonded cut-and-shifted potential energy of atom ri with a set of other atoms
    ! partial%vir is the corresponding virial of atom ri
    ! partial%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the values of partial%pot_c etc should not be used
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: j, j1, j2, ncut
    REAL               :: r_cut_box, r_cut_box_sq, box_sq
    REAL               :: sr2, sr6, sr12, rij_sq, potij, virij
    REAL, DIMENSION(3) :: rij
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in potential_1'
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
          STOP 'Impossible error in potential_1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box**2
    box_sq       = box**2

    partial = potential_type ( pot_c=0.0, pot_s=0.0, vir=0.0, overlap=.FALSE. ) ! Initialize
    ncut    = 0

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( i == j ) CYCLE ! Skip self

       rij(:) = ri(:) - r(:,j)            ! Separation vector
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )            ! Squared separation in box=1 units

       IF ( rij_sq < r_cut_box_sq ) THEN ! Check within range

          rij_sq = rij_sq * box_sq ! now in sigma=1 units
          sr2    = 1.0 / rij_sq    ! (sigma/rij)**2

          IF ( sr2 > sr2_overlap ) THEN
             partial%ovr = .TRUE. ! Overlap detected
             RETURN                ! Return immediately
          END IF

          sr6   = sr2**3
          sr12  = sr6**2
          potij = sr12 - sr6   ! LJ potential (cut but not shifted)
          virij = potij + sr12 ! LJ virial

          partial%pot_c = partial%pot_c + potij 
          partial%vir   = partial%vir   + virij 

          ncut = ncut + 1

       END IF ! End check within range

    END DO ! End loop over selected range of partners

    ! Calculate shifted potential
    sr2   = 1.0 / r_cut**2 ! in sigma=1 units
    sr6   = sr2**3
    sr12  = sr6**2
    potij = sr12 - sr6
    partial%pot_s = partial%pot_c - REAL ( ncut ) * potij

    ! Include numerical factors
    partial%pot_c   = partial%pot_c * 4.0
    partial%pot_s   = partial%pot_s * 4.0
    partial%vir     = partial%vir   * 24.0 / 3.0
    partial%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential_1

  FUNCTION potential_lrc ( density, r_cut )
    IMPLICIT NONE
    REAL                :: potential_lrc ! Returns long-range correction to potential/atom
    REAL,    INTENT(in) :: density       ! Number density N/V
    REAL,    INTENT(in) :: r_cut         ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones potential per atom
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3 = 1.0 / r_cut**3

    potential_lrc = pi * ( (8.0/9.0)  * sr3**3  - (8.0/3.0)  * sr3 ) * density

  END FUNCTION potential_lrc

  FUNCTION pressure_lrc ( density, r_cut )
    IMPLICIT NONE
    REAL                :: pressure_lrc ! Returns long-range correction to pressure
    REAL,    INTENT(in) :: density      ! Number density N/V
    REAL,    INTENT(in) :: r_cut        ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones pressure
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3 = 1.0 / r_cut**3

    pressure_lrc = pi * ( (32.0/9.0) * sr3**3  - (16.0/3.0) * sr3 ) * density**2

  END FUNCTION pressure_lrc

  FUNCTION pressure_delta ( density, r_cut )
    IMPLICIT NONE
    REAL                :: pressure_delta ! Returns delta correction to pressure
    REAL,    INTENT(in) :: density        ! Number density N/V
    REAL,    INTENT(in) :: r_cut          ! Cutoff distance

    ! Calculates correction for Lennard-Jones pressure
    ! due to discontinuity in the potential at r_cut
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3 = 1.0 / r_cut**3

    pressure_delta = pi * (8.0/3.0) * ( sr3**3  - sr3 ) * density**2

  END FUNCTION pressure_delta

  SUBROUTINE move ( i, ri )
    IMPLICIT NONE
    INTEGER,               INTENT(in) :: i
    REAL,    DIMENSION(3), INTENT(in) :: ri

    r(:,i) = ri ! New position

  END SUBROUTINE move

  SUBROUTINE create ( ri )
    IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(in) :: ri

    n      = n+1 ! Increase number of atoms
    r(:,n) = ri  ! Add new atom at the end

  END SUBROUTINE create

  SUBROUTINE destroy ( i )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: i

    r(:,i) = r(:,n) ! Replace atom i coordinates with atom n
    n      = n - 1  ! Reduce number of atoms

  END SUBROUTINE destroy

END MODULE mc_module
