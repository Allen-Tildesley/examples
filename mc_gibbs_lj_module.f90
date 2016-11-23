! mc_gibbs_lj_module.f90
! Energy and move routines for Gibbs MC, LJ potential
MODULE mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: potential_1, potential, potential_lrc, pressure_lrc, pressure_delta, force_sq
  PUBLIC :: move, create, destroy

  ! Public data
  INTEGER,                              PUBLIC :: n1,n2 ! Number of atoms in each box
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r     ! Positions (3,n1+n2)

  ! Private data
  REAL, DIMENSION(:,:), ALLOCATABLE :: f ! Forces for force_sq calculation (3,n1+n2)

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range

  ! Public derived type
  TYPE, PUBLIC :: potential_type   ! A composite variable for interactions comprising
     REAL    :: pot ! the potential energy cut at r_cut and
     REAL    :: vir ! the virial and
     REAL    :: lap ! the Laplacian and
     LOGICAL :: ovr ! a flag indicating overlap (i.e. pot too high to use)
   CONTAINS
     PROCEDURE :: add_potential_type
     PROCEDURE :: subtract_potential_type
     GENERIC   :: OPERATOR(+) => add_potential_type
     GENERIC   :: OPERATOR(-) => subtract_potential_type
  END TYPE potential_type

CONTAINS

  FUNCTION add_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the sum of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot  +   b%pot
    c%vir = a%vir  +   b%vir
    c%lap = a%lap  +   b%lap
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  FUNCTION subtract_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the difference of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot  -   b%pot
    c%vir = a%vir  -   b%vir
    c%lap = a%lap  -   b%lap
    c%ovr = a%ovr .OR. b%ovr ! This is meaningless, but inconsequential
  END FUNCTION subtract_potential_type

  SUBROUTINE introduction ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut (but not shifted)'
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
    REAL, DIMENSION(2), INTENT(in) :: box   ! Simulation box lengths
    REAL,               INTENT(in) :: r_cut ! Potential cutoff distance

    REAL :: r_cut_box

    ALLOCATE ( r(3,n1+n2), f(3,n1+n2) )

    r_cut_box = r_cut / minval ( box )
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)') 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, f )

  END SUBROUTINE deallocate_arrays

  FUNCTION potential ( i1, i2, box, r_cut ) RESULT ( total )
    IMPLICIT NONE
    TYPE(potential_type) :: total  ! Returns a composite of pot, vir etc
    INTEGER, INTENT(in)  :: i1, i2 ! Index range defining system: (1,n1) or (n1+1,n2)
    REAL,    INTENT(in)  :: box    ! Simulation box length
    REAL,    INTENT(in)  :: r_cut  ! Potential cutoff distance

    ! total%pot is the nonbonded cut (not shifted) potential energy for whole system
    ! total%vir is the corresponding virial for whole system
    ! total%lap is the corresponding Laplacian for whole system
    ! total%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of total%pot etc should not be used
    ! Actual calculation is performed by function potential_1

    TYPE(potential_type) :: partial ! Atomic contribution to total
    INTEGER              :: i

    IF ( n1+n2 > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n1+n2, SIZE(r,dim=2)
       STOP 'Error in potential'
    END IF

    total = potential_type ( pot=0.0, vir=0.0, lap=0.0, ovr=.FALSE. ) ! Initialize

    DO i = i1, i2

       partial = potential_1 ( i1, i2, r(:,i), i, box, r_cut, gt )

       IF ( partial%ovr ) THEN
          total%ovr = .TRUE. ! Overlap detected
          RETURN             ! Return immediately
       END IF

       total = total + partial

    END DO

    total%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential

  FUNCTION potential_1 ( i1, i2, ri, i, box, r_cut, j_range ) RESULT ( partial )
    IMPLICIT NONE
    TYPE(potential_type)            :: partial ! Returns a composite of pot, vir etc for given atom
    INTEGER,            INTENT(in)  :: i1, i2  ! Index range defining system: (1,n1) or (n1+1,n2)
    REAL, DIMENSION(3), INTENT(in)  :: ri      ! Coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i       ! Index of atom of interest
    REAL,               INTENT(in)  :: box     ! Simulation box length
    REAL,               INTENT(in)  :: r_cut   ! Potential cutoff distance
    INTEGER, OPTIONAL,  INTENT(in)  :: j_range ! Optional partner index range

    ! partial%pot is the nonbonded cut (not shifted) potential energy of atom ri with a set of other atoms
    ! partial%vir is the corresponding virial of atom ri
    ! partial%lap is the corresponding Laplacian of atom ri
    ! partial%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the values of partial%pot etc should not be used
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER              :: j, j1, j2
    REAL                 :: r_cut_box, r_cut_box_sq, box_sq
    REAL                 :: sr2, sr6, sr12, rij_sq
    REAL, DIMENSION(3)   :: rij
    REAL, PARAMETER      :: sr2_ovr = 1.77 ! overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    IF ( n1+n2 > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n1+n2, SIZE(r,dim=2)
       STOP 'Error in potential_1'
    END IF
    IF ( i < i1 .or. i > i2 ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,3i15)' ) 'Indexing error', i1, i2, i
       STOP 'Error in potential_1'
    END IF

    IF ( PRESENT ( j_range ) ) THEN
       SELECT CASE ( j_range )
       CASE ( lt ) ! j < i
          j1 = i1
          j2 = i-1
       CASE ( gt ) ! j > i
          j1 = i+1
          j2 = i2
       CASE default ! should never happen
          WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error ', j_range
          STOP 'Impossible error in potential_1'
       END SELECT
    ELSE
       j1 = i1
       j2 = i2
    END IF

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box**2
    box_sq       = box**2

    partial = potential_type ( pot=0.0, vir=0.0, lap=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( i == j ) CYCLE ! Skip self

       rij(:) = ri(:) - r(:,j)            ! Separation vector
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )            ! Squared separation in box=1 units

       IF ( rij_sq < r_cut_box_sq ) THEN ! Check within range

          rij_sq   = rij_sq * box_sq ! Now in sigma=1 units
          sr2      = 1.0 / rij_sq    ! (sigma/rij)**2
          pair%ovr = sr2 > sr2_ovr   ! Overlap if too close

          IF ( pair%ovr ) THEN
             partial%ovr = .TRUE. ! Overlap detected
             RETURN               ! Return immediately
          END IF

          sr6      = sr2**3
          sr12     = sr6**2
          pair%pot = sr12 - sr6                    ! LJ pair potential (cut but not shifted)
          pair%vir = pair%pot + sr12               ! LJ pair virial
          pair%lap = ( 22.0*sr12 - 5.0*sr6 ) * sr2 ! LJ pair Laplacian

          partial = partial + pair 

       END IF ! End check within range

    END DO ! End loop over selected range of partners

    ! Include numerical factors
    partial%pot = partial%pot * 4.0        ! 4*epsilon
    partial%vir = partial%vir * 24.0 / 3.0 ! 24*epsilon and divide virial by 3
    partial%lap = partial%lap * 24.0 * 2.0 ! 24*epsilon and factor 2 for ij and ji
    partial%ovr = .FALSE.                  ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential_1

  FUNCTION force_sq ( i1, i2, box, r_cut ) RESULT ( fsq )
    IMPLICIT NONE
    REAL                :: fsq    ! Returns total squared force
    INTEGER, INTENT(in) :: i1, i2 ! Index range defining system: (1,n1) or (n1+1,n2)
    REAL,    INTENT(in) :: box    ! Simulation box length
    REAL,    INTENT(in) :: r_cut  ! Potential cutoff distance

    ! Calculates total squared force (using array f)

    INTEGER            :: i, j
    REAL               :: r_cut_box, r_cut_box_sq, box_sq, rij_sq
    REAL               :: sr2, sr6, sr12
    REAL, DIMENSION(3) :: rij, fij

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box ** 2
    box_sq       = box ** 2

    f(:,i1:i2) = 0.0 ! Initialize

    DO i = i1, i2 - 1 ! Begin outer loop over atoms

       DO j = i + 1, i2 ! Begin inner loop over atoms

          rij(:) = r(:,i) - r(:,j)           ! Separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )            ! Squared separation

          IF ( rij_sq < r_cut_box_sq ) THEN ! Check within cutoff

             rij_sq = rij_sq * box_sq ! Now in sigma=1 units
             rij(:) = rij(:) * box    ! Now in sigma=1 units
             sr2    = 1.0 / rij_sq
             sr6    = sr2 ** 3
             sr12   = sr6 ** 2
             fij    = rij * (2.0*sr12 - sr6) *sr2 ! LJ pair forces
             f(:,i) = f(:,i) + fij
             f(:,j) = f(:,j) - fij

          END IF ! End check within cutoff

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    f(:,i1:i2) = f(:,i1:i2) * 24.0     ! Numerical factor 24*epsilon
    fsq        = SUM ( f(:,i1:i2)**2 ) ! Result

  END FUNCTION force_sq

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

  SUBROUTINE swap ( i, ri )
    IMPLICIT NONE
    INTEGER,            INTENT(in) :: i  ! Index of particle to be destroyed
    REAL, DIMENSION(3), INTENT(in) :: ri ! Position of particle to be created

    IF ( i <= n1 ) THEN ! Destroy in system 1, create in system 2
       r(:,i)    = r(:,n1) ! Replace coordinates of i with those of N1
       n1        = n1 - 1  ! Decrease N1
       n2        = n2 + 1  ! Increase N2
       r(:,n1+1) = ri      ! New particle coordinates at start of system 2
    ELSE ! Destroy in system 2, create in system 1
       r(:,i)  = r(:,n1+1) ! Replace coordinates of i with those of N1+1
       n1      = n1 + 1    ! Increase N1
       n2      = n2 - 1    ! Decrease N2
       r(:,n1) = ri        ! New particle coordinates at end of system 1
    END IF
     
  END SUBROUTINE swap

END MODULE mc_module
