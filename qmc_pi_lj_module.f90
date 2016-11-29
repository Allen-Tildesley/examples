! qmc_pi_lj_module.f90
! Energy and move routines for PIMC simulation, LJ potential
MODULE qmc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: potential_1, potential, spring_1, spring, potential_lrc

  ! Public data
  INTEGER,                                PUBLIC :: n ! Number of atoms
  INTEGER,                                PUBLIC :: p ! Number of beads
  REAL,    DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n,p)

  ! Private data
  INTEGER, PARAMETER :: lt = -1, gt = 1 ! j-range and l-range options

  ! Public derived type
  TYPE, PUBLIC :: potential_type ! A composite variable for interaction energies comprising
     REAL    :: pot ! the cut (but not shifted) potential energy and
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
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  FUNCTION subtract_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the difference of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot  -   b%pot
    c%ovr = a%ovr .OR. b%ovr ! This is meaningless, but inconsequential
  END FUNCTION subtract_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut (but not shifted)'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'   
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1' 

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
    IMPLICIT NONE
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

    REAL :: r_cut_box

    ALLOCATE ( r(3,n,p) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)') 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r )

  END SUBROUTINE deallocate_arrays

  FUNCTION potential ( box, r_cut ) RESULT ( total )
    IMPLICIT NONE
    TYPE(potential_type) :: total ! Returns a composite of pot and ovr
    REAL, INTENT(in)     :: box   ! Simulation box length
    REAL, INTENT(in)     :: r_cut ! Potential cutoff distance

    ! total%pot is the nonbonded cut (not shifted) classical potential energy for whole system
    ! total%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of total%pot etc should not be used
    ! Actual calculation is performed by function potential_1

    TYPE(potential_type) :: partial ! Atomic contribution to total
    INTEGER              :: i, k

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in potential'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in potential'
    END IF

    total = potential_type ( pot=0.0, ovr=.FALSE. ) ! Initialize

    DO k = 1, p ! Loop over ring polymers
       DO i = 1, n - 1 ! Loop over atoms within polymer

          partial = potential_1 ( r(:,i,k), i, k, box, r_cut, gt )

          IF ( partial%ovr ) THEN
             total%ovr = .TRUE. ! Overlap detected
             RETURN             ! Return immediately
          END IF

          total = total + partial

       END DO ! End loop over atoms within polymer
    END DO ! End loop over ring polymers

    total%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential

  FUNCTION potential_1 ( rik, i, k, box, r_cut, j_range ) RESULT ( partial )
    IMPLICIT NONE
    TYPE(potential_type)           :: partial ! Returns a composite of pot and ovr for given atom
    REAL, DIMENSION(3), INTENT(in) :: rik     ! Coordinates of atom of interest
    INTEGER,            INTENT(in) :: i, k    ! Index and polymer id of atom of interest
    REAL,               INTENT(in) :: box     ! Simulation box length
    REAL,               INTENT(in) :: r_cut   ! Potential cutoff distance
    INTEGER, OPTIONAL,  INTENT(in) :: j_range ! Optional partner index range

    ! partial%pot is the nonbonded cut (not shifted) classical potential energy of atom rik with a set of other atoms
    ! partial%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the value of partial%pot should not be used
    ! The coordinates in rik are not necessarily identical with those in r(:,i,k)
    ! The partner atoms always have the same polymer index k
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER              :: j, j1, j2
    REAL                 :: r_cut_box, r_cut_box_sq, box_sq
    REAL                 :: sr2, sr6, r_ik_jk_sq
    REAL, DIMENSION(3)   :: r_ik_jk
    REAL, PARAMETER      :: sr2_ovr = 1.77 ! overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in potential_1'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
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

    partial = potential_type ( pot=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( i == j ) CYCLE ! Skip self

       r_ik_jk(:) = rik(:) - r(:,j,k)                 ! Separation vector
       r_ik_jk(:) = r_ik_jk(:) - ANINT ( r_ik_jk(:) ) ! Periodic boundaries in box=1 units
       r_ik_jk_sq = SUM ( r_ik_jk**2 )                ! Squared separation in box=1 units

       IF ( r_ik_jk_sq < r_cut_box_sq ) THEN ! Check within range

          r_ik_jk_sq = r_ik_jk_sq * box_sq ! now in sigma=1 units
          sr2        = 1.0 / r_ik_jk_sq    ! (sigma/rikjk)**2
          pair%ovr   = sr2 > sr2_ovr       ! Overlap if too close

          IF ( pair%ovr ) THEN
             partial%ovr = .TRUE. ! Overlap detected
             RETURN               ! Return immediately
          END IF

          sr6      = sr2**3
          pair%pot = sr6**2 - sr6 ! LJ potential (cut but not shifted)

          partial = partial + pair

       END IF ! End check with range

    END DO ! End loop over selected range of partners

    ! Include numerical factors
    partial%pot = partial%pot * 4.0 / REAL ( p ) ! Classical potentials are weaker by a factor p
    partial%ovr = .FALSE.                        ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential_1

  FUNCTION spring ( box, k_spring ) RESULT ( total )
    IMPLICIT NONE
    REAL             :: total    ! Returns quantum spring potential 
    REAL, INTENT(in) :: box      ! Simulation box length
    REAL, INTENT(in) :: k_spring ! Spring potential constant

    ! total is the quantum spring potential for whole system
    ! Actual calculation is performed by subroutine spring_1

    INTEGER :: i, k
    REAL    :: partial_pot

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in spring'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in spring'
    END IF

    total = 0.0

    DO k = 1, p ! Loop over ring polymers
       DO i = 1, n ! Loop over atoms within polymer

          partial_pot = spring_1 ( r(:,i,k), i, k, box, k_spring, gt )
          total       = total + partial_pot

       END DO ! End loop over atoms within polymer
    END DO ! End loop over ring polymers

  END FUNCTION spring

  FUNCTION spring_1 ( rik, i, k, box, k_spring, l_range ) RESULT ( partial_pot )
    IMPLICIT NONE
    REAL                           :: partial_pot ! Returns quantum potential for given atom
    REAL, DIMENSION(3), INTENT(in) :: rik         ! Coordinates of atom of interest
    INTEGER,            INTENT(in) :: i, k        ! Index and polymer id of atom of interest
    REAL,               INTENT(in) :: box         ! Simulation box length
    REAL,               INTENT(in) :: k_spring    ! Spring potential constant
    INTEGER, OPTIONAL,  INTENT(in) :: l_range     ! Optional partner index range

    ! partial_pot is the quantum spring energy of atom (i,k) in rik for polymer k
    ! with its neighbours (i,k-1) and (i,k+1), in the ring 1..p
    ! The coordinates in rik are not necessarily identical with those in r(:,i,k)
    ! The partner atom index is always the same as i
    ! The optional argument l_range restricts partner polymer ids to l=k-1, or l=k+1

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: l, l1, l2
    REAL               :: box_sq
    REAL               :: r_ik_il_sq
    REAL, DIMENSION(3) :: r_ik_il

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in spring_1'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in spring_1'
    END IF

    IF ( PRESENT ( l_range ) ) THEN
       SELECT CASE ( l_range )
       CASE ( lt ) ! Look at l = k-1 only
          l1 = k-1
          l2 = k-1
       CASE ( gt ) ! Look at l = k+1 only
          l1 = k+1
          l2 = k+1
       CASE default ! should never happen
          WRITE ( unit = error_unit, fmt='(a,i10)') 'l_range error ', l_range
          STOP 'Impossible error in spring_1'
       END SELECT
    ELSE ! Look at both l = k-1 and k+1
       l1 = k-1
       l2 = k+1
    END IF

    box_sq = box**2
    partial_pot = 0.0 ! Initialize

    DO l = l1, l2, 2 ! Loop over neighbours skipping l=k

       IF ( l < 1 ) THEN
          r_ik_il(:) = rik(:) - r(:,i,p) ! link to end of ring polymer
       ELSE IF ( l > p ) THEN
          r_ik_il(:) = rik(:) - r(:,i,1) ! link to start of ring polymer
       ELSE
          r_ik_il(:) = rik(:) - r(:,i,l)
       END IF

       r_ik_il(:)  = r_ik_il(:) - ANINT ( r_ik_il(:) ) ! Periodic boundaries in box=1 units
       r_ik_il_sq  = SUM ( r_ik_il**2 ) * box_sq       ! Squared distance in sigma=1 units
       partial_pot = partial_pot + r_ik_il_sq          ! Spring potential

    END DO ! End loop over possibilities skipping l=k

    ! Include numerical factors
    partial_pot = partial_pot * 0.5 * k_spring 

  END FUNCTION spring_1

  FUNCTION potential_lrc ( density, r_cut )
    IMPLICIT NONE
    REAL                :: potential_lrc ! Returns long-range correction to potential/atom
    REAL,    INTENT(in) :: density       ! Number density N/V
    REAL,    INTENT(in) :: r_cut         ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones potential per atom
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3           = 1.0 / r_cut**3
    potential_lrc = pi * ( (8.0/9.0)  * sr3**3  - (8.0/3.0)  * sr3 ) * density

  END FUNCTION potential_lrc

END MODULE qmc_module
