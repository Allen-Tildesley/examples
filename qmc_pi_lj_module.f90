! qmc_pi_lj_module.f90
! Energy and move routines for PIMC simulation, LJ potential
MODULE qmc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, p, r
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: energy_cl_1, energy_cl, energy_qu_1, energy_qu
  PUBLIC :: move
  PUBLIC :: pot_type

  INTEGER                                :: n ! number of atoms
  INTEGER                                :: p ! number of beads
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: r ! positions (3,n,p)

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! j-range and l-range options
  REAL,    PARAMETER :: sigma = 1.0     ! Lennard-Jones diameter (unit of length)
  REAL,    PARAMETER :: epslj = 1.0     ! Lennard-Jones well depth (unit of energy)

  TYPE pot_type ! A composite variable for interaction energies comprising
     REAL    :: pot ! the potential energy and
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
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_pot_type

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

    ALLOCATE ( r(3,n,p) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r )
  END SUBROUTINE deallocate_arrays

  FUNCTION energy_cl ( box, r_cut ) RESULT ( energy )
    TYPE(pot_type)     :: energy ! Returns a composite of pot and ovr
    REAL, INTENT(in) :: box    ! Simulation box length
    REAL, INTENT(in) :: r_cut  ! Potential cutoff distance

    ! energy%pot is the nonbonded classical potential energy for whole system
    ! energy%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the value of energy%pot should not be used
    ! Actual calculation is performed by subroutine energy_cl_1

    TYPE(pot_type) :: energy_ik
    INTEGER      :: i, k

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in energy_cl'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in energy_cl'
    END IF

    energy = pot_type ( pot=0.0, ovr=.FALSE. ) ! Initialize

    DO k = 1, p ! Loop over ring polymers
       DO i = 1, n - 1 ! Loop over atoms within polymer
          energy_ik = energy_cl_1 ( r(:,i,k), i, k, box, r_cut, gt )
          IF ( energy_ik%ovr ) THEN
             energy%ovr = .TRUE. ! Overlap detected
             RETURN              ! Return immediately
          END IF
          energy = energy + energy_ik
       END DO ! End loop over atoms within polymer
    END DO ! End loop over ring polymers

    energy%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION energy_cl

  FUNCTION energy_cl_1 ( rik, i, k, box, r_cut, j_range ) RESULT ( energy )
    TYPE(pot_type)                    :: energy  ! Returns a composite of pot and ovr
    REAL, DIMENSION(3), INTENT(in)  :: rik     ! Coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i, k    ! Index and polymer id of atom of interest
    REAL,               INTENT(in)  :: box     ! Simulation box length
    REAL,               INTENT(in)  :: r_cut   ! Potential cutoff distance
    INTEGER, OPTIONAL,  INTENT(in)  :: j_range ! Optional partner index range

    ! energy%pot is the nonbonded potential energy of atom rik with a set of other atoms
    ! energy%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the value of energy%pot should not be used
    ! The coordinates in rik are not necessarily identical with those in r(:,i,k)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: j, j1, j2
    REAL               :: r_cut_box, r_cut_box_sq, box_sq
    REAL               :: sr2, sr6, r_ik_jk_sq
    REAL, DIMENSION(3) :: r_ik_jk
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in energy_cl_1'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in energy_cl_1'
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
          STOP 'Impossible error in energy_cl_1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box**2
    box_sq       = box**2

    energy = pot_type ( pot=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2

       IF ( i == j ) CYCLE

       r_ik_jk(:) = rik(:) - r(:,j,k)
       r_ik_jk(:) = r_ik_jk(:) - ANINT ( r_ik_jk(:) ) ! periodic boundaries in box=1 units
       r_ik_jk_sq = SUM ( r_ik_jk**2 )

       IF ( r_ik_jk_sq < r_cut_box_sq ) THEN

          r_ik_jk_sq = r_ik_jk_sq * box_sq ! now in sigma=1 units
          sr2        = 1.0 / r_ik_jk_sq    ! (sigma/rikjk)**2

          IF ( sr2 > sr2_overlap ) THEN
             energy%ovr = .TRUE. ! Overlap detected
             RETURN              ! Return immediately
          END IF

          sr6 = sr2**3
          energy%pot = energy%pot + sr6**2 - sr6

       END IF

    END DO

    ! Include numerical factors
    energy%pot = 4.0 * energy%pot        ! factor of 4*epsilon
    energy%pot = energy%pot / REAL ( p ) ! classical potentials are weaker by a factor p

    energy%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION energy_cl_1

  FUNCTION energy_qu ( box, k_spring ) RESULT ( pot )
    REAL             :: pot      ! Returns quantum spring potential 
    REAL, INTENT(in) :: box      ! Simulation box length
    REAL, INTENT(in) :: k_spring ! Spring potential constant

    ! Calculates quantum spring potential for whole system
    ! Actual calculation is performed by subroutine energy_qu_1

    INTEGER :: i, k

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in energy_qu'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in energy_qu'
    END IF

    pot = 0.0

    DO k = 1, p ! Loop over ring polymers
       DO i = 1, n ! Loop over atoms within polymer
          pot = pot + energy_qu_1 ( r(:,i,k), i, k, box, k_spring, gt )
       END DO ! End loop over atoms within polymer
    END DO ! End loop over ring polymers

  END FUNCTION energy_qu

  FUNCTION energy_qu_1 ( rik, i, k, box, k_spring, l_range ) RESULT ( pot )
    REAL                           :: pot      ! Returns quantum potential
    REAL, DIMENSION(3), INTENT(in) :: rik      ! Coordinates of atom of interest
    INTEGER,            INTENT(in) :: i, k     ! Index and polymer id of atom of interest
    REAL,               INTENT(in) :: box      ! Simulation box length
    REAL,               INTENT(in) :: k_spring ! Spring potential constant
    INTEGER, OPTIONAL,  INTENT(in) :: l_range  ! Optional partner index range

    ! Calculates quantum spring potential energy of atom in rik for given polymer k
    ! The coordinates in rik are not necessarily identical with those in r(:,i,k)
    ! The optional argument l_range restricts partner indices to l=k-1, or l=k+1

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: l, l1, l2
    REAL               :: box_sq
    REAL               :: r_ik_il_sq
    REAL, DIMENSION(3) :: r_ik_il

    IF ( n /= SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, n=', n, SIZE(r,dim=2)
       STOP 'Error in energy_qu_1'
    END IF
    IF ( p /= SIZE(r,dim=3) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r, p=', p, SIZE(r,dim=3)
       STOP 'Error in energy_qu_1'
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
          STOP 'Impossible error in energy_qu_1'
       END SELECT
    ELSE ! Look at both l = k-1 and k+1
       l1 = k-1
       l2 = k+1
    END IF

    box_sq = box**2
    pot    = 0.0

    DO l = l1, l2, 2 ! Loop over possibilities skipping l=k

       IF ( l < 1 ) THEN
          r_ik_il(:) = rik(:) - r(:,i,p) ! link around to top
       ELSE IF ( l > p ) THEN
          r_ik_il(:) = rik(:) - r(:,i,1) ! link around to bottom
       ELSE
          r_ik_il(:) = rik(:) - r(:,i,l)
       END IF

       r_ik_il(:) = r_ik_il(:) - ANINT ( r_ik_il(:) ) ! Periodic boundaries in box=1 units
       r_ik_il_sq = SUM ( r_ik_il**2 ) * box_sq       ! Squared distance in sigma=1 units
       pot        = pot + 0.5 * k_spring * r_ik_il_sq ! Spring potential

    END DO ! End loop over possibilities skipping l=k

  END FUNCTION energy_qu_1

  SUBROUTINE move ( i, k, rik )
    INTEGER,               INTENT(in) :: i, k
    REAL,    DIMENSION(3), INTENT(in) :: rik

    r(:,i,k) = rik

  END SUBROUTINE move

END MODULE qmc_module
