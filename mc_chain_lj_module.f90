! mc_chain_lj_module.f90
! Monte Carlo, single chain, LJ atoms
MODULE mc_module
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: model_description, allocate_arrays, deallocate_arrays
  PUBLIC :: regrow, energy
  PUBLIC :: bond, n, r

  INTEGER                             :: n     ! number of atoms
  REAL                                :: bond  ! bond length
  REAL,   DIMENSION(:,:), ALLOCATABLE :: r     ! atomic positions (3,n)
  REAL,   DIMENSION(:,:), ALLOCATABLE :: r_old ! working array (3,n)
  REAL,   DIMENSION(:,:), ALLOCATABLE :: r_new ! working array (3,n)

  REAL,    PARAMETER :: sigma = 1.0             ! LJ diameter
  REAL,    PARAMETER :: epslj = 1.0             ! LJ well depth
  INTEGER, PARAMETER :: lt = -1, gt = 1, ne = 0 ! options for j-range

CONTAINS

  SUBROUTINE model_description ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'LJ chain with harmonic bond length'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ', sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epslj = ', epslj    
  END SUBROUTINE model_description

  SUBROUTINE allocate_arrays
    ALLOCATE ( r(3,n), r_old(3,n), r_new(3,n) ) 
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, r_old, r_new )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE regrow ( temperature, m_max, k_max, k_spring, pot, ratio ) ! regrow polymer
    USE maths_module, ONLY : random_integer, random_orientation_vector

    REAL,    INTENT(in)    :: temperature ! temperature
    INTEGER, INTENT(in)    :: m_max       ! max atoms to regrow
    INTEGER, INTENT(in)    :: k_max       ! number of random tries per atom in regrow
    REAL,    INTENT(in)    :: k_spring    ! harmonic bond spring constant
    REAL,    INTENT(inout) :: pot         ! energy (to be updated)
    REAL,    INTENT(out)   :: ratio       ! acceptance ratio of moves

    REAL                         :: w_old   ! old weight
    REAL                         :: w_new   ! new weight
    REAL                         :: pot_old ! old potential energy
    REAL                         :: pot_new ! new potential energy
    INTEGER                      :: m       ! number of atoms to regrow
    INTEGER                      :: c       ! growth option
    INTEGER                      :: k       ! try index
    REAL,    DIMENSION(k_max)    :: w       ! Rosenbluth weights (involving non-bonded interactions)
    REAL,    DIMENSION(3,k_max)  :: r_try   ! coordinates of trial atoms

    INTEGER              :: i
    LOGICAL              :: overlap
    REAL                 :: d, d_range, stddev, zeta
    REAL,   DIMENSION(3) :: u

    stddev  = SQRT(temperature / k_spring) ! spring bond standard deviation
    d_range = 3.0*stddev                   ! impose a limit on variation
    IF ( d_range > bond ) THEN             ! must not be too large
       WRITE ( unit=error_unit, fmt='(a,2f15.5)' ) 'Spring bond strength error', d_range, bond
       STOP 'Error in regrow'
    END IF

    ratio   = 0.0
    r_old   = r   ! store old configuration
    pot_old = pot ! store old potential

    m = random_integer ( 1, m_max ) ! number of atoms to regrow
    c = random_integer ( 1, 4 )     ! growth option

    ! PART 1: CONSTRUCT NEW CONFIGURATION WITH NEW WEIGHT

    SELECT CASE ( c )

    CASE ( 1 ) ! remove from end and add to end
       r(:,1:n-m) = r_old(:,1:n-m) ! copy first n-m atoms

    CASE ( 2 ) ! remove from end and add to start
       r(:,1:n-m) = r_old(:,n-m:1:-1) ! copy and reverse first n-m atoms

    CASE ( 3 ) ! remove from start and add to start
       r(:,1:n-m) = r_old(:,n:m+1:-1) ! copy and reverse last n-m atoms

    CASE ( 4 ) ! remove from start and add to end
       r(:,1:n-m) = r_old(:,m+1:n) ! copy last n-m atoms

    END SELECT

    w_new = 1.0
    DO i = n-m+1, n ! Loop to regrow last m atoms, computing new weight

       DO k = 1, k_max ! loop over k_max tries
          CALL random_orientation_vector ( u )
          d = random_bond ( bond, stddev, d_range )         ! generate random bond length around d=bond
          r_try(:,k) = r(:,i-1) + d * u                     ! generate trial position
          CALL energy_1 ( r_try(:,k), i, lt, overlap, pot ) ! non-bonded interactions with earlier atoms
          IF ( overlap ) THEN
             w(k) = 0.0
          ELSE
             w(k) = EXP(-pot/temperature)  ! store overlap weight for this try
          END IF
       END DO ! end loop over k_max tries

       k      = pick(w)              ! pick winning position
       r(:,i) = r_try(:,k)           ! store winning position
       w_new  = w_new * REAL(SUM(w)) ! accumulate total weight

    END DO ! end loop to regrow last m atoms, computing new weight

    ! New configuration and weight are complete
    CALL energy ( overlap, pot_new ) ! compute full energy
    r_new = r

    ! PART 2: RECONSTRUCT OLD CONFIGURATION WITH OLD WEIGHT

    SELECT CASE ( c )

    CASE ( 1, 2 ) ! remove and hence reconstruct at end
       r = r_old(:,1:n) ! copy all n atoms

    CASE ( 3, 4 ) ! remove and reconstruct at start
       r = r_old(:,n:1:-1) ! copy and reverse all n atoms

    END SELECT

    w_old = 1.0
    DO i = n-m+1, n ! Loop to regrow last m atoms calculating old weight

       r_try(:,1) = r(:,i) ! store old position in try 1
       CALL energy_1 ( r_try(:,1), i, lt, overlap, pot ) ! nonbonded energy with earlier atoms
       IF ( overlap ) THEN
          w(1) = 0.0 ! this should not happen
       ELSE
          w(1) = EXP(-pot/temperature) ! store overlap weight in try 1
       END IF

       DO k = 2, k_max ! loop over k_max-1 other tries
          CALL random_orientation_vector ( u )
          d = random_bond ( bond, stddev, d_range )         ! generate random bond length around d=bond
          r_try(:,k) = r(:,i-1) + d * u                     ! generate trial position
          CALL energy_1 ( r_try(:,k), i, lt, overlap, pot ) ! nonbonded energy with earlier atoms
          IF ( overlap ) THEN
             w(k) = 0.0
          ELSE
             w(k) = EXP(-pot/temperature) ! store overlap weight
          END IF
       END DO ! end loop over k_max-1 other tries

       r(:,i) = r_try(:,1)           ! restore winning position (always the original one)
       w_old  = w_old * REAL(SUM(w)) ! accumulate total weight

    END DO ! end loop to regrow last m atoms calculating old weight

    ! Choose either old or new configuration
    CALL RANDOM_NUMBER(zeta)
    IF ( zeta < (w_new/w_old) ) THEN ! accept
       ratio = 1.0
       r     = r_new
       pot   = pot_new
    ELSE ! reject
       r   = r_old
       pot = pot_old
    END IF

  END SUBROUTINE regrow

  FUNCTION pick ( w ) RESULT ( k ) !  select one of the tries according to overlap weight
    IMPLICIT NONE

    ! Argument
    REAL, DIMENSION(:), INTENT(in) :: w
    INTEGER                        :: k ! Result

    ! Local variables
    REAL :: cumw, ran

    CALL RANDOM_NUMBER(ran)
    ran = ran*SUM(w)
    k=1
    cumw = w(1)
    DO WHILE ( cumw < ran )
       k = k+1
       cumw = cumw+w(k)
    ENDDO

  END FUNCTION pick

  SUBROUTINE energy ( overlap, pot )
    LOGICAL, INTENT(out) :: overlap ! shows if an overlap was detected
    REAL,    INTENT(out) :: pot     ! potential

    ! Calculates potential for whole system
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the value of pot should not be used
    ! Actual calculation is performed by subroutine energy_1

    REAL               :: pot_i
    INTEGER            :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy'
    END IF

    overlap  = .FALSE.
    pot      = 0.0

    DO i = 1, n - 1
       CALL energy_1 ( r(:,i), i, gt, overlap, pot_i )
       IF ( overlap ) EXIT ! jump out of loop
       pot  = pot  + pot_i
    END DO

  END SUBROUTINE energy

  SUBROUTINE energy_1 ( ri, i, j_range, overlap, pot )

    REAL, DIMENSION(3), INTENT(in)  :: ri         ! coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i, j_range ! index, and partner index range
    LOGICAL,            INTENT(out) :: overlap    ! shows if an overlap was detected
    REAL,               INTENT(out) :: pot        ! potential

    ! Calculates potential energy and virial of atom in ri
    ! with j/=i, j>i, or j<i depending on j_range
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the value of pot should not be used
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: j, j1, j2
    REAL               :: sr2, sr6, rij_sq
    REAL, DIMENSION(3) :: rij
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy_1'
    END IF

    pot     = 0.0
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

       IF ( ABS(j-i) <= 1 ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij_sq = SUM ( rij**2 )
       sr2    = 1.0 / rij_sq    ! (sigma/rij)**2

       IF ( sr2 > sr2_overlap ) THEN
          overlap = .TRUE.
          EXIT ! jump out of loop
       END IF

       sr6 = sr2**3
       pot = pot + sr6**2 - sr6

    END DO
    pot = 4.0 * pot

  END SUBROUTINE energy_1

  FUNCTION random_bond ( bond, stddev, d_range ) RESULT ( d ) ! generate random bond length around d=bond
    USE maths_module, ONLY : random_normal
    REAL, INTENT(in) :: bond    ! reference bond length
    REAL, INTENT(in) :: stddev  ! standard deviation in Gaussian
    REAL, INTENT(in) :: d_range ! max allowed range of d
    REAL             :: d       ! result

    REAL :: zeta, d_max

    ! Uses von Neumann's rejection method to sample (d**2)*exp(-0.5*(d-bond)**2/stddev**2)
    ! The sampled distribution is the same but with d replaced by the constant d_max
    ! Hence, the range must be restricted to d<d_max, for the rejection method to work
    ! and we actually restrict d symmetrically about d=bond
    ! It will be reasonably efficient provided stddev is small compared with bond
    ! This is essentially the same method as an example in
    ! Understanding Molecular Simulation by D Frenkel and B Smit

    d_max = bond + d_range
    DO
       d = random_normal ( 0.0, stddev )
       IF ( ABS(d) > d_range ) CYCLE ! reject outside range
       d = d + bond
       CALL RANDOM_NUMBER ( zeta )
       IF ( zeta <= (d/d_max)**2 ) EXIT ! accept
    END DO
  END FUNCTION random_bond

END MODULE mc_module

