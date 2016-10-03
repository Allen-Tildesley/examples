! mc_chain_lj_module.f90
! Monte Carlo, single chain, LJ atoms
MODULE mc_module
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: regrow, energy, spring_pot
  PUBLIC :: n, r

  INTEGER                             :: n     ! number of atoms
  REAL,   DIMENSION(:,:), ALLOCATABLE :: r     ! atomic positions (3,n)
  REAL,   DIMENSION(:,:), ALLOCATABLE :: r_old ! working array (3,n)
  REAL,   DIMENSION(:,:), ALLOCATABLE :: r_new ! working array (3,n)

  REAL,    PARAMETER :: sigma = 1.0             ! LJ diameter
  REAL,    PARAMETER :: epslj = 1.0             ! LJ well depth
  INTEGER, PARAMETER :: lt = -1, gt = 1, ne = 0 ! options for j-range

CONTAINS

  SUBROUTINE introduction ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'LJ chain with harmonic bond potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ',   sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epslj = ', epslj    
  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output
    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'
  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays
    ALLOCATE ( r(3,n), r_old(3,n), r_new(3,n) ) 
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, r_old, r_new )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE regrow ( temperature, m_max, k_max, bond, k_spring, pot, ratio ) ! Routine to regrow polymer
    USE maths_module, ONLY : random_integer, random_vector

    REAL,    INTENT(in)    :: temperature ! Specified temperature
    INTEGER, INTENT(in)    :: m_max       ! Max atoms to regrow
    INTEGER, INTENT(in)    :: k_max       ! Number of random tries per atom in regrow
    REAL,    INTENT(in)    :: bond        ! Bond length
    REAL,    INTENT(in)    :: k_spring    ! Harmonic bond spring constant
    REAL,    INTENT(inout) :: pot         ! Total nonbonded potential energy (to be updated)
    REAL,    INTENT(out)   :: ratio       ! Acceptance ratio of moves

    ! A short sequence of m atoms (m<=m_max) is deleted and regrown in the CBMC manner
    ! We randomly select which end of the chain to apply each of these operations to
    ! At each stage, k_max different atom positions are tried
    ! Function random_bond selects bond lengths according to the internal (harmonic) potential
    ! Rosenbluth weights are computed using the external (nonbonded) potential
    ! Acceptance/rejection is determined using these weights
    
    REAL                       :: w_old   ! Old weight
    REAL                       :: w_new   ! New weight
    REAL                       :: pot_old ! Old potential energy
    REAL                       :: pot_new ! New potential energy
    INTEGER                    :: m       ! Number of atoms to regrow
    INTEGER                    :: c       ! Growth option
    INTEGER                    :: k       ! Try index
    REAL,   DIMENSION(k_max)   :: w       ! Rosenbluth weights (involving nonbonded interactions)
    REAL,   DIMENSION(3,k_max) :: r_try   ! Coordinates of trial atoms

    INTEGER :: i
    LOGICAL :: overlap
    REAL    :: d, d_max, std, zeta

    std   = SQRT(temperature/k_spring) ! Spring bond standard deviation
    d_max = 3.0*std                    ! Impose a limit on variation, say 3*std
    IF ( d_max > 0.5*bond ) THEN       ! must not be too large, say 0.5*bond
       WRITE ( unit=error_unit, fmt='(a,2f15.5)' ) 'Spring bond strength error', d_max, bond
       STOP 'Error in regrow'
    END IF
    d_max = d_max + bond ! This is the actual max d allowed

    ratio   = 0.0
    r_old   = r   ! Store old configuration
    pot_old = pot ! Store old nonbonded potential

    m = random_integer ( 1, m_max ) ! Number of atoms to regrow
    c = random_integer ( 1, 4 )     ! Growth option

    ! PART 1: CONSTRUCT NEW CONFIGURATION WITH NEW WEIGHT

    SELECT CASE ( c )

    CASE ( 1 ) ! Remove from end and add to end
       r(:,1:n-m) = r_old(:,1:n-m) ! copy first n-m atoms

    CASE ( 2 ) ! Remove from end and add to start
       r(:,1:n-m) = r_old(:,n-m:1:-1) ! copy and reverse first n-m atoms

    CASE ( 3 ) ! Remove from start and add to start
       r(:,1:n-m) = r_old(:,n:m+1:-1) ! copy and reverse last n-m atoms

    CASE ( 4 ) ! Remove from start and add to end
       r(:,1:n-m) = r_old(:,m+1:n) ! copy last n-m atoms

    END SELECT

    w_new = 1.0
    DO i = n-m+1, n ! Loop to regrow last m atoms, computing new weight

       DO k = 1, k_max ! Loop over k_max tries
          d          = random_bond ( bond, std, d_max )     ! Generate random bond length around d=bond
          r_try(:,k) = r(:,i-1) + d * random_vector()       ! Trial position in random direction from i-1
          CALL energy_1 ( r_try(:,k), i, lt, overlap, pot ) ! Nonbonded interactions with earlier atoms
          IF ( overlap ) THEN
             w(k) = 0.0 ! Weight for this try is zero
          ELSE
             w(k) = EXP(-pot/temperature) ! Weight for this try given by Boltzmann factor
          END IF
       END DO ! End loop over k_max tries

       k      = pick ( w )           ! Pick winning try according to weights
       r(:,i) = r_try(:,k)           ! Store winning position
       w_new  = w_new * REAL(SUM(w)) ! Accumulate total weight

    END DO ! End loop to regrow last m atoms, computing new weight

    CALL energy ( overlap, pot_new ) ! Compute full nonbonded energy
    r_new = r                        ! Store new configuration

    ! END OF PART 1: NEW CONFIGURATION AND WEIGHT ARE COMPLETE

    ! PART 2: RECONSTRUCT OLD CONFIGURATION WITH OLD WEIGHT

    SELECT CASE ( c )

    CASE ( 1, 2 ) ! Remove and hence reconstruct at end
       r = r_old(:,1:n) ! Copy all n atoms

    CASE ( 3, 4 ) ! Remove and reconstruct at start
       r = r_old(:,n:1:-1) ! Copy and reverse all n atoms

    END SELECT

    w_old = 1.0
    DO i = n-m+1, n ! Loop to regrow last m atoms computing old weight

       ! Old position and weight are stored as try 1

       r_try(:,1) = r(:,i)
       CALL energy_1 ( r_try(:,1), i, lt, overlap, pot ) ! Nonbonded energy with earlier atoms
       IF ( overlap ) THEN
          w(1) = 0.0 ! Current weight is zero; this should not happen
       ELSE
          w(1) = EXP(-pot/temperature) ! Current weight given by Boltzmann factor
       END IF

       ! Remaining tries only required to compute weight
       
       DO k = 2, k_max ! Loop over k_max-1 other tries
          d          = random_bond ( bond, std, d_max )     ! Generate random bond length around d=bond
          r_try(:,k) = r(:,i-1) + d * random_vector ( )     ! Trial position in random direction from i-1
          CALL energy_1 ( r_try(:,k), i, lt, overlap, pot ) ! Nonbonded interactions with earlier atoms
          IF ( overlap ) THEN
             w(k) = 0.0 ! Weight for this try is zero
          ELSE
             w(k) = EXP(-pot/temperature) ! Weight for this try given by Boltzmann factor
          END IF
       END DO ! End loop over k_max-1 other tries

       r(:,i) = r_try(:,1)           ! Restore winning position (always the original one)
       w_old  = w_old * REAL(SUM(w)) ! Accumulate total weight

    END DO ! End loop to regrow last m atoms computing old weight

    ! END OF PART 2: OLD CONFIGURATION AND WEIGHT ARE COMPLETE

    ! Choose either old or new configuration
    CALL RANDOM_NUMBER(zeta)
    IF ( zeta < ( w_new / w_old ) ) THEN ! accept
       ratio = 1.0
       r     = r_new
       pot   = pot_new
    ELSE ! reject
       r   = r_old
       pot = pot_old
    END IF

  END SUBROUTINE regrow

  FUNCTION pick ( w ) RESULT ( k )
    INTEGER                        :: k ! Returns one of the options with probability proportional to
    REAL, DIMENSION(:), INTENT(in) :: w ! the supplied weights

    REAL :: cumw, zeta

    CALL RANDOM_NUMBER ( zeta ) ! Random number between 0 and 1
    zeta = zeta*SUM(w)          ! Scale up to total weight
    k    = 1
    cumw = w(1)
    DO ! Loop over possible outcomes
       IF ( zeta <= cumw ) EXIT ! Random number less than cumulative weight up to k
       k = k+1
       IF ( k > SIZE(w) ) STOP 'Error in pick' ! Should never happen
       cumw = cumw+w(k)
    END DO ! End loop over possible outcomes

  END FUNCTION pick

  SUBROUTINE energy ( overlap, pot )
    LOGICAL, INTENT(out) :: overlap ! Shows if an overlap was detected
    REAL,    INTENT(out) :: pot     ! Nonbonded potential

    ! Calculates nonbonded potential for whole system
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

    REAL, DIMENSION(3), INTENT(in)  :: ri      ! Coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i       ! Index of atom of interest
    INTEGER,            INTENT(in)  :: j_range ! Partner index range
    LOGICAL,            INTENT(out) :: overlap ! Shows if an overlap was detected
    REAL,               INTENT(out) :: pot     ! Nonbonded potential

    ! Calculates potential energy of atom ri (not necessarily identical with r(:,i))
    ! due to nonbonded interactions with j/=i, j>i, or j<i depending on j_range
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the value of pot should not be used
    ! Coordinates are assumed to be in LJ units where sigma = 1, no periodic boundaries
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: j, j1, j2
    REAL               :: sr2, sr6, rij_sq
    REAL, DIMENSION(3) :: rij
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

    IF ( n > SIZE(r,dim=2) ) THEN ! Should never happen
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

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( ABS(j-i) <= 1 ) CYCLE ! Skip self and bonded neighbours

       rij(:) = ri(:) - r(:,j) ! Separation vector
       rij_sq = SUM ( rij**2 ) ! Squared separation
       sr2    = 1.0 / rij_sq   ! (sigma/rij)**2

       IF ( sr2 > sr2_overlap ) THEN
          overlap = .TRUE.
          EXIT ! Jump out of loop
       END IF

       sr6 = sr2**3
       pot = pot + sr6**2 - sr6

    END DO ! End loop over selected range of partners

    pot = 4.0 * pot ! Numerical factor

  END SUBROUTINE energy_1

  FUNCTION random_bond ( b, std, d_max ) RESULT ( d )
    USE maths_module, ONLY : random_normal
    REAL             :: d     ! Returns random bond length, generated around a
    REAL, INTENT(in) :: b     ! specified mean bond length, Gaussian-distributed with
    REAL, INTENT(in) :: std   ! specified standard deviation, but subject to a
    REAL, INTENT(in) :: d_max ! maximum allowed value of d

    ! Uses von Neumann's rejection method to sample (d**2)*exp(-0.5*(d-b)**2/std**2)
    ! The sampled distribution is the same but with d in the prefactor (which arises from
    ! the Jacobian in 3D) replaced by the constant d_max, which makes it a Gaussian function
    ! Hence, the range must be restricted to d<d_max, for the rejection method to work
    ! It will be reasonably efficient provided std is small compared with b
    ! This is essentially the same method as an example in
    ! Understanding Molecular Simulation by D Frenkel and B Smit

    REAL :: zeta

    DO ! Loop over attempts
       d = random_normal ( b, std )        ! Generate usual Gaussian distribution
       IF ( d < 0.0 .OR. d > d_max ) CYCLE ! Reject if outside range
       CALL RANDOM_NUMBER ( zeta )         ! Uniform in range 0..1 for rejection method
       IF ( zeta <= (d/d_max)**2 ) EXIT    ! Compare with ratio of distributions to accept
    END DO ! End loop over attempts
  END FUNCTION random_bond

  FUNCTION spring_pot ( b, k ) RESULT ( pot )
    REAL             :: pot ! Internal spring potential energy, given
    REAL, INTENT(in) :: b   ! specified bond length and
    REAL, INTENT(in) :: k   ! specified force constant

    INTEGER            :: i
    REAL               :: d
    REAL, DIMENSION(3) :: rij

    pot = 0.0
    
    DO i = 1, n-1 ! Loop over atoms
       rij = r(:,i) - r(:,i+1)     ! Bond vector
       d = SQRT ( SUM ( rij**2 ) ) ! Bond distance
       pot = pot + (d-b)**2        ! Squared displacement
    END DO ! End loop over atoms

    ! Numerical factors
    pot = pot * 0.5 * k
    
  END FUNCTION spring_pot
  
END MODULE mc_module

