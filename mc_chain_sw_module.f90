! mc_chain_sw_module.f90
! Monte Carlo, single chain, square wells
MODULE mc_module
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: regrow, cranks, pivots, weight, qcount, write_histogram, histogram_flat, zero_histogram
  PUBLIC :: n, nq, s, r, ds, wl

  INTEGER                              :: n       ! Number of atoms
  INTEGER                              :: nq      ! Maximum anticipated energy
  LOGICAL                              :: wl      ! Flag for Wang-Landau method
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r       ! Atomic positions (3,n)
  REAL,    DIMENSION(:),   ALLOCATABLE :: h       ! Histogram of q values (0:nq)
  REAL,    DIMENSION(:),   ALLOCATABLE :: s       ! "Entropy" histogram used in acceptance (0:nq)
  LOGICAL, DIMENSION(:),   ALLOCATABLE :: hit     ! Identifies q-values that have been visited (0:nq)
  REAL                                 :: ds      ! Entropy increment for Wang-Landau method

  REAL, DIMENSION(:,:), ALLOCATABLE :: r_old, r_new ! Working arrays (3,n)

  REAL,    PARAMETER :: sigma = 1.0             ! core diameter
  INTEGER, PARAMETER :: lt = -1, gt = 1, ne = 0 ! options for j-range

CONTAINS

  SUBROUTINE introduction ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Hard-sphere chain with fixed bond length'
    WRITE ( unit=output_unit, fmt='(a)'           ) 'Square-well attractive potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ', sigma    
  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output
    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'
  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays
    nq = 6*n ! Anticipated maximum number of pair interactions within range
    ALLOCATE ( r(3,n), r_old(3,n), r_new(3,n) ) 
    ALLOCATE ( h(0:nq), s(0:nq), hit(0:nq) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, r_old, r_new, h, s, hit )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE regrow ( m_max, k_max, bond, range, q, ratio ) ! Routine to regrow polymer
    USE maths_module, ONLY : random_integer, random_vector, pick

    INTEGER, INTENT(in)    :: m_max  ! Max atoms to regrow
    INTEGER, INTENT(in)    :: k_max  ! Number of random tries per atom in regrow
    REAL,    INTENT(in)    :: bond   ! Bond length
    REAL,    INTENT(in)    :: range  ! Range of attractive well
    INTEGER, INTENT(inout) :: q      ! Energy, negative of (to be updated)
    REAL,    INTENT(out)   :: ratio  ! Acceptance ratio of moves

    ! This routine carries out a single regrowth move
    ! A short sequence of m atoms (m<=m_max) is deleted and regrown in the CBMC manner
    ! We randomly select which end of the chain to apply each of these operations to
    ! At each stage, k_max different atom positions are tried
    ! The bond length is fixed throughout
    ! Weights used in the regrowth are athermal, computed only on the basis of the
    ! hard-core overlap part of the non-bonded interactions: essentially they count non-overlaps
    ! Hence they are suitable for use in both NVT and Wang-Landau simulations
    
    INTEGER                      :: q_old, q_new ! Old and new energies
    REAL                         :: w_old, w_new ! Old and new weights
    INTEGER                      :: m            ! Number of atoms to regrow
    INTEGER                      :: c            ! Growth option
    INTEGER                      :: k            ! Try index
    INTEGER, DIMENSION(k_max)    :: w            ! Rosenbluth weights (here integers, 0 or 1)
    REAL,    DIMENSION(3,k_max)  :: r_try        ! Coordinates of trial atoms
    INTEGER                      :: i            ! Atom index

    REAL, PARAMETER :: w_tol = 0.5 ! Min weight tolerance (w takes integer values)

    ratio = 0.0
    r_old = r ! Store old configuration
    q_old = q ! Store old q

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
          r_try(:,k) = r(:,i-1) + bond * random_vector() ! Trial position in random direction from i-1
          w(k)       = weight_1 ( r_try(:,k), i, lt )    ! Store overlap weight for this try
       END DO ! End loop over k_max tries

       IF ( SUM(w) == 0 ) THEN ! Early exit if this happens at any stage
          r = r_old                   ! Restore original configuration
          q = q_old                   ! Redundant, but for clarity
          CALL update_histogram ( q ) ! Always do this
          RETURN                      ! No possible move: reject
       END IF

       k      = pick ( w )           ! Pick winning try according to weights
       r(:,i) = r_try(:,k)           ! Store winning position
       w_new  = w_new * REAL(SUM(w)) ! Accumulate total weight

    END DO ! End loop to regrow last m atoms, computing new weight

    IF ( w_new < w_tol ) THEN ! This should be a redundant test
       r = r_old                   ! Restore original configuration
       q = q_old                   ! Redundant, but for clarity
       CALL update_histogram ( q ) ! Always do this
       RETURN                      ! No possible move: reject
    END IF

    q_new = qcount ( range ) ! Compute full new nonbonded energy
    r_new = r                ! Store new configuration

    ! END OF PART 1: NEW CONFIGURATION AND WEIGHT ARE COMPLETE

    ! PART 2: RECONSTRUCT OLD CONFIGURATION WITH OLD WEIGHT

    SELECT CASE ( c )

    CASE ( 1, 2 ) ! Remove and hence reconstruct at end
       r = r_old(:,1:n) ! copy all n atoms

    CASE ( 3, 4 ) ! Remove and reconstruct at start
       r = r_old(:,n:1:-1) ! copy and reverse all n atoms

    END SELECT

    w_old = 1.0
    DO i = n-m+1, n ! Loop to regrow last m atoms calculating old weight

       ! Old position and weight are stored as try 1

       r_try(:,1) = r(:,i) ! Current position
       w(1) = 1            ! Current weight must be 1

       ! Remaining tries only required to compute weight

       DO k = 2, k_max ! Loop over k_max-1 other tries
          r_try(:,k) = r(:,i-1) + bond * random_vector() ! Trial position in random direction from i-1
          w(k)       = weight_1 ( r_try(:,k), i, lt )    ! Store overlap weight for this try
       END DO ! End loop over k_max-1 other tries

       r(:,i) = r_try(:,1)           ! Restore winning position (always the original one)
       w_old  = w_old * REAL(SUM(w)) ! Accumulate total weight

    END DO ! End loop to regrow last m atoms calculating old weight

    ! END OF PART 2: OLD CONFIGURATION AND WEIGHT ARE COMPLETE

    ! Choose either old or new configuration
    IF ( accept ( q_old, q_new, w_old, w_new ) ) THEN ! Accept
       ratio = 1.0
       r     = r_new
       q     = q_new
    ELSE ! Reject
       r = r_old
       q = q_old
    END IF

    CALL update_histogram ( q ) ! Always do this

  END SUBROUTINE regrow

  SUBROUTINE pivots ( n_try, phi_max, range, q, ratio )
    USE maths_module, ONLY : random_integer, random_vector, rotate_vector

    INTEGER, INTENT(in)    :: n_try   ! Number of pivots to try
    REAL,    INTENT(in)    :: phi_max ! Maximum angle of pivot
    REAL,    INTENT(in)    :: range   ! Range of attractive well
    INTEGER, INTENT(inout) :: q       ! Energy, negative of (to be updated)
    REAL,    INTENT(out)   :: ratio   ! Move acceptance ratio

    ! This routine carries out a specified number of pivot moves
    ! An atom is picked at random, and the part of the chain lying to one side of it
    ! is rotated as a whole, by a random angle, about a randomly oriented axis
    ! There are no weights to take into account in the acceptance/rejection decision
    ! (the function weight_1 is used simply to indicate overlap / no overlap)
    
    INTEGER            :: q_old, q_new ! Old and new energies
    INTEGER            :: j            ! Pivot position
    INTEGER            :: k            ! Which part to pivot
    INTEGER            :: i            ! Atom index
    INTEGER            :: try, n_acc   ! Trial and acceptance counters
    REAL, DIMENSION(3) :: rj, rij      ! Local position vectors
    REAL, DIMENSION(3) :: axis         ! Pivot axis
    REAL               :: phi          ! Pivot angle
    REAL               :: zeta         ! Random number

    n_acc = 0 ! Zero acceptance counter

    tries: DO try = 1, n_try ! Loop over pivot move attempts

       r_old = r ! Store old configuration
       q_old = q ! Store old q

       j = random_integer ( 2, n-1 ) ! Pivot position (not at either end)
       k = random_integer ( 1, 2 )   ! Which part to pivot (actually redundant here)

       axis = random_vector ()            ! Pivot rotation axis
       CALL RANDOM_NUMBER ( zeta )        ! Uniform random number in range (0,1)
       phi = phi_max * ( 2.0*zeta - 1.0 ) ! Pivot rotation angle in desired range

       SELECT CASE ( k )
       CASE ( 1 )
          r = r_old ! Copy atoms
       CASE ( 2 )
          r = r_old(:,n:1:-1) ! Copy and reverse atoms
       END SELECT

       ! Pivot, and check for overlap in new configuration
       ! NB include overlap checks within rotated segment, because of roundoff

       rj = r(:,j) ! Pivot point

       DO i = j+1, n ! Loop over moving atoms
          rij = r(:,i) - rj                      ! Relative vector of atom
          rij = rotate_vector ( phi, axis, rij ) ! Rotate relative vector
          r(:,i) = rj + rij                      ! New absolute position

          IF ( weight_1 ( r(:,i), i, lt ) == 0 ) THEN ! Check for overlaps
             r = r_old                    ! Restore original configuration
             q = q_old                    ! Redundant but for clarity
             CALL update_histogram ( q )  ! Always do this
             CYCLE tries                  ! This move is rejected
          END IF ! End check for overlaps

       END DO ! End loop over moving atoms

       q_new = qcount ( range ) ! Compute full new nonbonded energy
       r_new = r                ! Store new configuration

       IF ( accept ( q_old, q_new ) ) THEN ! Accept
          n_acc = n_acc + 1
          r     = r_new
          q     = q_new
       ELSE ! Reject
          r = r_old
          q = q_old
       END IF

       CALL update_histogram ( q ) ! Always do this

    END DO tries ! End loop over pivot move attempts

    ratio = REAL(n_acc) / REAL(n_try)

  END SUBROUTINE pivots

  SUBROUTINE cranks ( n_try, phi_max, range, q, ratio )
    USE maths_module, ONLY : random_integer, rotate_vector

    INTEGER, INTENT(in)    :: n_try   ! Number of crankshaft moves to try
    REAL,    INTENT(in)    :: phi_max ! Maximum crank angle
    REAL,    INTENT(in)    :: range   ! Range of attractive well
    INTEGER, INTENT(inout) :: q       ! Energy, negative of (to be updated)
    REAL,    INTENT(out)   :: ratio   ! Move acceptance ratio

    ! This routine carries out a specified number of crankshaft moves
    ! An atom is picked at random. Unless it is an end-atom, a rotation axis
    ! is defined as the line joining the two atoms on either side, and the
    ! chosen atom is rotated about that axis by a random angle.
    ! In the case of end-atoms, the axis is defined by the line joining its
    ! nearest and next-nearest neighbours.
    ! There are no weights to take into account in the acceptance/rejection decision
    ! (the function weight_1 is used simply to indicate overlap / no overlap)

    INTEGER            :: q_old, q_new ! Old and new energies
    INTEGER            :: i            ! Atom index
    INTEGER            :: try, n_acc   ! Trial and acceptance counters
    REAL, DIMENSION(3) :: rj, rij      ! Local position vectors
    REAL, DIMENSION(3) :: axis         ! Crankshaft axis
    REAL               :: phi          ! Crankshaft angle
    REAL               :: zeta         ! Random number
    REAL               :: norm

    n_acc = 0 ! Zero acceptance counter

    tries: DO try = 1, n_try ! Loop over crankshaft move attempts

       r_old = r ! Store old configuration
       q_old = q ! Store old q

       r = r_old

       i     = random_integer ( 1, n )                   ! Pick random atom to move
       q_new = q_old - qcount_1 ( r(:,i), i, ne, range ) ! Subtract old energies for moving atom

       CALL RANDOM_NUMBER ( zeta )        ! Uniform random number in range (0,1)
       phi = ( 2.0*zeta - 1.0 ) * phi_max ! Rotation angle in desired range

       IF ( i == 1 ) THEN            ! Rotate about 2--3 bond
          rj   = r(:,2)              !   reference position
          axis = r(:,2) - r(:,3)     !   axis of rotation
       ELSE IF ( i == n ) THEN       ! Rotate about (n-1)--(n-2) bond
          rj   = r(:,n-1)            !   reference position
          axis = r(:,n-1) - r(:,n-2) !   axis of rotation
       ELSE                          ! Rotate about (i-1)--(i+1) bond
          rj   = r(:,i-1)            !   reference position
          axis = r(:,i+1) - r(:,i-1) !   axis of rotation
       END IF

       norm = SQRT(SUM(axis**2)) ! Squared length of rotation axis
       axis = axis / norm        ! Unit vector along rotation axis

       rij = r(:,i) - rj                      ! Relative vector of atom
       rij = rotate_vector ( phi, axis, rij ) ! Rotate relative vector
       r(:,i) = rj + rij                      ! New absolute position

       IF ( weight_1 ( r(:,i), i, ne ) == 0 ) THEN ! Check for overlaps
          r = r_old                   ! Restore original configuration
          q = q_old                   ! Redundant but for clarity
          CALL update_histogram ( q ) ! Always do this
          CYCLE tries                 ! This move is rejected
       END IF ! End check for overlaps

       q_new = q_new + qcount_1 ( r(:,i), i, ne, range ) ! Add new energies for moving atom
       r_new = r

       IF ( accept ( q_old, q_new ) ) THEN ! Accept
          n_acc = n_acc + 1
          r     = r_new
          q     = q_new
       ELSE ! Reject
          r = r_old
          q = q_old
       END IF

       CALL update_histogram ( q ) ! Always do this

    END DO tries ! End loop over crankshaft move attempts

    ratio = REAL(n_acc) / REAL(n_try)

  END SUBROUTINE cranks

  FUNCTION weight () RESULT ( w )
    INTEGER :: w ! Returns configuration weight = 0 (overlap) or 1 (no overlap)

    ! Arithmetically w is the product of all the pair Boltzmann factors
    ! Here, each is 0 or 1, but in the more general case we would multiply them
    
    INTEGER :: i

    w = 0 ! Weight is zero unless we get through the following loop

    DO i = 1, n-1 ! Loop over atoms checking upwards for weight
       IF ( weight_1 ( r(:,i), i, gt ) == 0 ) RETURN ! Overlap detected, no need to continue
    END DO ! End loop over atoms checking upwards for weight

    w = 1 ! No overlaps detected, so weight is one

  END FUNCTION weight

  FUNCTION weight_1 ( ri, i, j_range ) RESULT ( w )
    INTEGER                           :: w       ! Returns weight (0 or 1) associated with single atom
    REAL,    DIMENSION(3), INTENT(in) :: ri      ! Coordinates of atom of interest
    INTEGER,               INTENT(in) :: i       ! Index of atom of interest
    INTEGER,               INTENT(in) :: j_range ! Partner index range

    ! Calculates weight associated with atom ri (not necessarily identical with r(:,i))
    ! due to nonbonded interactions with j/=i, j>i, or j<i depending on j_range
    ! This is athermal, based on hard core nonbonded overlaps only
    ! The result is either 0 (overlap) or 1 (non overlap)
    ! Arithmetically w is the product of all the pair Boltzmann factors
    ! Here, each is 0 or 1, but in the more general case we would multiply them

    INTEGER            :: j, j1, j2
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq

    SELECT CASE ( j_range )
    CASE ( ne )
       j1 = 1
       j2 = n
    CASE ( lt )
       j1 = 1
       j2 = i-1
    CASE ( gt )
       j1 = i+1
       j2 = n
    CASE default
       WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error should never happen ', j_range
       STOP 'Impossible error in weight_1'
    END SELECT

    w = 0 ! Weight is zero unless we get through the following loop

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( ABS(j-i) <= 1 ) CYCLE ! Skip self and bonded neighbours

       rij    = ri(:) - r(:,j)     ! Separation vector
       rij_sq = SUM(rij**2)        ! Squared separation
       IF (  rij_sq < 1.0 ) RETURN ! Overlap detected, no need to continue

    END DO ! End loop over selected range of partners

    w = 1 ! No overlaps detected, so weight is one

  END FUNCTION weight_1

  FUNCTION qcount ( range ) RESULT ( q )
    INTEGER          :: q     ! Counts all square-well interactions
    REAL, INTENT(in) :: range ! Range of attractive well

    INTEGER :: i

    q = 0
    DO i = 1, n-1 ! Loop over atoms, checking upwards for partners
       q = q + qcount_1 ( r(:,i), i, gt, range )
    END DO ! End loop over atoms, checking upwards for partners

  END FUNCTION qcount

  FUNCTION qcount_1 ( ri, i, j_range, range ) RESULT ( q )
    INTEGER                           :: q       ! Counts square-well interactions for single atom
    REAL,    DIMENSION(3), INTENT(in) :: ri      ! Coordinates of atom of interest
    INTEGER,               INTENT(in) :: i       ! Index of atom of interest
    INTEGER,               INTENT(in) :: j_range ! Partner index range
    REAL,                  INTENT(in) :: range   ! Range of attractive well

    INTEGER            :: j, j1, j2
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq, range_sq

    range_sq = range**2

    SELECT CASE ( j_range )
    CASE ( ne )
       j1 = 1
       j2 = n
    CASE ( lt )
       j1 = 1
       j2 = i-1
    CASE ( gt )
       j1 = i+1
       j2 = n
    CASE default
       WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error should never happen ', j_range
       STOP 'Impossible error in qcount_1'
    END SELECT

    q = 0

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( ABS(j-i) <= 1 ) CYCLE ! Skip self and bonded neighbours

          rij    = ri(:) - r(:,j)            ! Separation vector
          rij_sq = SUM(rij**2)               ! Squared separation
          IF ( rij_sq < range_sq ) q = q + 1 ! Accumulate if within range

       END DO ! End loop over selected range of partners

  END FUNCTION qcount_1

  FUNCTION accept ( q_old, q_new, w_old, w_new )
    LOGICAL                       :: accept       ! Returns accept/reject decision based on
    INTEGER, INTENT(in)           :: q_old, q_new ! old & new energies, and if present,
    REAL,    INTENT(in), OPTIONAL :: w_old, w_new ! old & new weights

    ! This routine is essentially the Metropolis-Hastings formula, generalized
    ! to use a specified function of the energy in the exponent
    ! For NVT MC, s(q) is just E(q)/kT = -q/kT since energy is -q
    ! For Wang-Landau, s(q) is the entropy function, which is updated through the run
    ! In either case, weights may appear from the regrowth moves
    
    REAL :: zeta, delta

    delta = s(q_new) - s(q_old)

    CALL RANDOM_NUMBER(zeta)

    IF ( PRESENT ( w_new ) ) THEN ! include weights
       accept = (w_new/w_old)*EXP(-delta) > zeta
    ELSE ! do not include weights
       accept = EXP(-delta) > zeta
    END IF

  END FUNCTION accept

  SUBROUTINE zero_histogram
    h(:)   = 0.0
    hit(:) = .FALSE.
  END SUBROUTINE zero_histogram
  
  SUBROUTINE update_histogram ( q )
    INTEGER, INTENT(in) :: q ! Histogram bin to be updated

    ! This routine simply updates the probability histogram of (negative) energies
    ! In the case of Wang-Landau, it also updates the entropy histogram by ds
    ! and records a hit to show that this energy has been visited at least once
    ! The value of ds is stored in this module, and updated in the main program
    ! at the start of each block
    
    IF ( q <= nq .AND. q >= 0 ) THEN

       h(q) = h(q) + 1.0

       IF ( wl ) THEN ! only for Wang-Landau method
          s(q)   = s(q) + ds
          hit(q) = .TRUE.
       END IF

    ELSE
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'q out of range ', q, nq
       STOP 'Error in update_histogram'
    END IF

  END SUBROUTINE update_histogram

  FUNCTION histogram_flat ( flatness ) RESULT ( flat )
    LOGICAL            :: flat     ! Returns a signal that the histogram is "sufficiently flat"
    REAL,   INTENT(in) :: flatness ! Specified degree of flatness to pass the test 

    ! This routine is part of the Wang-Landau algorithm
    ! We only look at parts of the histogram that have been visited ("hit")
    ! The flatness should be in (0,1), higher values corresponding to flatter histograms

    REAL :: norm, avg

    IF ( flatness <= 0.0 .OR. flatness >= 1.0 ) THEN ! Ensure sensible value
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'Flatness error ', flatness
       STOP 'Error in histogram_flat'
    END IF

    norm = REAL(COUNT(hit))     ! How many entries
    avg  = SUM(h/norm,mask=hit) ! Average over just those entries

    IF ( avg < 0.0 ) THEN ! This should never happen, unless h somehow overflows
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'Error in h ', norm, avg
       STOP 'Error in histogram_flat'
    END IF

    flat = REAL(MINVAL(h,mask=hit)) > flatness*avg

  END FUNCTION histogram_flat

  SUBROUTINE write_histogram ( filename )
    USE, INTRINSIC :: iso_fortran_env, ONLY : iostat_end, iostat_eor

    CHARACTER(len=*), INTENT(in) :: filename

    INTEGER :: q, ioerr, his_unit

    OPEN ( newunit=his_unit, file=filename, status='replace', action='write', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)' ) 'Error opening ', filename, ioerr
       STOP 'Error in write_histogram'
    END IF

    DO q = 0, nq
       WRITE ( unit=his_unit, fmt='(i10,2es20.8)') q, h(q), s(q)
    END DO

    CLOSE ( unit=his_unit )

  END SUBROUTINE write_histogram

END MODULE mc_module

