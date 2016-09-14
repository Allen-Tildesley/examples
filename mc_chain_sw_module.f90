! mc_chain_sw_module.f90
! Monte Carlo, single chain, square wells
MODULE mc_module
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: model_description, allocate_arrays, deallocate_arrays
  PUBLIC :: regrow, cranks, pivots, weight, qcount, write_histogram, histogram_flat
  PUBLIC :: range, bond, n, nq, verbose, h, s, hit, r, ds, wl

  INTEGER                              :: n       ! number of atoms
  INTEGER                              :: nq      ! maximum anticipated energy
  LOGICAL                              :: verbose ! flag for verbose output
  LOGICAL                              :: wl      ! flag for Wang-Landau method
  REAL                                 :: bond    ! bond length
  REAL                                 :: range   ! range of attractive well
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r       ! atomic positions (3,n)
  REAL,    DIMENSION(:),   ALLOCATABLE :: h       ! histogram of q values (0:nq)
  REAL,    DIMENSION(:),   ALLOCATABLE :: s       ! histogram used in acceptance (0:nq)
  LOGICAL, DIMENSION(:),   ALLOCATABLE :: hit     ! identifies q-values that have been visited (0:nq)
  REAL                                 :: ds      ! entropy increment for Wang-Landau method

  REAL, DIMENSION(:,:), ALLOCATABLE :: r_old, r_new ! working arrays (3,n)

  REAL,    PARAMETER :: sigma = 1.0             ! core diameter
  INTEGER, PARAMETER :: lt = -1, gt = 1, ne = 0 ! options for j-range

CONTAINS

  SUBROUTINE model_description ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Hard-sphere chain with fixed bond length'
    WRITE ( unit=output_unit, fmt='(a)'           ) 'Square-well attractive potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ', sigma    
  END SUBROUTINE model_description

  SUBROUTINE allocate_arrays
    nq = 6*n ! anticipated maximum possible pair interactions
    ALLOCATE ( r(3,n), r_old(3,n), r_new(3,n) ) 
    ALLOCATE ( h(0:nq), s(0:nq), hit(0:nq) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, r_old, r_new, h, s, hit )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE regrow ( m_max, k_max, q, ratio ) ! regrow polymer
    USE maths_module, ONLY : random_integer, random_orientation_vector

    INTEGER, INTENT(in)    :: m_max  ! max atoms to regrow
    INTEGER, INTENT(in)    :: k_max  ! number of random tries per atom in regrow
    INTEGER, INTENT(inout) :: q      ! energy (to be updated)
    REAL,    INTENT(out)   :: ratio  ! acceptance ratio of moves

    INTEGER                      :: q_old, q_new ! old and new energies
    REAL                         :: w_old, w_new ! old and new weights
    INTEGER                      :: m            ! number of atoms to regrow
    INTEGER                      :: c            ! growth option
    INTEGER                      :: k            ! try index
    INTEGER, DIMENSION(k_max)    :: w            ! Rosenbluth weights (here integers, 0 or 1)
    REAL,    DIMENSION(3,k_max)  :: r_try        ! coordinates of trial atoms
    REAL,    DIMENSION(3)        :: u            ! vector
    INTEGER                      :: i            ! atom index
    REAL,              PARAMETER :: w_tol = 0.5  ! min weight tolerance (w takes integer values)

    ratio = 0.0
    r_old = r ! store old configuration
    q_old = q ! store old q

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
          r_try(:,k) = r(:,i-1) + bond*u ! generate trial position
          w(k)       = weight1 ( r_try(:,k), i, lt ) ! store overlap weight for this try
       END DO ! end loop over k_max tries

       IF ( SUM(w) == 0 ) THEN ! early exit
          r = r_old   ! restore original configuration
          q = q_old   ! redundant, but for clarity
          call update_histogram ( q )
          RETURN      ! and return
       END IF

       k      = pick(w)              ! pick winning position
       r(:,i) = r_try(:,k)           ! store winning position
       w_new  = w_new * REAL(SUM(w)) ! accumulate total weight

    END DO ! end loop to regrow last m atoms, computing new weight

    IF ( w_new < w_tol ) THEN ! this should be a redundant test
       r = r_old   ! restore original configuration
       q = q_old   ! redundant, but for clarity
       CALL update_histogram ( q )
       RETURN      ! no possible move
    END IF

    ! New configuration and weight are complete
    q_new = qcount() ! compute energy
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
       w(1) = 1            ! store overlap weight in try 1 (must be 1)

       DO k = 2, k_max ! loop over k_max-1 other tries
          CALL random_orientation_vector ( u )
          r_try(:,k) = r(:,i-1) + bond*u ! generate trial position
          w(k)   = weight1 ( r_try(:,k), i, lt ) ! store overlap weight
       END DO ! end loop over k_max-1 other tries

       r(:,i) = r_try(:,1)           ! restore winning position (always the original one)
       w_old  = w_old * REAL(SUM(w)) ! accumulate total weight

    END DO ! end loop to regrow last m atoms calculating old weight

    ! Choose either old or new configuration
    IF ( accept ( q_old, q_new, w_old, w_new ) ) THEN
       ratio = 1.0
       r     = r_new
       q     = q_new
    ELSE
       r = r_old
       q = q_old
    END IF

    CALL update_histogram ( q )

    IF ( verbose ) THEN
       IF ( q /= qcount() ) THEN
          WRITE ( unit = error_unit, fmt='(a,2i10)') 'Warning! ', q, qcount()
       END IF
    END IF

  END SUBROUTINE regrow

  SUBROUTINE pivots ( n_try, phi_max, q, ratio )
    USE maths_module, ONLY : random_integer, random_orientation_vector, rotate_vector

    INTEGER, INTENT(in)    :: n_try   ! number of pivots to try
    REAL,    INTENT(in)    :: phi_max ! maximum angle of pivot
    INTEGER, INTENT(inout) :: q       ! energy (to be updated)
    REAL,    INTENT(out)   :: ratio   ! move acceptance ratio

    INTEGER            :: q_old, q_new ! old and new energies
    INTEGER            :: j            ! pivot position
    INTEGER            :: k            ! which part to pivot
    INTEGER            :: i            ! atom index
    INTEGER            :: try, n_acc   ! trial and acceptance counters
    REAL, DIMENSION(3) :: u, rj, rij   ! local vectors
    REAL               :: phi          ! pivot angle

    n_acc = 0 ! zero acceptance counter

    tries: DO try = 1, n_try ! loop over pivot move attempts

       r_old = r ! store old configuration
       q_old = q ! store old q

       j = 1 + random_integer ( 1, n-2 ) ! pivot position (not at either end)
       k = random_integer ( 1, 2 )       ! which part to pivot (actually redundant here)

       CALL random_orientation_vector ( u ) ! rotation axis
       CALL RANDOM_NUMBER ( phi )
       phi = phi_max * ( 2.0*phi - 1.0 ) ! rotation angle

       SELECT CASE ( k )
       CASE ( 1 )
          r = r_old ! copy atoms
       CASE ( 2 )
          r = r_old(:,n:1:-1) ! copy and reverse atoms
       END SELECT

       ! Pivot, check for overlap in new configuration
       ! NB include overlap checks within rotated segment, because of roundoff
       rj = r(:,j)                            ! reference point
       DO i = j+1, n                          ! pivot these atoms
          rij = r(:,i) - rj                   ! relative vector of atom
          rij = rotate_vector ( phi, u, rij ) ! rotate relative vector
          r(:,i) = rj + rij                   ! new position

          IF ( weight1 ( r(:,i), i, lt ) == 0 ) THEN
             r = r_old  ! restore original configuration
             q = q_old  ! redundant but for clarity
             CALL update_histogram ( q )
             CYCLE tries ! and escape immediately
          END IF
       END DO

       q_new = qcount() ! calculate new energy
       r_new = r        ! store new configuration

       IF ( accept ( q_old, q_new ) ) THEN
          n_acc = n_acc + 1
          r = r_new
          q = q_new
       ELSE
          r = r_old
          q = q_old
       END IF
          call update_histogram ( q )

    END DO tries ! end loop over pivot move attempts

    ratio = REAL(n_acc) / REAL(n_try)

    IF ( verbose ) THEN
       IF ( q /= qcount() ) THEN
          WRITE ( unit = error_unit, fmt='(a,2i10)') 'Warning! ', q, qcount()
       END IF
    END IF

  END SUBROUTINE pivots

  SUBROUTINE cranks ( n_try, phi_max, q, ratio )
    USE maths_module, ONLY : random_integer, rotate_vector

    INTEGER, INTENT(in)    :: n_try   ! number of crankshaft moves to try
    REAL,    INTENT(in)    :: phi_max ! maximum crank angle
    INTEGER, INTENT(inout) :: q       ! energy (to be updated)
    REAL,    INTENT(out)   :: ratio   ! move acceptance ratio

    INTEGER              :: q_old, q_new ! old and new energies
    REAL, DIMENSION(3)   :: u, rj, rij   ! local vectors
    REAL                 :: phi          ! pivot angle
    INTEGER              :: i            ! atom index
    INTEGER              :: try, n_acc   ! trial and acceptance counters

    n_acc = 0 ! zero acceptance counter

    tries: DO try = 1, n_try ! Loop over crankshaft move attempts

       r_old = r ! store old configuration
       q_old = q ! store old q

       r = r_old

       i     = random_integer ( 1, n )            ! pick random atom to move
       q_new = q_old - qcount_1 ( r(:,i), i, ne ) ! subtract energies for moving atom

       CALL RANDOM_NUMBER ( phi )
       phi = ( 2.0*phi - 1.0 ) * phi_max ! rotation angle

       IF ( i == 1 ) THEN          ! rotate about 2--3 bond
          rj = r(:,2)              !   reference position
          u  = r(:,2) - r(:,3)     !   axis of rotation
       ELSE IF ( i == n ) THEN     ! rotate about (n-1)--(n-2) bond
          rj = r(:,n-1)            !   reference position
          u  = r(:,n-1) - r(:,n-2) !   axis of rotation
       ELSE                        ! rotate about (i-1)--(i+1) bond
          rj = r(:,i-1)            !   reference position
          u  = r(:,i+1) - r(:,i-1) !   axis of rotation
       END IF
       u = u / SQRT(SUM(u**2))     ! unit vector along rotation axis

       rij = r(:,i) - rj                   ! relative vector of atom
       rij = rotate_vector ( phi, u, rij ) ! rotate relative vector
       r(:,i) = rj + rij                   ! new position

       ! Check for overlaps in new configuration
       IF ( weight1 ( r(:,i), i, ne ) == 0 ) THEN
          r = r_old ! restore original configuration
          q = q_old ! redundant but for clarity
          call update_histogram ( q )
          CYCLE tries
       END IF
       q_new = q_new + qcount_1 ( r(:,i), i, ne ) ! add energies for moving atom
       r_new = r

       IF ( accept ( q_old, q_new ) ) THEN
          n_acc = n_acc + 1
          r     = r_new
          q     = q_new
       ELSE
          r = r_old
          q = q_old
       END IF
          call update_histogram ( q )

    END DO tries ! end loop over crankshaft move attempts

    ratio = REAL(n_acc) / REAL(n_try)

    IF ( verbose ) THEN
       IF ( q /= qcount() ) THEN
          WRITE ( unit = error_unit, fmt='(a,2i10)') 'Warning! ', q, qcount()
       END IF
    END IF

  END SUBROUTINE cranks

  FUNCTION pick ( w ) RESULT ( k ) !  select one of the tries according to overlap weight
    IMPLICIT NONE

    ! Argument
    INTEGER, DIMENSION(:), INTENT(in) :: w
    INTEGER                           :: k ! Result

    ! Local variables
    INTEGER :: cumw
    REAL    :: ran

    CALL RANDOM_NUMBER(ran)
    ran = ran*REAL(SUM(w))
    k=1
    cumw = w(1)
    DO WHILE ( cumw < ran )
       k = k+1
       cumw = cumw+w(k)
    ENDDO

  END FUNCTION pick

  FUNCTION weight () RESULT ( w ) ! compute configuration weight
    INTEGER :: w ! Result

    INTEGER :: i

    w = 1
    DO i = 1, n-1
       w = w  * weight1 ( r(:,i), i, gt )
    END DO

  END FUNCTION weight

  FUNCTION weight1 ( ri, i, j_range ) RESULT ( w ) !  compute single atom weight
    REAL,    DIMENSION(3), intent(in) :: ri      ! coordinates of atom of interest
    INTEGER,               INTENT(in) :: i       ! atom of interest
    INTEGER,               INTENT(in) :: j_range ! range of atoms to be studied
    INTEGER             :: w       ! Result

    INTEGER            :: j, j_lo, j_hi
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq

    SELECT CASE ( j_range )
    CASE ( ne )
       j_lo = 1
       j_hi = n
    CASE ( lt )
       j_lo = 1
       j_hi = i-1
    CASE ( gt )
       j_lo = i+1
       j_hi = n
    CASE default
       WRITE ( unit = error_unit, fmt='(a,i10)') 'This should never happen ', j_range
       STOP 'Impossible error in weight1'
    END SELECT

    DO j = j_lo, j_hi
       IF ( ABS(j-i) > 1 ) THEN ! skip self and bonded neighbours
          rij = ri(:) - r(:,j)
          rij_sq = SUM(rij**2)
          IF (  rij_sq < 1.0 ) THEN ! test against core sigma = 1.0
             IF ( verbose ) THEN
                WRITE ( unit = error_unit, fmt='(a,2i10,f15.5)') 'Warning core overlap = ', i, j, SQRT(rij_sq)
             END IF
             w = 0
             RETURN
          END IF ! end test against core sigma
       END IF ! end skip self and bonded neighbours
    END DO
    w = 1
  END FUNCTION weight1

  FUNCTION qcount () RESULT ( q ) ! counts total energy
    INTEGER :: q

    INTEGER :: i

    q = 0
    DO i = 1, n-1
       q = q + qcount_1 ( r(:,i), i, gt )
    END DO
  END FUNCTION qcount

  FUNCTION qcount_1 ( ri, i, j_range ) RESULT ( q ) !  counts single atom energy
    REAL,    DIMENSION(3), intent(in) :: ri      ! coordinates of atom of interest
    INTEGER,               INTENT(in) :: i       ! atom of interest
    INTEGER,               INTENT(in) :: j_range ! range of atoms to be studied
    INTEGER                           :: q

    INTEGER            :: j, j_lo, j_hi
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq, range_sq

    range_sq = range**2

    SELECT CASE ( j_range )
    CASE ( ne )
       j_lo = 1
       j_hi = n
    CASE ( lt )
       j_lo = 1
       j_hi = i-1
    CASE ( gt )
       j_lo = i+1
       j_hi = n
    CASE default
       WRITE ( unit = error_unit, fmt='(a,i10)') 'This should never happen ', j_range
       STOP 'Impossible error in qcount_1'
    END SELECT

    q = 0
    DO j = j_lo, j_hi
       IF ( ABS(i-j) > 1 ) THEN  ! skip self and bonded neighbours
          rij = ri(:) - r(:,j)
          rij_sq = SUM(rij**2)
          IF ( rij_sq < range_sq ) q = q + 1
       END IF ! end skip self and bonded neighbours
    END DO

  END FUNCTION qcount_1

  FUNCTION accept ( q_old, q_new, w_old, w_new ) ! Uses s(q) to determine acceptance
    INTEGER, INTENT(in)           :: q_old, q_new ! old & new energies
    REAL,    INTENT(in), OPTIONAL :: w_old, w_new ! old & new weights
    LOGICAL                       :: accept       ! Result

    REAL             :: zeta, delta
    REAL, PARAMETER  :: w_tol = 0.5 ! min weight tolerance (w takes integer values)

    IF ( PRESENT ( w_new ) ) THEN
       IF ( w_new < w_tol ) THEN
          accept = .FALSE. ! no possible move
          RETURN
       END IF
    END IF

    delta = s(q_old)-s(q_new)
    CALL RANDOM_NUMBER(zeta)
    IF ( PRESENT (w_new) ) THEN
       accept = zeta < (w_new/w_old)*EXP(delta)
    ELSE
       accept = zeta < EXP(delta)
    END IF

  END FUNCTION accept

  SUBROUTINE update_histogram ( q )
    INTEGER, INTENT(in) :: q
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
    REAL,   INTENT(in) :: flatness
    LOGICAL            :: flat

    REAL :: norm, avg

    norm = REAL(COUNT(hit))
    avg  = SUM(REAL(h)/norm,mask=hit)
    IF ( avg < 0.0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'Error in h ', norm, avg
       STOP 'Overflow in histogram_flat'
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

