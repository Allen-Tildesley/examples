! mc_chain_sw_module.f90
! Monte Carlo, single chain, square wells
MODULE mc_module

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             !
  ! to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     !
  ! published by Oxford University Press ("the publishers").                                       !
  !                                                                                                !
  ! LICENCE                                                                                        !
  ! Creative Commons CC0 Public Domain Dedication.                                                 !
  ! To the extent possible under law, the authors have dedicated all copyright and related         !
  ! and neighboring rights to this software to the PUBLIC domain worldwide.                        !
  ! This software is distributed without any warranty.                                             !
  ! You should have received a copy of the CC0 Public Domain Dedication along with this software.  !
  ! If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               !
  !                                                                                                !
  ! DISCLAIMER                                                                                     !
  ! The authors and publishers make no warranties about the software, and disclaim liability       !
  ! for all uses of the software, to the fullest extent permitted by applicable law.               !
  ! The authors and publishers do not recommend use of this software for any purpose.              !
  ! It is made freely available, solely to clarify points made in the text. When using or citing   !
  ! the software, you should not imply endorsement by the authors or publishers.                   !
  !------------------------------------------------------------------------------------------------!

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: regrow, crank, pivot, weight, qcount

  ! Public data
  INTEGER,                             PUBLIC :: n ! Number of atoms
  REAL,   DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Atomic positions (3,n)

  ! Private data
  REAL, DIMENSION(:,:), ALLOCATABLE :: r_old ! Working array (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: r_new ! Working array (3,n)

  REAL,    PARAMETER :: sigma = 1.0     ! Core diameter
  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range

CONTAINS

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Hard-sphere chain with fixed bond length'
    WRITE ( unit=output_unit, fmt='(a)'           ) 'Square-well attractive potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Diameter, sigma = ', sigma

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays
    IMPLICIT NONE

    ALLOCATE ( r(3,n), r_old(3,n), r_new(3,n) )

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, r_old, r_new )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE regrow ( s, m_max, k_max, bond, range, q, accepted )
    USE maths_module, ONLY : random_integer, random_vector, pick
    IMPLICIT NONE

    REAL,    DIMENSION(0:), INTENT(in)    :: s        ! Entropy table used in accept function
    INTEGER,                INTENT(in)    :: m_max    ! Max atoms to regrow
    INTEGER,                INTENT(in)    :: k_max    ! Number of random tries per atom in regrow
    REAL,                   INTENT(in)    :: bond     ! Bond length
    REAL,                   INTENT(in)    :: range    ! Range of attractive well
    INTEGER,                INTENT(inout) :: q        ! Energy, negative of (to be updated)
    LOGICAL,                INTENT(out)   :: accepted ! Indicates success or failure of move

    ! This routine carries out a single regrowth move, modifying the array r
    ! It is assumed that q contains the (negative of the) energy on input
    ! and this is updated to give the new value, if the move is accepted, on return

    ! A short sequence of m atoms (m<=m_max) is deleted and regrown in the CBMC manner
    ! We randomly select which end of the chain to apply each of these operations to
    ! At each stage, k_max different atom positions are tried
    ! The bond length is fixed throughout
    ! Weights used in the regrowth are athermal, computed only on the basis of the
    ! hard-core overlap part of the non-bonded interactions: essentially they count non-overlaps
    ! Hence they are suitable for use in both NVT and Wang-Landau simulations

    INTEGER                      :: q_old ! Old energy
    INTEGER                      :: q_new ! New energy
    REAL                         :: w_old ! Old weight
    REAL                         :: w_new ! New weight
    INTEGER                      :: m     ! Number of atoms to regrow
    INTEGER                      :: c     ! Growth option
    INTEGER                      :: k     ! Try index
    INTEGER, DIMENSION(k_max)    :: w     ! Rosenbluth weights (here integers, 0 or 1)
    REAL,    DIMENSION(3,k_max)  :: r_try ! Coordinates of trial atoms

    INTEGER            :: i
    REAL, DIMENSION(3) :: r1
    REAL, PARAMETER    :: w_tol = 0.5 ! Min weight tolerance (w takes integer values)

    r_old = r ! Store old configuration
    q_old = q ! Store old q

    m = random_integer ( 1, m_max ) ! Number of atoms to regrow
    c = random_integer ( 1, 4 )     ! Growth option

    ! PART 1: CONSTRUCT NEW CONFIGURATION WITH NEW WEIGHT

    SELECT CASE ( c )

    CASE ( 1 ) ! Remove from end and add to end
       r(:,1:n-m) = r_old(:,1:n-m) ! Copy first n-m atoms

    CASE ( 2 ) ! Remove from end and add to start
       r(:,1:n-m) = r_old(:,n-m:1:-1) ! Copy and reverse first n-m atoms

    CASE ( 3 ) ! Remove from start and add to start
       r(:,1:n-m) = r_old(:,n:m+1:-1) ! Copy and reverse last n-m atoms

    CASE ( 4 ) ! Remove from start and add to end
       r(:,1:n-m) = r_old(:,m+1:n) ! Copy last n-m atoms

    END SELECT

    ! Take the opportunity to place atom 1 at the origin
    r1         = r(:,1)
    r(:,1:n-m) = r(:,1:n-m) - SPREAD ( r1, dim=2, ncopies=n-m )

    w_new = 1.0
    DO i = n-m+1, n ! Loop to regrow last m atoms, computing new weight

       DO k = 1, k_max ! Loop over k_max tries
          r_try(:,k) = r(:,i-1) + bond * random_vector() ! Trial position in random direction from i-1
          w(k)       = weight_1 ( r_try(:,k), i, lt )    ! Store overlap weight for this try
       END DO ! End loop over k_max tries

       IF ( SUM(w) == 0 ) THEN ! Early exit if this happens at any stage
          r        = r_old   ! Restore original configuration
          q        = q_old   ! Redundant, but for clarity
          accepted = .FALSE. ! No possible move: reject
          RETURN
       END IF

       k      = pick ( w )           ! Pick winning try according to weights
       r(:,i) = r_try(:,k)           ! Store winning position
       w_new  = w_new * REAL(SUM(w)) ! Accumulate total weight

    END DO ! End loop to regrow last m atoms, computing new weight

    IF ( w_new < w_tol ) THEN ! This should be a redundant test
       r        = r_old   ! Restore original configuration
       q        = q_old   ! Redundant, but for clarity
       accepted = .FALSE. ! No possible move: reject
       RETURN
    END IF

    q_new = qcount ( range ) ! Compute full new nonbonded energy
    r_new = r                ! Store new configuration

    ! END OF PART 1: NEW CONFIGURATION AND WEIGHT ARE COMPLETE

    ! PART 2: RECONSTRUCT OLD CONFIGURATION WITH OLD WEIGHT

    SELECT CASE ( c )

    CASE ( 1, 2 ) ! Remove and hence reconstruct at end
       r = r_old(:,1:n) ! Copy all n atoms

    CASE ( 3, 4 ) ! Remove and reconstruct at start
       r = r_old(:,n:1:-1) ! Copy and reverse all n atoms

    END SELECT

    w_old = 1.0
    DO i = n-m+1, n ! Loop to regrow last m atoms calculating old weight

       ! Old position and weight are stored as try 1

       r_try(:,1) = r(:,i) ! Current position
       w(1)       = 1      ! Current weight must be 1

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
    IF ( accept ( s, q_old, q_new, w_old, w_new ) ) THEN
       r        = r_new
       q        = q_new
       accepted = .TRUE.
    ELSE
       r        = r_old
       q        = q_old
       accepted = .FALSE.
    END IF

  END SUBROUTINE regrow

  SUBROUTINE pivot ( s, phi_max, range, q, accepted )
    USE maths_module, ONLY : random_integer, random_vector, rotate_vector
    IMPLICIT NONE

    REAL,    DIMENSION(0:), INTENT(in)    :: s        ! Entropy table used in accept function
    REAL,                   INTENT(in)    :: phi_max  ! Maximum angle of pivot
    REAL,                   INTENT(in)    :: range    ! Range of attractive well
    INTEGER,                INTENT(inout) :: q        ! Energy, negative of (to be updated)
    LOGICAL,                INTENT(out)   :: accepted ! Indicates success or failure of move

    ! This routine carries out a pivot move, modifying the array r
    ! It is assumed that q contains the (negative of the) energy on input
    ! and this is updated to give the new value, if the move is accepted, on return

    ! An atom is picked at random, and the part of the chain lying to one side of it
    ! is rotated as a whole, by a random angle, about a randomly oriented axis
    ! There are no weights to take into account in the acceptance/rejection decision
    ! (the function weight_1 is used simply to indicate overlap / no overlap)

    INTEGER              :: q_old ! Old energy
    INTEGER              :: q_new ! New energy
    INTEGER              :: j     ! Index of pivot
    REAL,   DIMENSION(3) :: rj    ! Position of pivot

    REAL,   DIMENSION(3) :: r1, rij, axis
    INTEGER              :: i, k
    REAL                 :: zeta, phi

    IF ( n < 3 ) THEN
       accepted = .FALSE. ! Makes no sense to pivot such a short chain
       RETURN             ! Return immediately
    END IF

    r_old = r ! Store old configuration
    q_old = q ! Store old q

    k = random_integer ( 1, 2 ) ! Which part to pivot (actually redundant here)

    SELECT CASE ( k )
    CASE ( 1 )
       r = r_old ! Copy atoms
    CASE ( 2 )
       r = r_old(:,n:1:-1) ! Copy and reverse atoms
    END SELECT

    ! Take the opportunity to place atom 1 at the origin
    r1     = r(:,1)
    r(:,:) = r(:,:) - SPREAD ( r1, dim=2, ncopies=n )

    j  = random_integer ( 2, n-1 ) ! Pivot position (not at either end)
    rj = r(:,j) ! Pivot point

    DO i = 2, j ! Loop over static atoms, redundant but for roundoff
       IF ( weight_1 ( r(:,i), i, lt ) == 0 ) THEN ! Check for overlaps
          r        = r_old   ! Restore original configuration
          q        = q_old   ! Redundant but for clarity
          accepted = .FALSE. ! This move is rejected
          RETURN
       END IF ! End check for overlaps
    END DO ! End loop over static atoms, redundant but for roundoff

    ! Pivot, and check for overlap in new configuration
    ! NB include overlap checks within rotated segment, because of roundoff

    axis = random_vector ()            ! Pivot rotation axis
    CALL RANDOM_NUMBER ( zeta )        ! Uniform random number in range (0,1)
    phi = phi_max * ( 2.0*zeta - 1.0 ) ! Pivot rotation angle in desired range

    DO i = j+1, n ! Loop over moving atoms

       rij    = r(:,i) - rj                      ! Relative vector of atom
       rij    = rotate_vector ( phi, axis, rij ) ! Rotate relative vector
       r(:,i) = rj + rij                         ! New absolute position

       IF ( weight_1 ( r(:,i), i, lt ) == 0 ) THEN ! Check for overlaps
          r        = r_old   ! Restore original configuration
          q        = q_old   ! Redundant but for clarity
          accepted = .FALSE. ! This move is rejected
          RETURN
       END IF ! End check for overlaps

    END DO ! End loop over moving atoms

    q_new = qcount ( range ) ! Compute full new nonbonded energy
    r_new = r                ! Store new configuration

    IF ( accept ( s, q_old, q_new ) ) THEN
       r        = r_new
       q        = q_new
       accepted = .TRUE.
    ELSE
       r        = r_old
       q        = q_old
       accepted = .FALSE.
    END IF

  END SUBROUTINE pivot

  SUBROUTINE crank ( s, phi_max, range, q, accepted )
    USE maths_module, ONLY : random_integer, rotate_vector
    IMPLICIT NONE

    REAL,    DIMENSION(0:), INTENT(in)    :: s        ! Entropy table used in accept function
    REAL,                   INTENT(in)    :: phi_max  ! Maximum crank angle
    REAL,                   INTENT(in)    :: range    ! Range of attractive well
    INTEGER,                INTENT(inout) :: q        ! Energy, negative of (to be updated)
    LOGICAL,                INTENT(out)   :: accepted ! Indicates success or failure of move

    ! This routine carries out a crankshaft move, modifying the array r
    ! It is assumed that q contains the (negative of the) energy on input
    ! and this is updated to give the new value, if the move is accepted, on return

    ! An atom is picked at random. Unless it is an end-atom, a rotation axis
    ! is defined as the line joining the two atoms on either side, and the
    ! chosen atom is rotated about that axis by a random angle.
    ! In the case of end-atoms, the axis is defined by the line joining its
    ! nearest and next-nearest neighbours.
    ! There are no weights to take into account in the acceptance/rejection decision
    ! (the function weight_1 is used simply to indicate overlap / no overlap)

    INTEGER :: q_old ! Old energy
    INTEGER :: q_new ! New energy

    INTEGER            :: i
    REAL, DIMENSION(3) :: rj, rij, axis
    REAL               :: phi, zeta, norm

    IF ( n < 3 ) THEN
       accepted = .FALSE. ! Makes no sense to crank such a short chain
       RETURN             ! Return immediately
    END IF

    r_old = r ! Store old configuration
    q_old = q ! Store old q

    r = r_old ! Copy old position (somewhat redundant here)

    i     = random_integer ( 1, n )               ! Pick random atom to move
    q_new = q_old - qcount_1 ( r(:,i), i, range ) ! Subtract old energies for moving atom

    CALL RANDOM_NUMBER ( zeta )        ! Uniform random number in range (0,1)
    phi = ( 2.0*zeta - 1.0 ) * phi_max ! Rotation angle in desired range

    IF ( i == 1 ) THEN            ! Rotate about 2--3 bond
       axis = r(:,2) - r(:,3)     !   axis of rotation
       rj   = r(:,2)              !   reference position on axis
    ELSE IF ( i == n ) THEN       ! Rotate about (n-1)--(n-2) bond
       axis = r(:,n-1) - r(:,n-2) !   axis of rotation
       rj   = r(:,n-1)            !   reference position on axis
    ELSE                          ! Rotate about (i-1)--(i+1) bond
       axis = r(:,i+1) - r(:,i-1) !   axis of rotation
       rj   = r(:,i-1)            !   reference position on axis
    END IF

    norm = SQRT(SUM(axis**2)) ! Squared length of rotation axis
    axis = axis / norm        ! Unit vector along rotation axis

    rij    = r(:,i) - rj                      ! Relative vector of atom
    rij    = rotate_vector ( phi, axis, rij ) ! Rotate relative vector
    r(:,i) = rj + rij                         ! New absolute position

    IF ( weight_1 ( r(:,i), i ) == 0 ) THEN ! Check for overlaps
       r        = r_old   ! Restore original configuration
       q        = q_old   ! Redundant but for clarity
       accepted = .FALSE. ! This move is rejected
       RETURN
    END IF ! End check for overlaps

    q_new = q_new + qcount_1 ( r(:,i), i, range ) ! Add new energies for moving atom
    r_new = r

    IF ( accept ( s, q_old, q_new ) ) THEN
       r        = r_new
       q        = q_new
       accepted = .TRUE.
    ELSE
       r        = r_old
       q        = q_old
       accepted = .FALSE.
    END IF

  END SUBROUTINE crank

  FUNCTION weight () RESULT ( w )
    IMPLICIT NONE
    INTEGER :: w ! Returns configuration weight = 0 (overlap) or 1 (no overlap)

    ! Arithmetically w is the product of all the pair Boltzmann factors
    ! Here, each is 0 or 1, which allows us to treat them as Boolean variables
    ! but in the more general case they would be real and we would multiply them all together

    INTEGER :: i

    w = 0 ! Weight is zero unless we get through the following loop

    DO i = 1, n-1 ! Loop over atoms checking upwards for weight
       IF ( weight_1 ( r(:,i), i, gt ) == 0 ) RETURN ! Overlap detected, no need to continue
    END DO ! End loop over atoms checking upwards for weight

    w = 1 ! No overlaps detected, so weight is one

  END FUNCTION weight

  FUNCTION weight_1 ( ri, i, j_range ) RESULT ( w )
    IMPLICIT NONE
    INTEGER                           :: w       ! Returns weight (0 or 1) associated with single atom
    REAL,    DIMENSION(3), INTENT(in) :: ri      ! Coordinates of atom of interest
    INTEGER,               INTENT(in) :: i       ! Index of atom of interest
    INTEGER, OPTIONAL,     INTENT(in) :: j_range ! Partner index range

    ! Calculates weight associated with atom ri due to nonbonded interactions with other atoms
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! This is athermal, based on hard core nonbonded overlaps only
    ! The result is either 0 (overlap) or 1 (non overlap)
    ! Arithmetically w is the product of all the pair Boltzmann factors
    ! Here, each is 0 or 1, which allows us to treat them as Boolean variables
    ! but in the more general case they would be real and we would multiply them all together

    INTEGER            :: j, j1, j2
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq

    IF ( PRESENT ( j_range ) ) THEN
       SELECT CASE ( j_range )
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
    ELSE
       j1 = 1
       j2 = n
    END IF

    w = 0 ! Weight is zero unless we get through the following loop

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( ABS(j-i) <= 1 ) CYCLE ! Skip self and bonded neighbours

       rij    = ri(:) - r(:,j)    ! Separation vector
       rij_sq = SUM(rij**2)       ! Squared separation
       IF ( rij_sq < 1.0 ) RETURN ! Overlap detected, no need to continue

    END DO ! End loop over selected range of partners

    w = 1 ! No overlaps detected, so weight is one

  END FUNCTION weight_1

  FUNCTION qcount ( range ) RESULT ( q )
    IMPLICIT NONE
    INTEGER          :: q     ! Returns a count of all square-well interactions
    REAL, INTENT(in) :: range ! Range of attractive well

    ! Actual calculation is performed by function qcount_1

    INTEGER :: i

    q = 0

    DO i = 1, n-1 ! Loop over atoms, checking upwards for partners
       q = q + qcount_1 ( r(:,i), i, range, gt )
    END DO ! End loop over atoms, checking upwards for partners

  END FUNCTION qcount

  FUNCTION qcount_1 ( ri, i, range, j_range ) RESULT ( q )
    IMPLICIT NONE
    INTEGER                           :: q       ! Returns a count of square-well interactions
    REAL,    DIMENSION(3), INTENT(in) :: ri      ! Coordinates of atom of interest
    INTEGER,               INTENT(in) :: i       ! Index of atom of interest
    REAL,                  INTENT(in) :: range   ! Range of attractive well
    INTEGER, OPTIONAL,     INTENT(in) :: j_range ! Optional partner index range

    ! Counts square-well interactions between specified atom and a set of other atoms
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    INTEGER            :: j, j1, j2
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq, range_sq

    range_sq = range**2

    IF ( PRESENT ( j_range ) ) THEN
       SELECT CASE ( j_range )
       CASE ( lt )
          j1 = 1
          j2 = i-1
       CASE ( gt )
          j1 = i+1
          j2 = n
       CASE default ! should never happen
          WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error ', j_range
          STOP 'Impossible error in qcount_1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    q = 0

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( ABS(j-i) <= 1 ) CYCLE ! Skip self and bonded neighbours

       rij    = ri(:) - r(:,j)            ! Separation vector
       rij_sq = SUM(rij**2)               ! Squared separation
       IF ( rij_sq < range_sq ) q = q + 1 ! Accumulate if within range

    END DO ! End loop over selected range of partners

  END FUNCTION qcount_1

  FUNCTION accept ( s, q_old, q_new, w_old, w_new )
    USE maths_module, ONLY : metropolis
    IMPLICIT NONE
    LOGICAL                            :: accept       ! Returns accept/reject decision based on
    REAL,    DIMENSION(0:), INTENT(in) :: s            ! supplied entropy table,
    INTEGER,                INTENT(in) :: q_old, q_new ! old & new energies, and if present,
    REAL,    OPTIONAL,      INTENT(in) :: w_old, w_new ! old & new weights

    ! This routine is essentially the Metropolis-Hastings formula, generalized
    ! to use a specified tabulated function of the energy in the exponent
    ! For NVT MC, s(q) is just E(q)/kT = -q/kT since energy is -q
    ! For Wang-Landau, s(q) is the entropy function, which is updated through the run
    ! In either case, weights (real, but taking positive integer values here) may optionally appear

    INTEGER         :: nq
    REAL            :: delta
    REAL, PARAMETER :: w_tol = 0.5 ! Min weight tolerance (w takes integer values)

    ! Check that we are within bounds for the entropy table
    nq = UBOUND(s,1)
    IF ( q_new < 0 .OR. q_old < 0 .OR. q_new > nq .OR. q_old > nq ) THEN
       WRITE ( unit = error_unit, fmt='(a,3i10)') 'q out of bounds ', q_new, q_old, nq
       STOP 'nq error in accept'
    END IF

    delta = s(q_new) - s(q_old) ! Change in entropy function

    IF ( PRESENT ( w_new ) ) THEN ! Inclusion of weights

       IF ( w_old < w_tol ) THEN ! Should never happen
          WRITE ( unit = error_unit, fmt='(a,es15.8)') 'w_old error ', w_old
          STOP 'Impossible error in accept'
       END IF

       IF ( w_new < w_tol ) THEN ! Should really have been detected before
          accept = .FALSE. ! This move is rejected
          RETURN
       END IF

       ! Equivalent to exp(-delta) -> (w_new/w_old)*exp(-delta)
       delta = delta - LOG(w_new/w_old)

    END IF ! End inclusion of weights

    accept = metropolis ( delta )

  END FUNCTION accept

END MODULE mc_module
