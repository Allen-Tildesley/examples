! mc_chain_lj_module.f90
! Monte Carlo, single chain, LJ atoms
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
  PUBLIC :: regrow, potential, spring_pot

  ! Public data
  INTEGER,                             PUBLIC :: n ! Number of atoms
  REAL,   DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Atomic positions (3,n)

  ! Private data
  REAL, DIMENSION(:,:), ALLOCATABLE :: r_old ! Working array (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: r_new ! Working array (3,n)

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range

  ! Public derived type
  TYPE, PUBLIC :: potential_type ! A composite variable for interactions comprising
     REAL    :: pot ! the potential energy and
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
    c%pot = a%pot   +  b%pot
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  FUNCTION subtract_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the difference of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot   -  b%pot
    c%ovr = a%ovr .OR. b%ovr ! This is meaningless, but inconsequential
  END FUNCTION subtract_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    ! This model, specifically its collapse behaviour, is discussed in detail by
    ! F Calvo, JPK Doye, DJ Wales, J Chem Phys 116, 2642 (2002)
    ! A similar model, with rigid bond lengths, is discussed by
    ! A Irback, E Sandelin, J Chem Phys 110, 12256 (1999)
    ! Both types of model are also studied by
    ! P Grassberger, R Hegger, J Phys Cond Matt 7, 3089 (1995)
    
    WRITE ( unit=output_unit, fmt='(a)' ) 'LJ chain, no cutoff, no shift, no periodic box'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'   
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'   
    WRITE ( unit=output_unit, fmt='(a)' ) 'Harmonic spring bond potential'

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

  SUBROUTINE regrow ( temperature, m_max, k_max, bond, k_spring, accepted )
    USE maths_module, ONLY : random_integer, random_vector, pick
    IMPLICIT NONE

    REAL,    INTENT(in)    :: temperature ! Specified temperature
    INTEGER, INTENT(in)    :: m_max       ! Max atoms to regrow
    INTEGER, INTENT(in)    :: k_max       ! Number of random tries per atom in regrow
    REAL,    INTENT(in)    :: bond        ! Bond length
    REAL,    INTENT(in)    :: k_spring    ! Harmonic bond spring constant
    LOGICAL, INTENT(out)   :: accepted    ! Indicates acceptance or rejection of moves

    ! This routine carries out a single regrowth move, modifying the array r

    ! A short sequence of m atoms (m<=m_max) is deleted and regrown in the CBMC manner
    ! We randomly select which end of the chain to apply each of these operations to
    ! At each stage, k_max different atom positions are tried
    ! Function random_bond selects bond lengths according to the internal (harmonic) potential
    ! Rosenbluth weights are computed using the external (nonbonded) potential
    ! Acceptance/rejection is determined using these weights

    ! r_old and r_new are used as working arrays

    REAL                       :: w_old ! Old weight
    REAL                       :: w_new ! New weight
    INTEGER                    :: m     ! Number of atoms to regrow
    INTEGER                    :: c     ! Growth option
    INTEGER                    :: k     ! Try index
    REAL,   DIMENSION(k_max)   :: w     ! Rosenbluth weights (involving nonbonded interactions)
    REAL,   DIMENSION(3,k_max) :: r_try ! Coordinates of trial atoms

    INTEGER              :: i
    TYPE(potential_type) :: partial
    REAL                 :: d, d_max, std, zeta
    REAL, DIMENSION(3)   :: r1
    REAL, PARAMETER      :: w_tol = 1.0e-10 ! Min weight tolerance

    std   = SQRT(temperature/k_spring) ! Spring bond standard deviation
    d_max = 3.0*std                    ! Impose a limit on variation, say 3*std
    IF ( d_max > 0.75*bond ) THEN      ! Must not be too large, say 0.75*bond
       WRITE ( unit=error_unit, fmt='(a,2f15.6)' ) 'Spring bond strength error', d_max, bond
       STOP 'Error in regrow'
    END IF
    d_max = d_max + bond ! This is the actual max d allowed

    r_old = r ! Store old configuration

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

          d          = random_bond ( bond, std, d_max )  ! Generate random bond length around d=bond
          r_try(:,k) = r(:,i-1) + d * random_vector()    ! Trial position in random direction from i-1
          partial    = potential_1 ( r_try(:,k), i, lt ) ! Nonbonded interactions with earlier atoms

          IF ( partial%ovr ) THEN                      ! Overlap test
             w(k) = 0.0                                ! Store weight for this try (zero)
          ELSE                                         ! Safe to calculate weight
             w(k) = EXP ( -partial%pot / temperature ) ! Store weight for this try (Boltzmann factor)
          END IF                                       ! End overlap test

       END DO ! End loop over k_max tries

       IF ( SUM(w) < w_tol ) THEN ! Early exit if this happens at any stage
          r        = r_old   ! Restore original configuration
          accepted = .FALSE. ! No possible move: reject
          RETURN
       END IF

       k      = pick ( w )     ! Pick winning try according to weights
       r(:,i) = r_try(:,k)     ! Store winning position
       w_new  = w_new * SUM(w) ! Accumulate total weight

    END DO ! End loop to regrow last m atoms, computing new weight

    IF ( w_new < w_tol ) THEN ! Overall weight is too small
       r        = r_old   ! Restore original configuration
       accepted = .FALSE. ! No possible move: reject
       RETURN
    END IF

    r_new = r ! Store new configuration

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
       partial    = potential_1 ( r_try(:,1), i, lt ) ! Nonbonded energy with earlier atoms

       IF ( partial%ovr ) THEN                      ! Overlap test
          w(1) = 0.0                                ! Current weight is zero; should not happen
       ELSE                                         ! Safe to calculate weight
          w(1) = EXP ( -partial%pot / temperature ) ! Current weight given by Boltzmann factor
       END IF                                       ! End overlap test

       ! Remaining tries only required to compute weight

       DO k = 2, k_max ! Loop over k_max-1 other tries

          d          = random_bond ( bond, std, d_max )  ! Generate random bond length around d=bond
          r_try(:,k) = r(:,i-1) + d * random_vector ( )  ! Trial position in random direction from i-1
          partial    = potential_1 ( r_try(:,k), i, lt ) ! Nonbonded interactions with earlier atoms

          IF ( partial%ovr ) THEN                      ! Overlap test
             w(k) = 0.0                                ! Store weight for this try (zero)
          ELSE                                         ! Safe to calculate weight
             w(k) = EXP ( -partial%pot / temperature ) ! Store weight for this try (Boltzmann factor)
          END IF                                       ! End overlap test

       END DO ! End loop over k_max-1 other tries

       r(:,i) = r_try(:,1)     ! Restore winning position (always the original one)
       w_old  = w_old * SUM(w) ! Accumulate total weight

    END DO ! End loop to regrow last m atoms computing old weight

    IF ( w_old < w_tol ) THEN ! The old weight really should be non-zero
       WRITE ( unit=error_unit, fmt='(a,es20.8)' ) 'Old weight error', w_old
       STOP 'Impossible error in regrow'
    END IF

    ! END OF PART 2: OLD CONFIGURATION AND WEIGHT ARE COMPLETE

    ! Choose either old or new configuration according to weight
    ! All non-bonded Boltzmann factors are incorporated into the weights
    ! All spring-bond Boltzmann factors are included in the selection of d

    CALL RANDOM_NUMBER(zeta)
    IF ( zeta < ( w_new / w_old ) ) THEN
       r        = r_new
       accepted = .TRUE.
    ELSE
       r        = r_old
       accepted = .FALSE.
    END IF

  END SUBROUTINE regrow

  FUNCTION potential ( ) RESULT ( total )
    IMPLICIT NONE
    TYPE(potential_type) :: total ! Returns a composite of pot and ovr

    ! total%pot is the nonbonded potential energy for whole system
    ! total%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the value of total%pot should not be used
    ! Actual calculation is performed by function potential_1

    TYPE(potential_type) :: partial
    INTEGER              :: i

    IF ( n > SIZE(r,dim=2) ) THEN ! Should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Impossible error in potential'
    END IF

    total = potential_type ( pot=0.0, ovr=.FALSE. ) ! Initialize

    DO i = 1, n - 1
       partial = potential_1 ( r(:,i), i, gt )
       IF ( partial%ovr ) THEN
          total%ovr = .TRUE. ! Overlap detected
          RETURN             ! Return immediately
       END IF
       total = total + partial
    END DO

    total%ovr = .FALSE. ! No overlaps detected (redundant but for clarity)

  END FUNCTION potential

  FUNCTION potential_1 ( ri, i, j_range ) RESULT ( partial )
    IMPLICIT NONE
    TYPE(potential_type)            :: partial ! Returns a composite of pot and overlap for given atom
    REAL, DIMENSION(3), INTENT(in)  :: ri      ! Coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i       ! Index of atom of interest
    INTEGER, OPTIONAL,  INTENT(in)  :: j_range ! Optional partner index range

    ! partial%pot is the nonbonded potential energy of atom ri with a set of other atoms
    ! partial%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the value of partial%pot should not be used
    ! The coordinates in ri are not necessarily identical with those in r(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! Coordinates are assumed to be in LJ units where sigma = 1, no box, no periodic boundaries
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER              :: j, j1, j2
    REAL                 :: sr2, sr6, sr12, rij_sq
    REAL, DIMENSION(3)   :: rij
    REAL, PARAMETER      :: sr2_ovr = 1.77 ! Overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    IF ( n > SIZE(r,dim=2) ) THEN ! Should never happen
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

    partial = potential_type ( pot=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2 ! Loop over selected range of partners

       IF ( ABS(j-i) <= 1 ) CYCLE ! Skip self and bonded neighbours

       rij(:)   = ri(:) - r(:,j) ! Separation vector
       rij_sq   = SUM ( rij**2 ) ! Squared separation
       sr2      = 1.0 / rij_sq   ! (sigma/rij)**2
       pair%ovr = sr2 > sr2_ovr  ! Overlap if too close

       IF ( pair%ovr ) THEN
          partial%ovr = .TRUE. ! Overlap detected
          RETURN               ! Return immediately
       END IF

       sr6      = sr2**3
       sr12     = sr6**2
       pair%pot = sr12 - sr6 ! LJ pair potential (neither cut nor shifted)

       partial = partial + pair

    END DO ! End loop over selected range of partners
    
    partial%pot = partial%pot * 4.0 ! Numerical factor of 4*epsilon
    partial%ovr = .FALSE.           ! No overlaps detected (redundant but for clarity)

  END FUNCTION potential_1

  FUNCTION random_bond ( b, std, d_max ) RESULT ( d )
    USE maths_module, ONLY : random_normal
    IMPLICIT NONE
    REAL             :: d     ! Returns random bond length, generated around
    REAL, INTENT(in) :: b     ! specified mean bond length, Gaussian-distributed with
    REAL, INTENT(in) :: std   ! specified standard deviation, but subject to
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
    IMPLICIT NONE
    REAL             :: pot ! Returns internal spring potential energy
    REAL, INTENT(in) :: b   ! Specified bond length
    REAL, INTENT(in) :: k   ! Specified force constant

    INTEGER            :: i
    REAL               :: d
    REAL, DIMENSION(3) :: rij

    pot = 0.0

    DO i = 1, n-1 ! Loop over atoms
       rij = r(:,i) - r(:,i+1)       ! Bond vector
       d   = SQRT ( SUM ( rij**2 ) ) ! Bond distance
       pot = pot + (d-b)**2          ! Squared displacement
    END DO ! End loop over atoms

    ! Numerical factor
    pot = pot * 0.5 * k

  END FUNCTION spring_pot

END MODULE mc_module

