! maths_module.f90
! Routines for maths, random numbers, order parameters
MODULE maths_module

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

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit, iostat_end, iostat_eor

  IMPLICIT NONE
  PRIVATE

  ! Public random number routines
  PUBLIC :: init_random_seed, random_integer, random_normal, random_normals, pick
  PUBLIC :: random_vector
  PUBLIC :: random_vector_1, random_vector_2, random_vector_3
  PUBLIC :: random_perpendicular_vector
  PUBLIC :: random_rotate_vector, random_translate_vector
  PUBLIC :: random_rotate_vector_1, random_rotate_vector_2, random_rotate_vector_3, random_rotate_vector_4
  PUBLIC :: random_quaternion, random_rotate_quaternion
  PUBLIC :: metropolis

  ! Public low-level mathematical routines and string operations
  PUBLIC :: rotate_vector, rotate_quaternion, cross_product, outer_product, q_to_a, lowercase, solve

  ! Public order parameter calculations
  PUBLIC :: orientational_order, translational_order, nematic_order

  ! Generic interface for the pick functions
  INTERFACE pick
     MODULE PROCEDURE pick_i ! for integer weights
     MODULE PROCEDURE pick_r ! for real weights
  END INTERFACE pick

  ! Generic interface for the outer_product functions
  INTERFACE outer_product
     MODULE PROCEDURE outer_product_2 ! for 2 vectors giving a rank-2 result
     MODULE PROCEDURE outer_product_3 ! for 3 vectors giving a rank-3 result
     MODULE PROCEDURE outer_product_4 ! for 4 vectors giving a rank-4 result
     MODULE PROCEDURE outer_product_5 ! for 5 vectors giving a rank-5 result
  END INTERFACE outer_product

  ! Generic interface for the random_normals functions
  INTERFACE random_normals
     MODULE PROCEDURE random_normals_1 ! for rank 1 vector of normals
     MODULE PROCEDURE random_normals_2 ! for rank 2 vector of normals
  END INTERFACE random_normals

  ! Interface to select one of the random_vector algorithms
  INTERFACE random_vector
     MODULE PROCEDURE random_vector_1 ! choose an alternative one instead if you prefer
  END INTERFACE random_vector

  ! Interface to select one of the random_rotate_vector algorithms
  INTERFACE random_rotate_vector
     MODULE PROCEDURE random_rotate_vector_1 ! choose an alternative one instead if you prefer
  END INTERFACE random_rotate_vector

  ! Private data
  REAL, PARAMETER :: pi = 4.0*ATAN(1.0), twopi = 2.0*pi
  REAL, PARAMETER :: tol = 1.e-6

CONTAINS

  ! Routines associated with random number generation

  SUBROUTINE init_random_seed
    IMPLICIT NONE

    ! Initializes random number generator differently each time
    
    ! Prior to gfortran v7, calling the intrinsic RANDOM_SEED() routine initializes the
    ! random number generator with the same random seed to a default state, 
    ! which may result in the same sequence being generated every time.
    ! If you have gfortran v6, you may like to replace this routine init_random_seed 
    ! with the contents of file gnu_v6_init_random_seed.f90

    ! In gfortran v7 the random number generator was changed.
    ! Now, calling RANDOM_SEED() initializes the random number generator with random data
    ! retrieved from the operating system. Various other compilers behave the same way.
    ! We assume that this applies here.

    ! YOU SHOULD INVESTIGATE THE BEHAVIOUR FOR YOUR OWN COMPILER AND MACHINE IMPLEMENTATION 

    CALL RANDOM_SEED()
  END SUBROUTINE init_random_seed

  FUNCTION random_integer ( k1, k2 ) RESULT ( k )
    IMPLICIT NONE
    INTEGER             :: k      ! Returns a uniformly distributed random integer
    INTEGER, INTENT(in) :: k1, k2 ! in range [k1,k2] inclusive

    INTEGER :: k_min, k_max
    REAL    :: zeta

    CALL RANDOM_NUMBER ( zeta )
    k_min = MIN ( k1, k2 )
    k_max = MAX ( k1, k2 )
    k     = k_min + FLOOR ( (k_max-k_min+1)*zeta )

    ! Guard against small danger of roundoff
    IF ( k < k_min ) k = k_min 
    IF ( k > k_max ) k = k_max

  END FUNCTION random_integer

  FUNCTION random_normal ( mean, std ) RESULT ( r )
    IMPLICIT NONE
    REAL             :: r    ! Returns a normally-distributed random number with
    REAL, INTENT(in) :: mean ! specified mean and
    REAL, INTENT(in) :: std  ! specified standard deviation

    ! Box-Muller transform produces numbers in pairs
    ! We alternate between generating them, saving one for next time,
    ! and using the saved one from last time

    REAL, DIMENSION(2)      :: zeta
    REAL,              SAVE :: s
    LOGICAL,           SAVE :: saved = .FALSE.

    IF ( saved ) THEN     ! Saved number is available
       r = s              ! Normal, with mean=0, std=1
       r = mean + std * r ! Normal, with desired mean, std
       saved = .FALSE.    ! Flag to generate fresh numbers next time

    ELSE                                              ! Saved number is not available
       CALL RANDOM_NUMBER (zeta)                      ! Two uniformly distributed random numbers
       r = SQRT(-2.0*LOG(zeta(1)))*COS(twopi*zeta(2)) ! Normal, with mean=0, std=1
       s = SQRT(-2.0*LOG(zeta(1)))*SIN(twopi*zeta(2)) ! Also normal, with mean=0, std=1
       r = mean + std * r                             ! Normal, with desired mean, std
       saved = .TRUE.                                 ! Flag to use saved number next time

    END IF

  END FUNCTION random_normal

  SUBROUTINE random_normals_1 ( mean, std, r ) ! Normal random numbers
    IMPLICIT NONE
    REAL,               INTENT(in)  :: mean ! Specified mean and 
    REAL,               INTENT(in)  :: std  ! specified standard deviation, used to return
    REAL, DIMENSION(:), INTENT(out) :: r    ! vector of normal random numbers

    INTEGER :: i

    DO i = 1, SIZE(r)
       r(i) = random_normal ( mean, std )
    END DO

  END SUBROUTINE random_normals_1

  SUBROUTINE random_normals_2 ( mean, std, r ) ! Normal random numbers
    IMPLICIT NONE
    REAL,                 INTENT(in)  :: mean  ! Specified mean and 
    REAL,                 INTENT(in)  :: std   ! specified standard deviation, used to return
    REAL, DIMENSION(:,:), INTENT(out) :: r     ! array of normal random numbers

    INTEGER :: i, j

    DO j = 1, SIZE(r,dim=2)
       DO i = 1, SIZE(r,dim=1)
          r(i,j) = random_normal ( mean, std )
       END DO
    END DO

  END SUBROUTINE random_normals_2

  FUNCTION pick_r ( w ) RESULT ( k ) ! Pick amongst options with real weights
    IMPLICIT NONE
    INTEGER                        :: k ! Returns one of the options with probability proportional to
    REAL, DIMENSION(:), INTENT(in) :: w ! the supplied weights

    REAL :: cumw, zeta

    CALL RANDOM_NUMBER ( zeta ) ! Random number between 0 and 1
    zeta = zeta*SUM(w)          ! Scale up to total weight
    k    = 1
    cumw = w(1)
    DO ! Loop over possible outcomes
       IF ( zeta <= cumw ) EXIT ! Random number less than cumulative weight up to k
       k = k + 1
       IF ( k > SIZE(w) ) THEN ! Should never happen
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'k error ', k, SIZE(w)
          STOP 'Error in pick_r' 
       END IF
       cumw = cumw+w(k)
    END DO ! End loop over possible outcomes

  END FUNCTION pick_r

  FUNCTION pick_i ( w ) RESULT ( k ) ! Pick amongst options with integer weights
    IMPLICIT NONE
    INTEGER                           :: k ! Returns one of the options with probability proportional to
    INTEGER, DIMENSION(:), INTENT(in) :: w ! the supplied weights

    INTEGER :: cumw
    REAL    :: zeta

    CALL RANDOM_NUMBER ( zeta ) ! Random number between 0 and 1
    zeta = zeta*REAL(SUM(w))    ! Scale up to total weight
    k    = 1
    cumw = w(1)
    DO ! Loop over possible outcomes
       IF ( zeta <= REAL(cumw) ) EXIT ! Random number less than cumulative weight up to k
       k = k + 1
       IF ( k > SIZE(w) ) THEN ! Should never happen
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'k error ', k, SIZE(w)
          STOP 'Error in pick_i'
       END IF
       cumw = cumw+w(k)
    END DO ! End loop over possible outcomes

  END FUNCTION pick_i

  FUNCTION random_vector_1 () RESULT ( e ) ! 1st alternative algorithm
    IMPLICIT NONE
    REAL, DIMENSION(3) :: e ! Returns a uniformly sampled unit vector

    ! The vector is chosen uniformly within the cube surrounding the unit sphere
    ! Vectors lying outside the unit sphere are rejected
    ! Having found a vector within the unit sphere, it is normalized
    ! Essentially the same routine will work in 2d, or for quaternions in 4d

    REAL :: norm

    DO ! Loop until within unit sphere
       CALL RANDOM_NUMBER ( e ) ! 3 random numbers uniformly sampled in range (0,1)
       e    = 2.0 * e - 1.0     ! Now in range (-1,+1) i.e. within containing cube
       norm = SUM ( e**2 )      ! Square modulus
       IF ( norm < 1.0 ) EXIT   ! Within unit sphere
    END DO ! End loop until within unit sphere

    e = e / SQRT ( norm ) ! Normalize

  END FUNCTION random_vector_1

  FUNCTION random_vector_2 () RESULT ( e ) ! 2nd alternative algorithm
    IMPLICIT NONE
    REAL, DIMENSION(3) :: e ! Returns a uniformly sampled unit vector

    ! The polar angles are chosen from the correct distribution
    ! sampling phi and cos(theta) uniformly
    ! Then these are used to compute components of e

    REAL               :: c, s, phi
    REAL, DIMENSION(2) :: zeta

    CALL RANDOM_NUMBER ( zeta ) ! Two uniformly sampled random numbers in range (0,1)
    c   = 2.0 * zeta(1) - 1.0   ! Random cos(theta) uniformly sampled in range (-1,+1)

    IF ( c >= 1.0 ) THEN         ! Guard against very small chance of roundoff error
       s = 0.0                   ! Set sin(theta) to zero
    ELSE
       s   = SQRT ( 1.0 - c**2 ) ! Calculate sin(theta) from cos(theta), always positive
    END IF

    phi = zeta(2) * twopi               ! Random angle uniformly sampled in range (0,2*pi)
    e   = [ s*COS(phi), s*SIN(phi), c ] ! Random unit vector

  END FUNCTION random_vector_2

  FUNCTION random_vector_3 () RESULT ( e ) ! 3rd alternative algorithm
    IMPLICIT NONE
    REAL, DIMENSION(3) :: e ! Returns a uniformly sampled unit vector

    REAL, DIMENSION(2) :: zeta
    REAL               :: norm, f

    DO ! Loop until within unit disk
       CALL RANDOM_NUMBER ( zeta ) ! Two uniform random numbers between 0 and 1
       zeta = 2.0 * zeta - 1.0     ! Now each between -1 and 1
       norm = SUM ( zeta**2 )      ! Squared magnitude
       IF ( norm < 1.0 ) EXIT      ! Test whether inside unit disk
    END DO ! End loop until within unit disk

    f = 2.0 * SQRT ( 1.0 - norm )
    e = [ zeta(1) * f, zeta(2) * f, 1.0 - 2.0 * norm ] ! On surface of unit sphere

  END FUNCTION random_vector_3

  FUNCTION random_perpendicular_vector ( old ) RESULT ( e )
    IMPLICIT NONE
    REAL, DIMENSION(3)             :: e   ! Returns a uniformly sampled unit vector perpendicular to
    REAL, DIMENSION(3), INTENT(in) :: old ! the old vector 

    ! Note that we do not require the reference vector to be of unit length
    ! However we do require its length to be greater than a small tolerance!

    REAL, DIMENSION(3) :: n
    REAL               :: proj, norm

    norm = SUM ( old**2 ) ! Old squared length
    IF ( norm < tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es15.3)' ) 'old normalization error', norm, tol
       STOP 'Error in random_perpendicular_vector' 
    END IF
    n = old / SQRT(norm) ! Normalized old vector

    DO ! Loop until generated vector is not too small
       e    = random_vector ()     ! Randomly oriented unit vector
       proj = DOT_PRODUCT ( e, n ) ! Projection along old
       e    = e - proj * n         ! Make e perpendicular to old
       norm = SUM ( e**2 )         ! Squared length
       IF ( norm > tol ) EXIT      ! Accept, unless e is too small (which is unlikely)
    END DO ! End loop until generated vector is not too small

    e = e / SQRT ( norm ) ! Normalize

  END FUNCTION random_perpendicular_vector

  FUNCTION random_translate_vector ( dr_max, old ) RESULT ( r )
    IMPLICIT NONE
    REAL, DIMENSION(3)             :: r      ! Returns a vector translated by a random amount with
    REAL,               INTENT(in) :: dr_max ! maximum displacement relative to
    REAL, DIMENSION(3), INTENT(in) :: old    ! the old vector

    ! A randomly chosen vector is added to the old one

    REAL, DIMENSION(3) :: zeta ! Random numbers

    CALL RANDOM_NUMBER ( zeta )   ! Three uniform random numbers in range (0,1)
    zeta = 2.0*zeta - 1.0         ! Now in range (-1,+1)
    r(:) = old(:) + zeta * dr_max ! Move to new position

  END FUNCTION random_translate_vector

  FUNCTION random_rotate_vector_1 ( angle_max, old ) RESULT ( e ) ! 1st alternative algorithm
    IMPLICIT NONE
    REAL, DIMENSION(3)             :: e         ! Returns a unit vector rotated by a
    REAL,               INTENT(in) :: angle_max ! maximum angle (in radians) relative to
    REAL, DIMENSION(3), INTENT(in) :: old       ! the old vector

    ! A small randomly chosen vector is added to the old one, and the result renormalized
    ! Provided angle_max is << 1, it is approximately the maximum rotation angle (in radians)
    ! The magnitude of the rotation is not uniformly sampled, but this should not matter

    ! Note that the old vector should be normalized and we test for this

    REAL :: norm

    norm = SUM ( old**2 ) ! Old squared length
    IF ( ABS ( norm - 1.0 ) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'old normalization error', norm, tol
       STOP 'Error in random_rotate_vector_1'
    END IF

    ! Choose new orientation by adding random small vector
    e    = old + angle_max * random_vector ()
    norm = SUM ( e**2 )
    e    = e / SQRT(norm) ! Normalize

  END FUNCTION random_rotate_vector_1

  FUNCTION random_rotate_vector_2 ( angle_max, old ) RESULT ( e ) ! 2nd alternative algorithm
    IMPLICIT NONE
    REAL, DIMENSION(3)             :: e         ! Returns a unit vector rotated by a
    REAL,               INTENT(in) :: angle_max ! maximum angle (in radians) relative to
    REAL, DIMENSION(3), INTENT(in) :: old       ! the old vector

    ! The magnitude of the rotation is uniformly sampled
    ! The rotation axis is chosen randomly, perpendicular to the old orientation

    ! Note that the old vector should be normalized and we test for this

    REAL, DIMENSION(3) :: perp
    REAL               :: norm, angle, zeta

    norm = SUM ( old**2 ) ! Old squared length
    IF ( ABS ( norm - 1.0 ) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'old normalization error', norm, tol
       STOP 'Error in random_rotate_vector_2'
    END IF

    perp = random_perpendicular_vector ( old ) ! Choose unit vector perpendicular to old

    CALL RANDOM_NUMBER ( zeta )              ! Uniform random number between 0 and 1
    angle = ( 2.0 * zeta - 1.0 ) * angle_max ! Uniform random angle in desired range

    e    = old * COS ( angle ) + perp * SIN ( angle ) ! Rotation in the (old,perp) plane
    norm = SUM ( e**2 )
    e    = e / SQRT(norm) ! Normalize (should be redundant)

  END FUNCTION random_rotate_vector_2

  FUNCTION random_rotate_vector_3 ( angle_max, old ) RESULT ( e ) ! 3rd alternative algorithm
    IMPLICIT NONE
    REAL, DIMENSION(3)             :: e         ! Returns a unit vector rotated by a
    REAL,               INTENT(in) :: angle_max ! maximum angle (in radians) relative to
    REAL, DIMENSION(3), INTENT(in) :: old       ! the old vector

    ! The magnitude of the rotation is uniformly sampled
    ! The rotation axis is a Cartesian axis selected at random
    ! Ref: Barker and Watts, Chem Phys Lett 3, 144 (1969)

    ! Note that the old vector should be normalized and we test for this

    INTEGER            :: k
    REAL, DIMENSION(3) :: axis
    REAL               :: angle, norm, zeta

    norm = SUM ( old**2 ) ! Old squared length
    IF ( ABS ( norm - 1.0 ) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'old normalization error', norm, tol
       STOP 'Error in random_rotate_vector_3'
    END IF

    axis    = 0.0
    k       = random_integer (1,3) ! Random axis choice x=1, y=2, z=3
    axis(k) = 1.0

    CALL RANDOM_NUMBER ( zeta )              ! Uniform random number between 0 and 1
    angle = ( 2.0 * zeta - 1.0 ) * angle_max ! Uniform random angle in desired range

    e = rotate_vector ( angle, axis, old ) ! General rotation function

  END FUNCTION random_rotate_vector_3

  FUNCTION random_rotate_vector_4 ( angle_max, old ) RESULT ( e ) ! 4th alternative algorithm
    IMPLICIT NONE
    REAL, DIMENSION(3)             :: e         ! Returns a unit vector rotated by a
    REAL, INTENT(in)               :: angle_max ! maximum angle (in radians) relative to
    REAL, DIMENSION(3), INTENT(in) :: old       ! the old vector

    ! Ref: Marsaglia, Ann Maths Stat 43, 645 (1972)
    ! Uses a rejection technique to create a trial orientation
    ! subject to the constraint that the cosine of the angle
    ! turned through is greater than cos(angle_max)
    ! Not very efficient if angle_max is small

    ! Note that the old vector should be normalized and we test for this

    REAL :: cos_min, norm

    norm = SUM ( old**2 ) ! Old squared length
    IF ( ABS ( norm - 1.0 ) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'old normalization error', norm, tol
       STOP 'Error in random_rotate_vector_4'
    END IF

    cos_min = COS ( angle_max )

    DO ! Loop until close enough
       e = random_vector ( )                        ! Choose completely random unit vector
       IF ( DOT_PRODUCT ( e, old ) > cos_min ) EXIT ! Close enough
    END DO ! End loop until close enough

  END FUNCTION random_rotate_vector_4

  FUNCTION random_quaternion () RESULT ( e )
    IMPLICIT NONE
    REAL, DIMENSION(0:3) :: e ! Returns a uniformly sampled unit quaternion

    REAL, DIMENSION(2) :: zeta
    REAL               :: norm1, norm2, f

    DO ! Loop until within unit disk
       CALL RANDOM_NUMBER ( zeta ) ! Two uniform random numbers between 0 and 1
       zeta = 2.0 * zeta - 1.0     ! Now each between -1 and 1
       norm1 = SUM ( zeta**2 )     ! Squared magnitude
       IF ( norm1 < 1.0 ) EXIT     ! Test for within unit disk
    END DO ! End loop until within unit disk

    e(0) = zeta(1)
    e(1) = zeta(2)

    DO ! Loop until within unit disk
       CALL RANDOM_NUMBER ( zeta ) ! Two uniform random numbers between 0 and 1
       zeta = 2.0 * zeta - 1.0     ! Now each between -1 and 1
       norm2 = SUM ( zeta**2 )     ! Squared magnitude
       IF ( norm2 < 1.0 ) EXIT     ! Test for within unit disk
    END DO ! End loop until within unit disk

    f = SQRT ( ( 1.0 - norm1 ) / norm2 )
    e(2) = zeta(1)*f
    e(3) = zeta(2)*f

  END FUNCTION random_quaternion

  FUNCTION random_rotate_quaternion ( angle_max, old ) RESULT ( e )
    IMPLICIT NONE
    REAL, DIMENSION(0:3)             :: e         ! Returns a unit quaternion rotated by a
    REAL,                 INTENT(in) :: angle_max ! maximum angle (in radians) relative to
    REAL, DIMENSION(0:3), INTENT(in) :: old       ! the old quaternion

    ! Note that the reference quaternion should be normalized and we test for this

    REAL, DIMENSION(3) :: axis
    REAL               :: zeta, angle, norm

    norm = SUM ( old**2 ) ! Old squared length
    IF ( ABS ( norm - 1.0 ) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'old normalization error', norm, tol
       STOP 'Error in random_rotate_quaternion'
    END IF

    axis = random_vector ( )               ! Choose random unit vector
    CALL RANDOM_NUMBER ( zeta )            ! Random number between 0 and 1
    angle = ( 2.0*zeta - 1.0 ) * angle_max ! Uniform random angle in desired range

    e = rotate_quaternion ( angle, axis, old ) ! General rotation function

  END FUNCTION random_rotate_quaternion

  FUNCTION metropolis ( delta ) RESULT ( accept ) ! Conduct Metropolis test, with safeguards
    IMPLICIT NONE
    LOGICAL          :: accept ! Returns decision
    REAL, INTENT(in) :: delta  ! Negative of argument of exponential

    REAL            :: zeta
    REAL, PARAMETER :: exponent_guard = 75.0

    IF ( delta > exponent_guard ) THEN ! Too high, reject without evaluating
       accept = .FALSE.
    ELSE IF ( delta < 0.0 ) THEN ! Downhill, accept without evaluating
       accept = .TRUE.
    ELSE
       CALL RANDOM_NUMBER ( zeta ) ! Uniform random number in range (0,1)
       accept = EXP(-delta) > zeta ! Metropolis test
    END IF

  END FUNCTION metropolis

  ! Low level mathematical operations and string manipulation

  FUNCTION rotate_vector ( angle, axis, old ) RESULT ( e )
    IMPLICIT NONE
    REAL, DIMENSION(3)             :: e     ! Returns a vector rotated by
    REAL,               INTENT(in) :: angle ! specified rotation angle (in radians) about
    REAL, DIMENSION(3), INTENT(in) :: axis  ! specified rotation axis relative to
    REAL, DIMENSION(3), INTENT(in) :: old   ! old vector

    ! Note that the axis vector should be normalized and we test for this
    ! In general, the old vector need not be normalized, and the same goes for the result
    ! although quite often in our applications they will be

    REAL :: proj, c, s, norm

    norm = SUM ( axis**2 ) ! Axis squared length
    IF ( ABS ( norm - 1.0 ) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'axis normalization error', norm, tol
       STOP 'Error in rotate_vector'
    END IF

    c    = COS ( angle )
    s    = SIN ( angle )
    proj = DOT_PRODUCT ( axis, old ) ! The two vectors need not be perpendicular

    ! Standard (Goldstein) rotation formula
    e = c * old + ( 1.0 - c ) * proj * axis + s * cross_product ( axis, old )

  END FUNCTION rotate_vector

  FUNCTION rotate_quaternion ( angle, axis, old ) RESULT ( e )
    IMPLICIT NONE
    REAL, DIMENSION(0:3)             :: e     ! Returns a quaternion rotated by
    REAL,                 INTENT(in) :: angle ! specified rotation angle (in radians) about
    REAL, DIMENSION(3),   INTENT(in) :: axis  ! specified rotation axis relative to
    REAL, DIMENSION(0:3), INTENT(in) :: old   ! old quaternion

    ! Note that the axis vector should be normalized and we test for this
    ! In general, the old quaternion need not be normalized, and the same goes for the result
    ! although in our applications we only ever use unit quaternions (to represent orientations)

    REAL                 :: norm
    REAL, DIMENSION(0:3) :: rot

    norm = SUM ( axis**2 ) ! Axis squared length
    IF ( ABS ( norm - 1.0 ) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'axis normalization error', norm, tol
       STOP 'Error in rotate_quaternion'
    END IF

    ! Standard formula for rotation quaternion, using half angles
    rot(0)   = COS(0.5*angle)
    rot(1:3) = SIN(0.5*angle)*axis

    e = quatmul ( rot, old ) ! Apply rotation to old quaternion

  END FUNCTION rotate_quaternion

  FUNCTION quatmul ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(0:3)             :: c    ! Returns quaternion product of
    REAL, DIMENSION(0:3), INTENT(in) :: a, b ! two supplied quaternions

    c(0) = a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
    c(1) = a(1)*b(0) + a(0)*b(1) - a(3)*b(2) + a(2)*b(3)
    c(2) = a(2)*b(0) + a(3)*b(1) + a(0)*b(2) - a(1)*b(3)
    c(3) = a(3)*b(0) - a(2)*b(1) + a(1)*b(2) + a(0)*b(3)

  END FUNCTION quatmul

  FUNCTION cross_product ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3)             :: c    ! Returns vector cross product of
    REAL, DIMENSION(3), INTENT(in) :: a, b ! two supplied vectors
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  END FUNCTION cross_product

  FUNCTION outer_product_2 ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(:),              INTENT(IN) :: a, b ! Given two supplied vectors,
    REAL, DIMENSION(SIZE(a),SIZE(b))            :: c    ! returns their rank-2 outer product

    INTEGER :: i, j

    DO j = 1, SIZE(b)
       DO i = 1, SIZE(a)
          c(i,j) = a(i) * b(j)
       END DO
    END DO

    ! The following one-line statement is equivalent, but the above loops are clearer
    ! c = SPREAD(a,dim=2,ncopies=SIZE(b)) * SPREAD(b,dim=1,ncopies=SIZE(a))

  END FUNCTION outer_product_2

  FUNCTION outer_product_3 ( a, b, c ) RESULT (d)
    IMPLICIT NONE
    REAL, DIMENSION(:),                      INTENT(IN) :: a, b, c ! Given three supplied vectors,
    REAL, DIMENSION(SIZE(a),SIZE(b),SIZE(c))            :: d       ! returns their rank-3 outer product

    INTEGER :: i, j, k

    DO k = 1, SIZE(c)
       DO j = 1, SIZE(b)
          DO i = 1, SIZE(a)
             d(i,j,k) = a(i) * b(j) * c(k)
          END DO
       END DO
    END DO

  END FUNCTION outer_product_3

  FUNCTION outer_product_4 ( a, b, c, d ) RESULT (e)
    IMPLICIT NONE
    REAL, DIMENSION(:),                   INTENT(IN) :: a, b, c, d ! Given four supplied vectors,
    REAL, DIMENSION(SIZE(a),SIZE(b),SIZE(c),SIZE(d)) :: e          ! returns their rank-4 outer product

    INTEGER :: i, j, k, l

    DO l = 1, SIZE(d)
       DO k = 1, SIZE(c)
          DO j = 1, SIZE(b)
             DO i = 1, SIZE(a)
                e(i,j,k,l) = a(i) * b(j) * c(k) * d(l)
             END DO
          END DO
       END DO
    END DO

  END FUNCTION outer_product_4

  FUNCTION outer_product_5 ( a, b, c, d, e ) RESULT (f)
    IMPLICIT NONE
    REAL, DIMENSION(:),                           INTENT(IN) :: a, b, c, d, e ! Given five  supplied vectors,
    REAL, DIMENSION(SIZE(a),SIZE(b),SIZE(c),SIZE(d),SIZE(e)) :: f             ! returns their rank-5 outer product

    INTEGER :: i, j, k, l, m

    DO m = 1, SIZE(e) 
       DO l = 1, SIZE(d)
          DO k = 1, SIZE(c)
             DO j = 1, SIZE(b)
                DO i = 1, SIZE(a)                
                   f(i,j,k,l,m) = a(i) * b(j) * c(k) * d(l) * e(m)
                END DO
             END DO
          END DO
       END DO
    END DO

  END FUNCTION outer_product_5

  FUNCTION lowercase ( oldstring ) RESULT ( newstring )
    IMPLICIT NONE
    CHARACTER(len=*),             INTENT(in) :: oldstring ! Given a supplied string,
    CHARACTER(len=LEN(oldstring))            :: newstring ! returns a copy converted to lowercase

    INTEGER :: i, k 

    ! Leaves non-alphabetic characters unchanged

    DO i = 1, LEN(oldstring) 
       k = IACHAR(oldstring(i:i)) 
       IF ( k >= IACHAR('A') .AND. k <= IACHAR('Z') ) THEN 
          k = k + IACHAR('a') - IACHAR('A') 
          newstring(i:i) = ACHAR(k)
       ELSE
          newstring(i:i) = oldstring(i:i)
       END IF
    END DO
  END FUNCTION lowercase

  ! Order parameter routines

  FUNCTION translational_order ( r, k ) RESULT ( order )
    IMPLICIT NONE
    REAL                                          :: order ! Returns a translational order parameter from
    REAL,    DIMENSION(:,:), INTENT(in)           :: r     ! set of molecular position vectors (3,n), and
    INTEGER, DIMENSION(3),   INTENT(in), OPTIONAL :: k     ! lattice reciprocal vector (integer)

    ! Calculate the "melting factor" for translational order 
    ! based on a single k-vector characterizing the original lattice
    ! and commensurate with the periodic box
    ! It is assumed that both r and k are in box=1 units
    ! k = (l,m,n) where l,m,n are integers, must be multiplied by 2*pi to get real k
    ! If optional argument k is omitted, we default to a choice
    ! based on the FCC lattice, if this makes sense based on the number of atoms
    ! order = 1 when all atoms are on their lattice positions
    ! order = 1/sqrt(n), approximately, for disordered positions

    INTEGER            :: i, n, nc
    REAL, DIMENSION(3) :: k_real
    REAL               :: kr
    COMPLEX            :: rho ! Fourier component of single-particle density

    IF ( SIZE(r,dim=1) /= 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error in r dimension ', SIZE(r,dim=1)
       STOP 'Error in translational_order'
    END IF
    n = SIZE(r,dim=2)

    IF ( PRESENT ( k ) ) THEN
       k_real = twopi * REAL ( k )
    ELSE                                          ! Make arbitrary choice assuming FCC
       nc = NINT ( ( REAL(n)/4.0 ) ** (1.0/3.0) ) ! number of FCC unit cells
       IF ( 4*nc**3 /= n ) THEN ! Check that n does indeed correspond to FCC lattice
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Error in value of n ', 4*nc**3, n
          STOP 'Error in translational_order'
       END IF
       k_real = twopi * REAL( [-nc,nc,-nc] )      ! arbitrary fcc reciprocal vector
    END IF

    rho = ( 0.0, 0.0 )

    DO i = 1, n
       kr  = DOT_PRODUCT ( k_real, r(:,i) )
       rho = rho + CMPLX ( COS(kr), SIN(kr) )
    END DO

    rho   = rho / REAL(n)
    order = REAL ( CONJG(rho)*rho )

  END FUNCTION translational_order

  FUNCTION orientational_order ( e ) RESULT ( order )
    IMPLICIT NONE
    REAL                             :: order ! Returns a crystal orientational order parameter from
    REAL, DIMENSION(:,:), INTENT(in) :: e     ! a set of molecular orientation vectors (3,n)

    ! Calculates an orientational order parameter to monitor "melting"
    ! The parameter depends completely on knowing the orientations of the molecules
    ! in the original crystal lattice, and here we assume a specific alpha-fcc crystal
    ! of the same kind as was set up in initialize.f90 and initialize_module.f90
    ! Four molecules per unit cell, each pointing along a body-diagonal
    ! Order parameter can be a low-ranking (e.g. 1st or 2nd) Legendre polynomial

    ! Note that this is not the same as the order parameter characterizing a nematic liquid crystal

    INTEGER :: n, nc, i, i0
    REAL    :: c

    REAL, DIMENSION(3,4), PARAMETER :: e0 = RESHAPE (  SQRT(1.0/3.0)*[ &
         &  1.0,  1.0,  1.0,    1.0, -1.0, -1.0,  &
         & -1.0,  1.0, -1.0,   -1.0, -1.0,  1.0 ],[3,4] ) ! Orientations in unit cell

    IF ( SIZE(e,dim=1) /= 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error in e dimension ', SIZE(e,dim=1)
       STOP 'Error in orientational_order'
    END IF
    n  = SIZE(e,dim=2)
    nc = NINT ( ( REAL(n)/4.0 ) ** (1.0/3.0) )
    IF ( 4*nc**3 /= n ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Error in value of n ', 4*nc**3, n
       STOP 'Error in orientational_order'
    END IF

    order = 0.0

    DO i = 1, n
       i0    = MODULO ( i, 4 ) + 1              ! Select appropriate original orientation
       c     = DOT_PRODUCT ( e(:,i), e0(:,i0) ) ! Cosine of angle
       order = order + 1.5*c**2 - 0.5           ! Second Legendre polynomial
    END DO

    order = order / REAL ( n )

  END FUNCTION orientational_order

  FUNCTION nematic_order ( e ) RESULT ( order )
    IMPLICIT NONE
    REAL                             :: order ! Returns a nematic orientational order parameter from
    REAL, DIMENSION(:,:), INTENT(in) :: e     ! a set of molecular orientation vectors (3,n)

    ! Calculate the nematic order parameter <P2(cos(theta))>
    ! where theta is the angle between a molecular axis and the director
    ! which is the direction that maximises the order parameter
    ! This is obtained by finding the largest eigenvalue of
    ! the 3x3 second-rank traceless order tensor

    ! Note that this is not the same as the order parameter characterizing a crystal

    INTEGER              :: i, n
    REAL, DIMENSION(3,3) :: q         ! order tensor
    REAL                 :: h, g, psi ! used in eigenvalue calculation

    IF ( SIZE(e,dim=1) /= 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error in e dimension ', SIZE(e,dim=1)
       STOP 'Error in nematic_order'
    END IF
    n = SIZE(e,dim=2)

    ! Order tensor: outer product of each orientation vector, summed over n molecules
    q = SUM ( SPREAD ( e, dim=2, ncopies=3) * SPREAD ( e, dim=1, ncopies=3 ), dim = 3 )
    q = 1.5 * q / REAL(n)                ! Normalize
    FORALL (i=1:3) q(i,i) = q(i,i) - 0.5 ! Make traceless

    ! Trigonometric solution of characteristic cubic equation, assuming real roots

    h =      q(1,1) * q(2,2) - q(1,2) * q(2,1) &
         & + q(2,2) * q(3,3) - q(2,3) * q(3,2) &
         & + q(3,3) * q(1,1) - q(3,1) * q(1,3)
    h = h / 3.0

    g =      q(1,1) * q(2,2) * q(3,3) - q(1,1) * q(2,3) * q(3,2) &
         & + q(1,2) * q(2,3) * q(3,1) - q(2,2) * q(3,1) * q(1,3) &
         & + q(2,1) * q(3,2) * q(1,3) - q(3,3) * q(1,2) * q(2,1)

    h = SQRT(-h)
    psi = -0.5 * g / h**3
    IF ( psi < -1.0 ) psi = -1.0
    IF ( psi >  1.0 ) psi =  1.0
    psi = ACOS(psi)
    h = -2.0*h

    ! Select largest root
    order = MAXVAL ( [ h*COS(psi/3.0), h*COS((psi+2.0*pi)/3.0), h*COS((psi+4.0*pi)/3.0) ] ) 

  END FUNCTION nematic_order

  FUNCTION solve ( a, b ) RESULT ( x )
    IMPLICIT NONE
    REAL, DIMENSION(3)               :: x ! Returns a vector x, the solution of ax=b, for
    REAL, DIMENSION(3,3), INTENT(in) :: a ! supplied 3x3 matrix, and
    REAL, DIMENSION(3),   INTENT(in) :: b ! supplied constant vector

    ! Use Cramer's rule to solve 3x3 linear system

    REAL                 :: d
    REAL, DIMENSION(3,3) :: ai
    INTEGER              :: i

    d = determinant ( a )
    IF ( ABS(d) < tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error, determinant zero ', d
       STOP 'Error in solve'
    END IF

    DO i = 1, 3
       ai      = a
       ai(:,i) = b ! Replace i-th column of a
       x(i)    = determinant ( ai ) / d
    END DO
  END FUNCTION solve

  FUNCTION determinant ( a ) RESULT ( d )
    IMPLICIT NONE
    REAL                 :: d ! Returns determinant of
    REAL, DIMENSION(3,3) :: a ! supplied 3x3 matrix

    d =      a(1,1) * a(2,2) * a(3,3) - a(1,1) * a(2,3) * a(3,2) &
         & + a(1,2) * a(2,3) * a(3,1) - a(2,2) * a(3,1) * a(1,3) &
         & + a(2,1) * a(3,2) * a(1,3) - a(3,3) * a(1,2) * a(2,1)

  END FUNCTION determinant

  FUNCTION q_to_a ( q ) RESULT ( a )
    IMPLICIT NONE
    REAL, DIMENSION(3,3)             :: a ! Returns a 3x3 rotation matrix calculated from
    REAL, DIMENSION(0:3), INTENT(in) :: q ! supplied quaternion

    ! The rows of the rotation matrix correspond to unit vectors of the molecule in the space-fixed frame
    ! The third row  a(3,:) is "the" axis of the molecule, for uniaxial molecules
    ! Use a to convert space-fixed to body-fixed axes thus: db = matmul(a,ds)
    ! Use transpose of a to convert body-fixed to space-fixed axes thus: ds = matmul(db,a)

    ! The supplied quaternion should be normalized and we check for this

    REAL :: norm

    norm = SUM ( q**2 ) ! Quaternion squared length
    IF ( ABS ( norm - 1.0 ) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,2es20.8)' ) 'quaternion normalization error', norm, tol
       STOP 'Error in q_to_a'
    END IF

    ! Write out row by row, for clarity

    a(1,:) = [ q(0)**2+q(1)**2-q(2)**2-q(3)**2,   2*(q(1)*q(2)+q(0)*q(3)),       2*(q(1)*q(3)-q(0)*q(2))     ] ! 1st row
    a(2,:) = [     2*(q(1)*q(2)-q(0)*q(3)),   q(0)**2-q(1)**2+q(2)**2-q(3)**2,   2*(q(2)*q(3)+q(0)*q(1))     ] ! 2nd row
    a(3,:) = [     2*(q(1)*q(3)+q(0)*q(2)),       2*(q(2)*q(3)-q(0)*q(1)),   q(0)**2-q(1)**2-q(2)**2+q(3)**2 ] ! 3rd row
  END FUNCTION q_to_a

END MODULE maths_module
