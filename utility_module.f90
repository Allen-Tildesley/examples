! utility_module.f90
! routines for I/O, random numbers, averages, order parameters
MODULE utility_module

  ! We use the standard error_unit for error messages
  ! but allow the output_unit to be passed in as an argument, when needed
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit, iostat_end, iostat_eor

  IMPLICIT NONE
  PRIVATE

  ! I/O routines, including some calculation of averages
  PUBLIC :: read_cnf_atoms, write_cnf_atoms, read_cnf_mols, write_cnf_mols
  PUBLIC :: run_begin, run_end, blk_begin, blk_end, blk_add, time_stamp

  ! Random number routines
  PUBLIC :: init_random_seed, random_integer, random_normal
  PUBLIC :: random_orientation_vector, random_perpendicular_vector
  PUBLIC :: random_orientation_vector_alt1, random_orientation_vector_alt2
  PUBLIC :: random_rotate_vector, random_rotate_vector_alt1, random_rotate_vector_alt2, random_rotate_vector_alt3
  PUBLIC :: random_quaternion, random_rotate_quaternion
  PUBLIC :: metropolis

  ! Low-level mathematical routines and string operations
  PUBLIC :: rotate_vector, cross_product, outer_product, q_to_a, lowercase

  ! Order parameter calculations
  PUBLIC :: orientational_order, translational_order, nematic_order

  ! These variables (run averages etc) are only accessed within this module
  INTEGER,                                      SAVE :: nvariables
  CHARACTER(len=15), DIMENSION(:), ALLOCATABLE, SAVE :: variable_names
  REAL,              DIMENSION(:), ALLOCATABLE, SAVE :: blk_averages, run_averages, errors
  REAL,                                         SAVE :: run_norm, blk_norm

  ! Define a generic interface for the outer_product functions
  INTERFACE outer_product
     MODULE PROCEDURE outer_product_2
     MODULE PROCEDURE outer_product_3
  END INTERFACE outer_product

CONTAINS

  ! Routines associated with input/output and quantities to be averaged

  SUBROUTINE time_stamp ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit
    CHARACTER(len=8)    :: date
    CHARACTER(len=10)   :: time
    REAL                :: cpu

    CALL DATE_AND_TIME ( date, time )
    CALL CPU_TIME ( cpu )
    WRITE ( unit=output_unit, fmt='(a,t45,a4,a1,a2,a1,a2)' ) 'Date: ', date(1:4), '/', date(5:6), '/', date(7:8)
    WRITE ( unit=output_unit, fmt='(a,t47,a2,a1,a2,a1,a2)' ) 'Time: ', time(1:2), ':', time(3:4), ':', time(5:6)
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)'          ) 'CPU time: ', cpu

  END SUBROUTINE time_stamp

  SUBROUTINE read_cnf_atoms ( filename, n, box, r, v ) ! Read in atomic configuration
    CHARACTER(len=*),               INTENT(in)    :: filename
    INTEGER,                        INTENT(inout) :: n
    REAL,                           INTENT(out)   :: box
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r, v

    INTEGER :: cnf_unit, ioerr, i

    ! Open given filename, will terminate on any errors

    OPEN ( newunit=cnf_unit, file=filename, status='old', action='read', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
       STOP 'Error in read_cnf_atoms'
    END IF
    READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) n
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading n from ', filename, ioerr
       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
       STOP 'Error in read_cnf_atoms'
    END IF
    READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) box
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading box from ', filename, ioerr
       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
       STOP 'Error in read_cnf_atoms'
    END IF

    ! The first call of this routine is used just to get n and box
    ! The second call attempts to read in the atomic positions and optionally velocities
    ! Expected format is one line per atom containing either r(:,i), v(:,i) or just r(:,i)

    IF ( PRESENT ( r ) ) THEN
       IF ( n /= SIZE ( r, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
          STOP 'Error in read_cnf_atoms'
       END IF

       IF ( PRESENT ( v ) ) THEN

          IF ( n /= SIZE ( v, dim=2 ) ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
             STOP 'Error in read_cnf_atoms'
          END IF

          ! Read positions, velocities
          DO i = 1, n
             READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i), v(:,i)
             IF ( ioerr /= 0 ) THEN
                WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r, v from ', filename, ioerr
                IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                STOP 'Error in read_cnf_atoms'
             END IF
          END DO

       ELSE

          ! Read positions
          DO i = 1, n
             READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i)
             IF ( ioerr /= 0 ) THEN
                WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r from ', filename, ioerr
                IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                STOP 'Error in read_cnf_atoms'
             END IF
          END DO

       END IF

    END IF

    CLOSE ( unit=cnf_unit )

  END SUBROUTINE read_cnf_atoms

  SUBROUTINE write_cnf_atoms ( filename, n, box, r, v ) ! Write out atomic configuration
    CHARACTER(len=*),               INTENT(in) :: filename
    INTEGER,                        INTENT(in) :: n
    REAL,                           INTENT(in) :: box
    REAL, DIMENSION(:,:),           INTENT(in) :: r
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: v

    INTEGER :: cnf_unit, ioerr, i

    ! Open given filename, replacing it if it already exists, will terminate on any errors
    OPEN ( newunit=cnf_unit, file=filename, status='replace', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
       STOP 'Error in write_cnf_atoms'
    END IF
    WRITE ( unit=cnf_unit, fmt='(i15)'  ) n
    WRITE ( unit=cnf_unit, fmt='(f15.8)') box

    IF ( n /= SIZE ( r, dim=2 ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
       STOP 'Error in write_cnf_atoms'
    END IF

    IF ( PRESENT ( v ) ) THEN

       IF ( n /= SIZE ( v, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
          STOP 'Error in write_cnf_atoms'
       END IF

       ! Write positions, velocities
       DO i = 1, n
          WRITE ( unit=cnf_unit, fmt='(*(f15.10))' ) r(:,i), v(:,i)
       END DO

    ELSE

       ! Write positions
       DO i = 1, n
          WRITE ( unit=cnf_unit, fmt='(*(f15.10))' ) r(:,i)
       END DO

    END IF

    CLOSE ( unit=cnf_unit )

  END SUBROUTINE write_cnf_atoms

  SUBROUTINE read_cnf_mols ( filename, n, box, r, e, v, w ) ! Read in molecular configuration
    CHARACTER(len=*),               INTENT(in)    :: filename
    INTEGER,                        INTENT(inout) :: n
    REAL,                           INTENT(out)   :: box
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r, e, v, w

    INTEGER :: cnf_unit, ioerr, i

    ! Open given filename, will terminate on any errors

    OPEN ( newunit=cnf_unit, file=filename, status='old', action='read', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
       STOP 'Error opening file in read_cnf_mols'
    END IF
    READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) n
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading n from ', filename, ioerr
       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
       STOP 'Error in read_cnf_mols'
    END IF
    READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) box
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading box from ', filename, ioerr
       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
       STOP 'Error in read_cnf_mols'
    END IF

    ! The first call of this routine is used just to get n and box
    ! The second call attempts to read in the atomic positions, orientations and optionally velocities, angular velocities
    ! Expected format is one line per atom containing either r(:,i), e(:,i), v(:,i), w(:,i)  or just r(:,i), e(:,i)
    ! The first dimension of the e array can be 3 (vector) or 4 (quaternion)

    IF ( PRESENT ( r ) ) THEN

       IF ( .NOT. PRESENT ( e )    ) THEN
          WRITE ( unit=error_unit, fmt='(a,a,i15)') 'r and e arguments must be present together'
          STOP 'Error in read_cnf_mols'
       END IF
       IF ( n /= SIZE ( r, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
          STOP 'Error in read_cnf_mols'
       END IF
       IF ( n /= SIZE ( e, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of e ', n, SIZE ( e, dim=2 )
          STOP 'Error in read_cnf_mols'
       END IF

       IF ( PRESENT ( v ) ) THEN

          IF ( .NOT. PRESENT ( w )    ) THEN
             WRITE ( unit=error_unit, fmt='(a,a,i15)') 'v and w arguments must be present together'
             STOP 'Error in read_cnf_mols'
          END IF
          IF ( n /= SIZE ( v, dim=2 ) ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
             STOP 'Error in read_cnf_mols'
          END IF
          IF ( n /= SIZE ( w, dim=2 ) ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of w ', n, SIZE ( w, dim=2 )
             STOP 'Error in read_cnf_mols'
          END IF

          ! Read positions, orientation vectors or quaternions, velocities, angular velocities
          DO i = 1, n
             READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i), e(:,i), v(:,i), w(:,i)
             IF ( ioerr /= 0 ) THEN
                WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r, e, v, w from ', filename, ioerr
                IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                STOP 'Error in read_cnf_mols'
             END IF
          END DO

       ELSE

          ! Read positions, orientation vectors or quaternions
          DO i = 1, n
             READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i), e(:,i)
             IF ( ioerr /= 0 ) THEN
                WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r, e from ', filename, ioerr
                IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                STOP 'Error in read_cnf_mols'
             END IF
          END DO

       END IF

    END IF

    CLOSE ( unit=cnf_unit )

  END SUBROUTINE read_cnf_mols

  SUBROUTINE write_cnf_mols ( filename, n, box, r, e, v, w ) ! Write out molecular configuration
    CHARACTER(len=*),               INTENT(in) :: filename
    INTEGER,                        INTENT(in) :: n
    REAL,                           INTENT(in) :: box
    REAL, DIMENSION(:,:),           INTENT(in) :: r, e
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: v, w

    INTEGER :: cnf_unit, ioerr, i

    ! Open given filename, replacing it if it already exists, will terminate on any errors
    OPEN ( newunit=cnf_unit, file=filename, status='replace', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
       STOP 'Error in write_cnf_mols'
    END IF
    WRITE(cnf_unit,'(i15)'  ) n
    WRITE(cnf_unit,'(f15.8)') box

    IF ( n /= SIZE ( r, dim=2 ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
       STOP 'Error in write_cnf_mols'
    END IF
    IF ( n /= SIZE ( e, dim=2 ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of e ', n, SIZE ( e, dim=2 )
       STOP 'Error in write_cnf_mols'
    END IF

    IF ( PRESENT ( v ) ) THEN
       IF ( .NOT. PRESENT ( w )    ) THEN
          WRITE ( unit=error_unit, fmt='(a,a,i15)') 'v and w arguments must be present together'
          STOP 'Error in write_cnf_mols'
       END IF
       IF ( n /= SIZE ( v, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
          STOP 'Error in write_cnf_mols'
       END IF
       IF ( n /= SIZE ( w, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of w ', n, SIZE ( w, dim=2 )
          STOP 'Error in write_cnf_mols'
       END IF

       ! Write positions, orientation vectors or quaternions, velocities, angular velocities
       DO i = 1, n
          WRITE ( unit=cnf_unit, fmt='(*(f15.10))') r(:,i), e(:,i), v(:,i), w(:,i)
       END DO

    ELSE

       ! Write positions, orientation vectors or quaternions
       DO i = 1, n
          WRITE ( unit=cnf_unit, fmt='(*(f15.10))') r(:,i), e(:,i)
       END DO

    END IF

    CLOSE ( unit=cnf_unit )

  END SUBROUTINE write_cnf_mols

  SUBROUTINE run_begin ( names ) ! Set up averaging variables based on supplied names
    CHARACTER(len=15), DIMENSION(:), INTENT(in) :: names

    nvariables = SIZE ( names )
    ALLOCATE ( variable_names(nvariables) )
    ALLOCATE ( blk_averages(nvariables) )
    ALLOCATE ( run_averages(nvariables) )
    ALLOCATE ( errors(nvariables) )

    variable_names = names
    run_norm       = 0.0
    run_averages   = 0.0
    errors         = 0.0

  END SUBROUTINE run_begin

  SUBROUTINE blk_begin ! Zero averaging variables at start of each block
    blk_norm     = 0.0
    blk_averages = 0.0
  END SUBROUTINE blk_begin

  SUBROUTINE blk_add ( variables ) ! Increment block-average variables
    REAL, DIMENSION(:), INTENT(in) :: variables

    IF ( SIZE(variables) /= nvariables ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Mismatched variable arrays', nvariables, SIZE(variables)
       STOP 'Error in blk_add'
    END IF

    blk_averages = blk_averages + variables ! Increment block averages
    blk_norm     = blk_norm + 1.0           ! Increment block normalizer
  END SUBROUTINE blk_add

  SUBROUTINE blk_end ( blk, output_unit ) ! Write out block averages
    INTEGER, INTENT(in) :: blk, output_unit

    LOGICAL, SAVE :: first_call = .TRUE.

    blk_averages = blk_averages / blk_norm     ! Normalize block averages
    run_averages = run_averages + blk_averages ! Increment run averages
    errors       = errors + blk_averages**2    ! Increment error accumulators
    run_norm     = run_norm + 1.0              ! Increment run normalizer

    IF ( first_call ) THEN  ! Write headings
       WRITE ( unit=output_unit, fmt='(*(a16))'   ) REPEAT ( '=', 16*(nvariables+1) ) 
       WRITE ( unit=output_unit, fmt='(*(1x,a15))') 'Block', ADJUSTR ( variable_names )
       WRITE ( unit=output_unit, fmt='(*(a16))'   ) REPEAT ( '=', 16*(nvariables+1) )
       first_call = .FALSE.
    END IF

    ! Write out block averages
    WRITE ( unit=output_unit, fmt='(1x,i15,*(1x,f15.5))') blk, blk_averages

  END SUBROUTINE blk_end

  SUBROUTINE run_end ( output_unit ) ! Write out run averages
    INTEGER, INTENT(in) :: output_unit

    run_averages = run_averages / run_norm  ! Normalize run averages
    errors       = errors / run_norm        ! Normalize error estimates
    errors       = errors - run_averages**2 ! Compute fluctuations
    WHERE ( errors > 0.0 )
       errors = SQRT ( errors / run_norm ) ! Normalize and get estimated errors
    END WHERE

    WRITE ( unit=output_unit, fmt='(*(a16))'             ) REPEAT('-',16*(nvariables+1))
    WRITE ( unit=output_unit, fmt='(1x,a15,*(1x,f15.5))' ) 'Run averages', run_averages
    WRITE ( unit=output_unit, fmt='(1x,a15,*(1x,f15.5))' ) 'Run errors', errors
    WRITE ( unit=output_unit, fmt='(*(a16))'             ) REPEAT('=',16*(nvariables+1))

    DEALLOCATE ( variable_names, blk_averages, run_averages, errors )

  END SUBROUTINE run_end

  ! Routines associated with random number generation

  ! This routine, and the next one, are taken from the online GNU documentation
  ! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
  ! and is specific to the gfortran compiler
  ! At the time of writing, calling RANDOM_SEED() initializes the random number generator
  ! with the same random seed to a default state, which may result in the same sequence
  ! being generated every time. The routines below are intended to generate different
  ! sequences on different calls.
  ! YOU SHOULD INVESTIGATE THE BEHAVIOUR FOR YOUR OWN COMPILER AND MACHINE IMPLEMENTATION 
  SUBROUTINE init_random_seed()
    USE iso_fortran_env, ONLY: int64
    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER :: i, n, un, istat, dt(8), pid
    INTEGER(int64) :: t

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    ! First try if the OS provides a random number generator
    OPEN(newunit=un, file='/dev/urandom', access='stream', &
         form='unformatted', action='read', status='old', iostat=istat)
    IF (istat == 0) THEN
       READ(un) seed
       CLOSE(un)
    ELSE
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       CALL SYSTEM_CLOCK(t)
       IF (t == 0) THEN
          CALL DATE_AND_TIME(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       END IF
       pid = getpid()
       t = IEOR(t, INT(pid, KIND(t)))
       DO i = 1, n
          seed(i) = lcg(t)
       END DO
    END IF
    CALL RANDOM_SEED(put=seed)
  END SUBROUTINE init_random_seed

  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  FUNCTION lcg(s)
    USE iso_fortran_env, ONLY: int64
    IMPLICIT NONE
    INTEGER :: lcg
    INTEGER(int64) :: s
    IF (s == 0) THEN
       s = 104729
    ELSE
       s = MOD(s, 4294967296_int64)
    END IF
    s = MOD(s * 279470273_int64, 4294967291_int64)
    lcg = INT(MOD(s, INT(HUGE(0), int64)), KIND(0))
  END FUNCTION lcg

  FUNCTION random_integer ( k1, k2 ) RESULT ( k )
    INTEGER             :: k      ! returns uniformly distributed random integer
    INTEGER, INTENT(in) :: k1, k2 ! in range [k1,k2] inclusive

    INTEGER :: k_lo, k_hi
    REAL    :: zeta

    CALL RANDOM_NUMBER ( zeta )
    k_lo = MIN(k1,k2)
    k_hi = MAX(k1,k2)
    k =  k_lo + FLOOR((k_hi-k_lo+1)*zeta)
    IF ( k < k_lo ) k = k_lo ! guard against small danger of roundoff
    IF ( k > k_hi ) k = k_hi ! guard against small danger of roundoff

  END FUNCTION random_integer

  FUNCTION random_normal ( mean, std ) RESULT ( r )
    REAL             :: r         ! returns normal random number
    REAL, INTENT(in) :: mean, std ! with required mean and standard deviation

    ! Box-Muller transform produces numbers in pairs, we save one for next time

    REAL, DIMENSION(2)      :: zeta
    REAL,              SAVE :: r_save
    LOGICAL,           SAVE :: saved = .FALSE.
    REAL,         PARAMETER :: pi = 4.0*ATAN(1.0)

    IF ( saved ) THEN
       r = r_save
       r = mean + std * r
       saved = .FALSE.
    ELSE
       CALL RANDOM_NUMBER (zeta)
       r      = SQRT(-2*LOG(zeta(1)))*COS(2*pi*zeta(2))
       r_save = SQRT(-2*LOG(zeta(1)))*SIN(2*pi*zeta(2))
       r      = mean + std * r
       saved = .TRUE.
    END IF
  END FUNCTION random_normal

  SUBROUTINE random_orientation_vector ( e )
    REAL, DIMENSION(3), INTENT(out) :: e ! Uniformly sampled orientation

    ! Firstly, the vector is chosen uniformly within the unit cube
    ! Vectors lying outside the unit sphere are rejected
    ! Having found a vector within the unit sphere, it is normalized
    ! Essentially the same routine will work in 2d, or for quaternions in 4d

    REAL :: e_sq

    DO
       CALL RANDOM_NUMBER ( e ) ! Random numbers uniformly sampled in range (0,1)
       e    = 2.0 * e - 1.0     ! Now in range (-1,+1)
       e_sq = SUM ( e**2 )
       IF ( e_sq <= 1.0 ) EXIT
    END DO

    e = e / SQRT ( e_sq )

  END SUBROUTINE random_orientation_vector

  SUBROUTINE random_orientation_vector_alt1 ( e )
    REAL, DIMENSION(3), INTENT(out) :: e ! Uniformly sampled orientation vector

    ! First alternative routine for choosing a random orientation in 3D

    REAL               :: c, s, phi
    REAL, PARAMETER    :: pi = 4.0*ATAN(1.0)
    REAL, DIMENSION(2) :: zeta ! random numbers

    CALL RANDOM_NUMBER ( zeta )        ! Two uniformly sampled random numbers in range (0,1)
    c   = 2.0*zeta(1) - 1.0            ! Random cosine uniformly sampled in range (-1,+1)
    s   = SQRT(1.0-c**2)               ! Sine
    phi = zeta(2) * 2.0*pi             ! Random angle uniform sampled in range (0,2*pi)
    e  = [ s*COS(phi), s*SIN(phi), c ] ! Random unit vector

  END SUBROUTINE random_orientation_vector_alt1

  SUBROUTINE random_orientation_vector_alt2 ( e )
    REAL, DIMENSION(3), INTENT(out) :: e ! Uniformly sampled orientation vector

    ! Second alternative routine for choosing a random orientation in 3D

    REAL, DIMENSION(2) :: zeta
    REAL               :: zeta_sq, f

    DO
       CALL RANDOM_NUMBER ( zeta ) ! Two uniform random numbers between 0 and 1
       zeta = 2.0 * zeta - 1.0     ! now each between -1 and 1
       zeta_sq = SUM ( zeta**2 )   ! squared magnitude
       IF ( zeta_sq < 1.0 ) EXIT   ! now inside unit disk
    END DO

    f = 2.0 * SQRT ( 1.0 - zeta_sq )
    e = [ zeta(1) * f, zeta(2) * f, 1.0 - 2.0 * zeta_sq ] ! on surface of unit sphere

  END SUBROUTINE random_orientation_vector_alt2

  FUNCTION random_rotate_vector ( delta_max, e_old ) RESULT ( e )
    REAL, DIMENSION(3)             :: e         ! orientation vector result
    REAL, INTENT(in)               :: delta_max ! maximum magnitude of rotation
    REAL, DIMENSION(3), INTENT(in) :: e_old     ! original orientation

    ! Function to generate random orientation for linear molecule in 3D
    ! by rotation through a small angle from a given orientation
    ! Provided delta_max is << 1, it is approximately the maximum rotation angle (in radians)
    ! The magnitude of the rotation is not uniformly sampled, but this should not matter

    REAL, DIMENSION(3) :: de   ! random vector
    REAL               :: e_sq

    CALL random_orientation_vector ( de ) ! Random unit vector
    de  = delta_max * de                  ! Random small vector

    e    = e_old + de     ! Choose new orientation by adding random small vector
    e_sq = SUM ( e**2 )
    e    = e / SQRT(e_sq) ! Normalize

  END FUNCTION random_rotate_vector

  FUNCTION random_perpendicular_vector ( e_old ) RESULT ( e )
    REAL, DIMENSION(3)             :: e         ! result is an orientation vector
    REAL, DIMENSION(3), INTENT(in) :: e_old     ! perpendicular to this vector 

    REAL            :: factor, e_sq
    REAL, PARAMETER :: tol = 1.e-6

    DO
       CALL random_orientation_vector ( e )                 ! Random unit vector
       factor = dot_PRODUCT ( e, e_old ) / SUM ( e_old**2 ) ! Projection along e_old
       e      = e - factor * e_old                          ! Make e perpendicular to e_old
       e_sq   = SUM ( e**2 )
       IF ( e_sq > tol ) EXIT ! Start again if e is too small
    END DO
    e = e / SQRT ( e_sq ) ! Random unit direction perpendicular to e_old
  END FUNCTION random_perpendicular_vector

  FUNCTION random_rotate_vector_alt1 ( delta_max, e_old ) RESULT ( e )
    REAL, DIMENSION(3)             :: e         ! orientation vector result
    REAL, INTENT(in)               :: delta_max ! maximum angle of rotation
    REAL, DIMENSION(3), INTENT(in) :: e_old     ! original orientation

    ! First alternative function to generate random orientation for linear molecule in 3D
    ! by rotation through a small angle from a given orientation
    ! delta_max is the maximum rotation angle (in radians)
    ! The magnitude of the rotation is uniformly sampled

    REAL, DIMENSION(3) :: e_perp
    REAL               :: e_sq, delta, zeta

    e_perp = random_perpendicular_vector ( e_old ) ! Choose unit vector perpendicular to e_old

    CALL RANDOM_NUMBER ( zeta ) ! Random number uniformly sampled in range (0,1)
    zeta  = 2.0 * zeta - 1.0    ! Now in range (-1,+1)
    delta = zeta * delta_max    ! Random rotation angle

    e    = e_old * COS ( delta ) + e_perp * SIN ( delta )
    e_sq = SUM ( e**2 )
    e    = e / SQRT(e_sq) ! Normalize

  END FUNCTION random_rotate_vector_alt1

  FUNCTION random_rotate_vector_alt2 ( delta_max, e_old ) RESULT ( e )
    REAL, DIMENSION(3)             :: e         ! orientation vector result
    REAL, INTENT(in)               :: delta_max ! maximum magnitude of rotation
    REAL, DIMENSION(3), INTENT(in) :: e_old     ! original orientation

    ! Second alternative function to generate random orientation in 3D without any angular bias
    ! by rotation through a small angle from a given orientation
    ! delta_max is the maximum rotation angle (in radians)
    ! The magnitude of the rotation is uniformly sampled
    ! The rotation axis is a Cartesian axis selected at random
    ! Ref: Barker and Watts, Chem Phys Lett 3, 144 (1969)

    INTEGER            :: k
    REAL, DIMENSION(3) :: axis  ! axis of rotation
    REAL               :: delta ! rotation angle
    REAL               :: zeta  ! random number

    k       = random_integer (1,3)! random axis choice 1 = x, 2 = y, 3 = z
    axis    = 0.0
    axis(k) = 1.0

    CALL RANDOM_NUMBER ( zeta )              ! uniform random number between 0 and 1
    delta = ( 2.0 * zeta - 1.0 ) * delta_max ! uniform random angle

    e = rotate_vector ( delta, axis, e_old )

  END FUNCTION random_rotate_vector_alt2

  FUNCTION random_rotate_vector_alt3 ( delta_max, e_old ) RESULT ( e )
    REAL, DIMENSION(3)             :: e         ! orientation vector result
    REAL, INTENT(in)               :: delta_max ! maximum magnitude of rotation
    REAL, DIMENSION(3), INTENT(in) :: e_old     ! original orientation

    ! Ref: Marsaglia, Ann Maths Stat 43, 645 (1972)
    ! Uses a rejection technique to create a trial orientation
    ! subject to the constraint that the cosine of the angle
    ! turned through is greater than cos(delta_max)

    REAL :: cos_min

    cos_min = COS ( delta_max )

    DO
       CALL random_orientation_vector_alt2 ( e )
       IF ( DOT_PRODUCT ( e, e_old ) > cos_min ) EXIT ! close enough
    END DO

  END FUNCTION random_rotate_vector_alt3

  SUBROUTINE random_quaternion ( e )
    REAL, DIMENSION(0:3), INTENT(out) :: e ! Uniformly sampled orientation quaternion

    REAL, DIMENSION(2) :: zeta
    REAL               :: s1, s2, f

    DO
       CALL RANDOM_NUMBER ( zeta ) ! Two uniform random numbers between 0 and 1
       zeta = 2.0 * zeta - 1.0     ! now each between -1 and 1
       s1 = SUM ( zeta**2 )        ! squared magnitude
       IF ( s1 < 1.0 ) EXIT        ! now inside unit disk
    END DO
    e(0) = zeta(1)
    e(1) = zeta(2)

    DO
       CALL RANDOM_NUMBER ( zeta ) ! Two uniform random numbers between 0 and 1
       zeta = 2.0 * zeta - 1.0     ! now each between -1 and 1
       s2 = SUM ( zeta**2 )        ! squared magnitude
       IF ( s2 < 1.0 ) EXIT        ! now inside unit disk
    END DO
    f = SQRT ( (1.0 - s1 ) / s2 )
    e(2) = zeta(1)*f
    e(3) = zeta(2)*f

  END SUBROUTINE random_quaternion

  FUNCTION random_rotate_quaternion ( delta_max, q_old ) RESULT ( q )
    REAL,                 INTENT(in) :: delta_max ! maximum rotation angle
    REAL, DIMENSION(0:3), INTENT(in) :: q_old     ! old quaternion
    REAL, DIMENSION(0:3)             :: q         ! result

    REAL, DIMENSION(3)   :: axis        ! rotation axis
    REAL                 :: delta, s, c ! rotation angle, sin, cos
    REAL, DIMENSION(0:3) :: q_rot       ! rotation quaternion

    CALL random_orientation_vector ( axis )
    CALL random_NUMBER ( delta )
    delta = ( 2.0*delta - 1.0 ) * delta_max
    s = SIN(0.5*delta)
    c = COS(0.5*delta)
    q_rot = [ c, s*axis(1), s*axis(2), s*axis(3) ]
    q = quatmul ( q_rot, q_old )

  END FUNCTION random_rotate_quaternion

  FUNCTION metropolis ( delta ) ! Conduct Metropolis test, with safeguards
    LOGICAL          :: metropolis
    REAL, INTENT(in) :: delta

    REAL            :: zeta
    REAL, PARAMETER :: exponent_guard = 75.0

    IF ( delta > exponent_guard ) THEN ! too high, reject without evaluating
       metropolis = .FALSE.
    ELSE IF ( delta < 0.0 ) THEN ! downhill, accept without evaluating
       metropolis = .TRUE.
    ELSE
       CALL RANDOM_NUMBER ( zeta )     ! Uniform random number in range (0,1)
       metropolis = EXP(-delta) > zeta ! Metropolis test
    END IF

  END FUNCTION metropolis

  ! Low level mathematical operations and string manipulation

  FUNCTION rotate_vector ( delta, axis, e_old ) RESULT ( e )
    REAL, DIMENSION(3)             :: e           ! orientation vector result
    REAL, INTENT(in)               :: delta       ! rotation angle
    REAL, DIMENSION(3), INTENT(in) :: axis, e_old ! rotation axis and original orientation

    REAL :: dot, c, s

    c   = COS ( delta )
    s   = SIN ( delta )
    dot = DOT_PRODUCT ( axis, e_old )

    e = c * e_old + (1.0-c)*dot*axis + s * cross_product ( axis, e_old )

  END FUNCTION rotate_vector

  FUNCTION quatmul ( a, b ) RESULT ( c ) ! multiply two quaternions
    REAL, DIMENSION(0:3), INTENT(in) :: a, b ! arguments
    REAL, DIMENSION(0:3)             :: c    ! result

    c(0) = a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
    c(1) = a(1)*b(0) + a(0)*b(1) - a(3)*b(2) + a(2)*b(3)
    c(2) = a(2)*b(0) + a(3)*b(1) + a(0)*b(2) - a(1)*b(3)
    c(3) = a(3)*b(0) - a(2)*b(1) + a(1)*b(2) + a(0)*b(3)

  END FUNCTION quatmul

  FUNCTION cross_product ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3)             :: c    ! result cross product
    REAL, DIMENSION(3), INTENT(in) :: a, b ! arguments
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  END FUNCTION cross_product

  FUNCTION outer_product_2 ( a, b ) RESULT (c)
    REAL, DIMENSION(:), INTENT(IN)   :: a, b
    REAL, DIMENSION(SIZE(a),SIZE(b)) :: c ! function result

    INTEGER :: i, j

    DO i = 1, SIZE(a)
       DO j = 1, SIZE(b)
          c(i,j) = a(i) * b(j)
       END DO
    END DO

    ! The following one-line statement is equivalent, but the above loops are clearer
    ! c = SPREAD(a,dim=2,ncopies=SIZE(b)) * SPREAD(b,dim=1,ncopies=SIZE(a))

  END FUNCTION outer_product_2

  FUNCTION outer_product_3 ( a, b, c ) RESULT (d)
    REAL, DIMENSION(:), INTENT(IN)           :: a, b, c
    REAL, DIMENSION(SIZE(a),SIZE(b),size(c)) :: d ! function result

    INTEGER :: i, j, k

    DO i = 1, SIZE(a)
       DO j = 1, SIZE(b)
          do k = 1, size(c)
             d(i,j,k) = a(i) * b(j) * c(k)
             end do
       END DO
    END DO

  END FUNCTION outer_product_3

  FUNCTION lowercase ( oldstring ) RESULT ( newstring )
    IMPLICIT NONE
    CHARACTER(len=*),             INTENT(in)    :: oldstring 
    CHARACTER(len=LEN(oldstring))               :: newstring 

    INTEGER :: i, k 

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
    REAL                                          :: order ! result order parameter
    REAL,    DIMENSION(:,:), INTENT(in)           :: r     ! set of molecular position vectors (3,n)
    INTEGER, DIMENSION(3),   INTENT(in), OPTIONAL :: k     ! Lattice reciprocal vector (integer)

    ! Calculate the "melting factor" for translational order 
    ! based on a single k-vector characterizing the original lattice
    ! and commensurate with the periodic box
    ! It is assumed that both r and k are in box=1 units
    ! k = (l,m,n) where l,m,n are integers
    ! If optional argument k is omitted, we default to a choice
    ! based on the fcc lattice, if this makes sense
    ! order = 1 when all atoms are on their lattice positions
    ! order = 1/sqrt(n), approximately, for disordered positions

    INTEGER            :: i, n, nc
    REAL, DIMENSION(3) :: k_real
    REAL               :: kr
    COMPLEX            :: rho ! Fourier component of single-particle density
    REAL, PARAMETER    :: pi = 4.0*ATAN(1.0), twopi = 2.0*pi

    IF ( SIZE(r,dim=1) /= 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error in r dimension ', SIZE(r,dim=1)
       STOP 'Error in translational_order'
    END IF
    n = SIZE(r,dim=2)

    IF ( PRESENT ( k ) ) THEN
       k_real = twopi * REAL ( k )
    ELSE                                          ! Make arbitrary choice assuming fcc
       nc = NINT ( ( REAL(n)/4.0 ) ** (1.0/3.0) ) ! number of fcc unit cells
       IF ( 4*nc**3 /= n ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Error in value of n ', 4*nc**3, n
          STOP 'Error in translational_order'
       END IF
       k_real = twopi * REAL( [-nc,nc,-nc] )      ! arbitrary fcc reciprocal vector
    END IF

    rho = ( 0.0, 0.0 )

    DO i = 1, n
       kr  = dot_PRODUCT ( k_real, r(:,i) )
       rho = rho + CMPLX ( COS(kr), SIN(kr) )
    END DO

    rho = rho / REAL(n)
    order = REAL ( CONJG(rho)*rho )

  END FUNCTION translational_order

  FUNCTION orientational_order ( e ) RESULT ( order )
    REAL                             :: order ! result order parameter
    REAL, DIMENSION(:,:), INTENT(in) :: e     ! set of molecular orientation vectors (3,n)

    ! Calculates an orientational order parameter to monitor "melting"
    ! The parameter depends completely on knowing the orientations of the molecules
    ! in the original crystal lattice, and here we assume a specific alpha-fcc crystal
    ! of the same kind as was set up in initialize.f90 and initialize_module.f90
    ! Four molecules per unit cell, each pointing along a body-diagonal
    ! Order parameter can be a low-ranking (e.g. 1st or 2nd) Legendre polynomial

    INTEGER :: n, nc, i, i0
    REAL    :: c

    REAL, PARAMETER :: rroot3 = 1.0 / SQRT ( 3.0 )

    REAL, DIMENSION(3,4), PARAMETER :: e0 = RESHAPE (  rroot3*[ &
         &  1.0,  1.0,  1.0,    1.0, -1.0, -1.0,  &
         & -1.0,  1.0, -1.0,   -1.0, -1.0,  1.0 ],[3,4] ) ! orientations in unit cell

    IF ( SIZE(e,dim=1) /= 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error in e dimension ', SIZE(e,dim=1)
       STOP 'Error in orientational_order'
    END IF
    n = SIZE(e,dim=2)
    nc = NINT ( ( REAL(n)/4.0 ) ** (1.0/3.0) )
    IF ( 4*nc**3 /= n ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Error in value of n ', 4*nc**3, n
       STOP 'Error in orientational_order'
    END IF
    order = 0.0
    DO i = 1, n
       i0 = MODULO ( i, 4 ) + 1             ! select appropriate original orientation
       c = dot_PRODUCT ( e(:,i), e0(:,i0) ) ! cosine of angle
       order = order + 1.5*c**2 - 0.5       ! Second Legendre polynomial
    END DO
    order = order / REAL ( n )

  END FUNCTION orientational_order

  FUNCTION nematic_order ( e ) RESULT ( order )
    REAL                             :: order
    REAL, DIMENSION(:,:), INTENT(in) :: e     ! set of molecular orientation vectors (3,n)

    ! Calculate the nematic order parameter <P2(cos(theta))>
    ! where theta is the angle between a molecular axis and the director
    ! which is the direction that maximises the order parameter
    ! This is obtained by finding the largest eigenvalue of
    ! the 3x3 second-rank traceless order tensor

    INTEGER              :: i, n
    REAL, DIMENSION(3,3) :: q         ! order tensor
    REAL                 :: h, g, psi ! used in eigenvalue calculation
    REAL, PARAMETER      :: pi = 4.0*ATAN(1.0)

    IF ( SIZE(e,dim=1) /= 3 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error in e dimension ', SIZE(e,dim=1)
       STOP 'Error in nematic_order'
    END IF
    n = SIZE(e,dim=2)

    ! Order tensor: outer product of each orientation vector, summed over molecules
    q = SUM ( SPREAD ( e, dim=2, ncopies=3) * SPREAD ( e, dim=1, ncopies=3 ), dim = 3 )
    q = 1.5 * q / REAL(n)                ! normalize
    FORALL (i=1:3) q(i,i) = q(i,i) - 0.5 ! make traceless

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

  FUNCTION q_to_a ( q ) RESULT ( a ) ! converts quaternion to rotation matrix
    IMPLICIT NONE

    ! Arguments
    REAL, DIMENSION(0:3), INTENT(in) :: q ! quaternion
    REAL, DIMENSION(3,3)             :: a ! rotation matrix

    ! The rows of the rotation matrix correspond to unit vectors of the molecule in the space-fixed frame
    ! The third row  a(3,:) is [2*(q(1)*q(3)+q(0)*q(2)),2*(q(2)*q(3)-q(0)*q(1)),q(0)**2-q(1)**2-q(2)**2+q(3)**2]
    ! which is "the" axis of the molecule, for uniaxial molecules
    ! use a to convert space-fixed to body-fixed axes thus: db = matmul(a,ds)
    ! use transpose of a to convert body-fixed to space-fixed axes thus: ds = matmul(db,a)

    a(1,:) = [ q(0)**2+q(1)**2-q(2)**2-q(3)**2,   2*(q(1)*q(2)+q(0)*q(3)),       2*(q(1)*q(3)-q(0)*q(2))     ] ! 1st row
    a(2,:) = [     2*(q(1)*q(2)-q(0)*q(3)),   q(0)**2-q(1)**2+q(2)**2-q(3)**2,   2*(q(2)*q(3)+q(0)*q(1))     ] ! 2nd row
    a(3,:) = [     2*(q(1)*q(3)+q(0)*q(2)),       2*(q(2)*q(3)-q(0)*q(1)),   q(0)**2-q(1)**2-q(2)**2+q(3)**2 ] ! 3rd row
  END FUNCTION q_to_a

END MODULE utility_module
