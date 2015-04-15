MODULE utility_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: metropolis
  PUBLIC :: read_cnf_atoms, write_cnf_atoms, read_cnf_molecules, write_cnf_molecules
  PUBLIC :: run_begin, run_end, blk_begin, blk_end, blk_add
  PUBLIC :: random_integer, random_normal
  PUBLIC :: random_orientation_vector, random_perpendicular_vector
  PUBLIC :: random_orientation_vector_alt1, random_orientation_vector_alt2
  PUBLIC :: random_rotate_vector, random_rotate_vector_alt1, random_rotate_vector_alt2, random_rotate_vector_alt3
  PUBLIC :: orientational_order

  INTEGER,                                      SAVE :: nvariables
  CHARACTER(len=15), DIMENSION(:), ALLOCATABLE, SAVE :: variable_names
  REAL,              DIMENSION(:), ALLOCATABLE, SAVE :: blk_averages, run_averages, errors
  REAL,                                         SAVE :: run_norm, blk_norm

CONTAINS

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

  SUBROUTINE read_cnf_atoms ( filename, n, box, r, v ) ! Read in atomic configuration
    CHARACTER(len=*),               INTENT(in)    :: filename
    INTEGER,                        INTENT(inout) :: n
    REAL,                           INTENT(out)   :: box
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r, v

    INTEGER :: cnf_unit, ioerr, i

    OPEN(newunit=cnf_unit,file=filename,status='old',action='read',iostat=ioerr)
    IF ( ioerr /= 0 ) STOP 'Error opening file in read_cnf_atoms'
    READ(cnf_unit,*) n
    READ(cnf_unit,*) box

    IF ( PRESENT ( r ) ) THEN
       IF ( n /= SIZE ( r, dim=2 ) ) STOP 'r size wrong in read_cnf_atoms'

       IF ( PRESENT ( v ) ) THEN
          IF ( n /= SIZE ( v, dim=2 ) ) STOP 'v size wrong in read_cnf_atoms'

          ! Read positions, velocities
          DO i = 1, n
             READ(cnf_unit,*) r(:,i), v(:,i)
          END DO

       ELSE

          ! Read positions
          DO i = 1, n
             READ(cnf_unit,*) r(:,i)
          END DO

       END IF

    END IF

    CLOSE(unit=cnf_unit)

  END SUBROUTINE read_cnf_atoms

  SUBROUTINE write_cnf_atoms ( filename, n, box, r, v )
    CHARACTER(len=*),               INTENT(in) :: filename
    INTEGER,                        INTENT(in) :: n
    REAL,                           INTENT(in) :: box
    REAL, DIMENSION(:,:),           INTENT(in) :: r
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: v

    INTEGER :: cnf_unit, ioerr, i

    OPEN(newunit=cnf_unit,file=filename,status='replace',iostat=ioerr)
    IF ( ioerr /= 0 ) STOP 'Error opening file in write_cnf_atoms'
    WRITE(cnf_unit,'(i15)'  ) n
    WRITE(cnf_unit,'(f15.8)') box

    IF ( n /= SIZE ( r, dim=2 ) ) STOP 'r size wrong in write_cnf_atoms'

    IF ( PRESENT ( v ) ) THEN

       IF ( n /= SIZE ( v, dim=2 ) ) STOP 'v size wrong in write_cnf_atoms'

       ! Write positions, velocities
       DO i = 1, n
          WRITE(cnf_unit,'(*(f15.10))') r(:,i), v(:,i)
       END DO

    ELSE

       ! Write positions
       DO i = 1, n
          WRITE(cnf_unit,'(*(f15.10))') r(:,i)
       END DO

    END IF

    CLOSE(unit=cnf_unit)

  END SUBROUTINE write_cnf_atoms

  SUBROUTINE read_cnf_molecules ( filename, n, box, r, e, v, w ) ! Read in molecular configuration
    CHARACTER(len=*),               INTENT(in)    :: filename
    INTEGER,                        INTENT(inout) :: n
    REAL,                           INTENT(out)   :: box
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r, e, v, w

    INTEGER :: cnf_unit, ioerr, i

    OPEN(newunit=cnf_unit,file=filename,status='old',action='read',iostat=ioerr)
    IF ( ioerr /= 0 ) STOP 'Error opening file in read_cnf_molecules'
    READ(cnf_unit,*) n
    READ(cnf_unit,*) box

    IF ( PRESENT ( r ) ) THEN

       IF ( .NOT. PRESENT ( e )    ) STOP 'Argument e missing in read_cnf_molecules'
       IF ( n /= SIZE ( r, dim=2 ) ) STOP 'r size wrong in read_cnf_molecules'
       IF ( n /= SIZE ( e, dim=2 ) ) STOP 'e size wrong in read_cnf_molecules'

       IF ( PRESENT ( v ) ) THEN

          IF ( .NOT. PRESENT ( w )    ) STOP 'Argument w missing in read_cnf_molecules'
          IF ( n /= SIZE ( v, dim=2 ) ) STOP 'v size wrong in read_cnf_molecules'
          IF ( n /= SIZE ( w, dim=2 ) ) STOP 'w size wrong in read_cnf_molecules'

          ! Read positions, orientation vectors or quaternions, velocities, angular velocities
          DO i = 1, n
             READ(cnf_unit,*) r(:,i), e(:,i), v(:,i), w(:,i)
          END DO

       ELSE

          ! Read positions, orientation vectors or quaternions
          DO i = 1, n
             READ(cnf_unit,*) r(:,i), e(:,i)
          END DO

       END IF

    END IF

    CLOSE(unit=cnf_unit)

  END SUBROUTINE read_cnf_molecules

  SUBROUTINE write_cnf_molecules ( filename, n, box, r, e, v, w )
    CHARACTER(len=*),               INTENT(in) :: filename
    INTEGER,                        INTENT(in) :: n
    REAL,                           INTENT(in) :: box
    REAL, DIMENSION(:,:),           INTENT(in) :: r, e
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: v, w

    INTEGER :: cnf_unit, ioerr, i

    OPEN(newunit=cnf_unit,file=filename,status='replace',iostat=ioerr)
    IF ( ioerr /= 0 ) STOP 'Error opening file in write_cnf_molecules'
    WRITE(cnf_unit,'(i15)'  ) n
    WRITE(cnf_unit,'(f15.8)') box

    IF ( n /= SIZE ( r, dim=2 ) ) STOP 'r size wrong in write_cnf_molecules'
    IF ( n /= SIZE ( e, dim=2 ) ) STOP 'e size wrong in write_cnf_molecules'

    IF ( PRESENT ( v ) ) THEN
       IF ( .NOT. PRESENT ( w )    ) STOP 'Argument w missing in write_cnf_molecules'
       IF ( n /= SIZE ( v, dim=2 ) ) STOP 'v size wrong in write_cnf_molecules'
       IF ( n /= SIZE ( w, dim=2 ) ) STOP 'w size wrong in write_cnf_molecules'

       ! Write positions, orientation vectors or quaternions, velocities, angular velocities
       DO i = 1, n
          WRITE(cnf_unit,'(*(f15.10))') r(:,i), e(:,i), v(:,i), w(:,i) ! positions and velocities
       END DO

    ELSE

       ! Write positions, orientation vectors or quaternions
       DO i = 1, n
          WRITE(cnf_unit,'(*(f15.10))') r(:,i), e(:,i)
       END DO

    END IF

    CLOSE(unit=cnf_unit)

  END SUBROUTINE write_cnf_molecules

  SUBROUTINE run_begin ( names )
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

  SUBROUTINE blk_begin
    blk_norm     = 0.0
    blk_averages = 0.0
  END SUBROUTINE blk_begin

  SUBROUTINE blk_add ( variables )
    REAL, DIMENSION(:), INTENT(in) :: variables

    IF ( SIZE(variables) /= nvariables ) STOP 'mismatched variable arrays in stp_end'
    blk_averages = blk_averages + variables ! Increment block averages
    blk_norm     = blk_norm + 1.0           ! Increment block normalizer
  END SUBROUTINE blk_add

  SUBROUTINE blk_end ( blk )
    INTEGER, INTENT(in) :: blk

    LOGICAL, SAVE :: first_call = .TRUE.

    blk_averages = blk_averages / blk_norm     ! Normalize block averages
    run_averages = run_averages + blk_averages ! Increment run averages
    errors       = errors + blk_averages**2    ! Increment error accumulators
    run_norm     = run_norm + 1.0              ! Increment run normalizer

    IF ( first_call ) THEN  ! Write headings
       WRITE(*,'(*(a16))') REPEAT ( '=', 16*(nvariables+1) ) 
       WRITE(*,'(*(1x,a15))') 'Block', variable_names
       WRITE(*,'(*(a16))') REPEAT ( '=', 16*(nvariables+1) )
       first_call = .FALSE.
    END IF

    ! Write out block averages
    WRITE(*,'(1x,i15,*(1x,f15.5))') blk, blk_averages

  END SUBROUTINE blk_end

  SUBROUTINE run_end

    run_averages = run_averages / run_norm  ! Normalize run averages
    errors       = errors / run_norm        ! Normalize error estimates
    errors       = errors - run_averages**2 ! Compute fluctuations
    WHERE ( errors > 0.0 )
       errors = SQRT ( errors / run_norm ) ! Normalize and get estimated errors
    END WHERE

    WRITE(*,'(*(a16))') REPEAT('-',16*(nvariables+1))
    WRITE(*,'(1x,a15,*(1x,f15.5))') 'Run averages', run_averages
    WRITE(*,'(1x,a15,*(1x,f15.5))') 'Run errors', errors
    WRITE(*,'(*(a16))') REPEAT('=',16*(nvariables+1))

    DEALLOCATE ( variable_names, blk_averages, run_averages, errors )

  END SUBROUTINE run_end

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
       CALL random_NUMBER ( e ) ! Random numbers uniformly sampled in range (0,1)
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

    CALL random_NUMBER ( zeta ) ! Random number uniformly sampled in range (0,1)
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

    INTEGER :: axis        ! axis of rotation
    REAL    :: delta, c, s ! rotation angle
    REAL    :: zeta        ! random number

    axis  = random_integer (1,3)             ! random axis choice 1 = x, 2 = y, 3 = z
    CALL RANDOM_NUMBER ( zeta )              ! uniform random number between 0 and 1
    delta = ( 2.0 * zeta - 1.0 ) * delta_max ! uniform random angle
    c     = COS ( delta )
    s     = SIN ( delta )

    SELECT CASE ( axis )

    CASE ( 1 ) ! rotate about x
       e(2) = c * e_old(2) + s * e_old(3)
       e(3) = c * e_old(3) - s * e_old(2)

    CASE ( 2 ) ! rotate about y
       e(3) = c * e_old(3) + s * e_old(1)
       e(1) = c * e_old(1) - s * e_old(3)

    CASE default ! rotate about z
       e(1) = c * e_old(1) + s * e_old(2)
       e(2) = c * e_old(2) - s * e_old(1)

    END SELECT

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

  FUNCTION orientational_order ( e ) RESULT ( order )
    REAL                             :: order
    REAL, DIMENSION(:,:), INTENT(in) :: e     ! set of molecular orientation vectors (3,n)

    ! Calculate the nematic order parameter <P2(cos(theta))>
    ! Where theta is the angle between a molecular axis and the director
    ! This is obtained by finding the largest eigenvalue of
    ! the 3x3 second-rank traceless order tensor

    INTEGER              :: i, n
    REAL, DIMENSION(3,3) :: q         ! order tensor
    REAL                 :: h, g, psi ! used in eigenvalue calculation
    REAL, PARAMETER      :: pi = 4.0*ATAN(1.0)

    IF ( SIZE(e,dim=1) /= 3 ) STOP 'Array error in orientational_order'
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

  END FUNCTION orientational_order

END MODULE utility_module
