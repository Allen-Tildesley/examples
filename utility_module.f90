MODULE utility_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: metropolis, read_cnf_atoms, write_cnf_atoms, read_cnf_molecules, write_cnf_molecules
  PUBLIC :: run_begin, run_end, blk_begin, blk_end, blk_add
  PUBLIC :: random_integer, orientational_order, random_rotate

  INTEGER,                                      SAVE :: nvariables
  CHARACTER(len=10), DIMENSION(:), ALLOCATABLE, SAVE :: variable_names
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
       CALL random_NUMBER ( zeta )     ! Uniform random number in range (0,1)
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
    CHARACTER(len=10), DIMENSION(:), INTENT(in) :: names

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
       WRITE(*,'(*(a15))') REPEAT ( '=', 15*(nvariables+1) ) 
       WRITE(*,'(*(5x,a10))') 'Block     ', variable_names
       WRITE(*,'(*(a15))') REPEAT ( '=', 15*(nvariables+1) )
       first_call = .FALSE.
    END IF

    ! Write out block averages
    WRITE(*,'(5x,i10,*(5x,f10.4))') blk, blk_averages

  END SUBROUTINE blk_end

  SUBROUTINE run_end

    run_averages = run_averages / run_norm  ! Normalize run averages
    errors       = errors / run_norm        ! Normalize error estimates
    errors       = errors - run_averages**2 ! Compute fluctuations
    WHERE ( errors > 0.0 )
       errors = SQRT ( errors / run_norm ) ! Normalize and get estimated errors
    END WHERE

    WRITE(*,'(*(a15))') REPEAT('-',15*(nvariables+1))
    WRITE(*,'(a15,*(5x,f10.4))') 'Run averages', run_averages
    WRITE(*,'(a15,*(5x,f10.4))') 'Run errors', errors
    WRITE(*,'(*(a15))') REPEAT('=',15*(nvariables+1))

  END SUBROUTINE run_end

  FUNCTION random_integer ( k1, k2 ) RESULT ( k ) ! returns random integer in range [k1,k2] inclusive
    INTEGER             :: k
    INTEGER, INTENT(in) :: k1, k2

    INTEGER :: k_lo, k_hi
    REAL    :: zeta

    CALL random_NUMBER ( zeta )
    k_lo = MIN(k1,k2)
    k_hi = MAX(k1,k2)
    k =  k_lo + FLOOR((k_hi-k_lo+1)*zeta)
    IF ( k < k_lo ) k = k_lo ! guard against small danger of roundoff
    IF ( k > k_hi ) k = k_hi ! guard against small danger of roundoff

  END FUNCTION random_integer

  FUNCTION orientational_order ( e ) RESULT ( order )
    REAL                             :: order
    REAL, DIMENSION(:,:), INTENT(in) :: e

    INTEGER              :: i, n
    REAL, DIMENSION(3,3) :: q
    REAL                 :: h, g, psi
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

  FUNCTION random_rotate ( de_max, e_old ) RESULT ( e )
    REAL, DIMENSION(3)             :: e      ! orientation vector result
    REAL, INTENT(in)               :: de_max ! maximum magnitude of rotation
    REAL, DIMENSION(3), INTENT(in) :: e_old  ! original orientation

    REAL, DIMENSION(2) :: zeta ! random numbers
    REAL, DIMENSION(3) :: de   ! random vector
    REAL               :: c, s, phi, norm
    REAL, PARAMETER    :: pi = 4.0*ATAN(1.0)

    CALL RANDOM_NUMBER ( zeta )          ! Two uniformly sampled random numbers in range (0,1)
    c   = 2.0*zeta(1) - 1.0             ! Random cosine uniformly sampled in range (-1,+1)
    s   = SQRT(1.0-c**2)                ! Sine
    phi = zeta(2) * 2.0*pi              ! Random angle uniform sampled in range (0,2*pi)
    de  = [ s*COS(phi), s*SIN(phi), c ] ! Random unit vector
    de  = de_max * de                   ! Random vector

    e    = e_old + de      ! Choose new orientation by adding random small vector
    norm = SQRT(SUM(e**2)) ! Normalizing factor
    e    = e / norm        ! Normalize

  END FUNCTION random_rotate

END MODULE utility_module
