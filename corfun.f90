! corfun.f90
! Time correlation function, directly and by FFT
PROGRAM corfun
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE, INTRINSIC :: iso_c_binding

  USE maths_module, ONLY : random_normal

  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

  ! Define underlying process by generalized Langevin equation
  ! with memory function expressed as a decaying exponential function
  ! see G Ciccotti and JP Ryckaert Mol Phys 40 141 (1980)
  ! and AD Baczewski and SD Bond J Chem Phys 139 044107 (2013)

  ! We assume that the program is linked with the FFTW library (version 3)
  ! but it could easily be adapted to use a different library
  ! With older FFT routines it was necessary to transform a number of points = exact power of 2
  ! which, for us, would mean nstep = 2**j (for integer j) and fft_len=2*nstep
  ! This is no longer true, so we do not restrict nstep in this way,
  ! but the transform should be more efficient if nstep happens to take such values
  ! Advantage can be taken of the fact that the data is real, but for clarity
  ! we just use the complex FFT with imaginary parts of the data set to zero

  REAL :: m     ! memory function coefficients
  REAL :: kappa ! memory function decay rates
  REAL :: zeta  ! random numbers
  REAL :: s     ! GLE auxiliary variables
  REAL :: delta ! time step
  REAL :: vt    ! velocity at time t

  REAL,    DIMENSION(:), ALLOCATABLE :: v     ! stored velocities (nstep)
  REAL,    DIMENSION(:), ALLOCATABLE :: v0    ! stored velocity origins (n0)
  INTEGER, DIMENSION(:), ALLOCATABLE :: t0    ! times of origins (n0)
  REAL,    DIMENSION(:), ALLOCATABLE :: c     ! velocity correlation function (direct method) (0:nt)
  REAL,    DIMENSION(:), ALLOCATABLE :: c_fft ! velocity correlation function (FFT method) (0:nt)
  REAL,    DIMENSION(:), ALLOCATABLE :: n     ! normalizing function (0:nt)

  INTEGER :: nt              ! number of timesteps to correlate
  INTEGER :: nstep           ! number of timesteps in run
  INTEGER :: nequil          ! number of equilibration timesteps
  INTEGER :: n0              ! number of time origins to store
  INTEGER :: origin_interval ! interval for time origins
  INTEGER :: dt              ! time difference
  INTEGER :: t               ! time (equivalent to step number in file)
  LOGICAL :: full
  INTEGER :: k, mk, nk
  INTEGER :: unit, ioerr
  REAL    :: temperature, stddev, cpu_1, cpu_2, cpu_3, cpu_4, x, e, b, d

  INTEGER(C_INT)                                       :: fft_len  ! the number of points for FFT
  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: fft_inp  ! data to be transformed (0:fft_len-1)
  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: fft_out  ! data to be transformed (0:fft_len-1)
  TYPE(C_PTR)                                          :: fft_plan ! plan needed for FFTW

  REAL,    PARAMETER :: b1 = 2.0, b2 = -2.0,     b3 = 4.0/3.0, b4 = -2.0/3.0
  REAL,    PARAMETER :: d1 = 1.0, d2 = -1.0/2.0, d3 = 1.0/6.0, d4 = -1.0/24.0

  NAMELIST /nml/ nt, origin_interval, nstep, nequil, delta, temperature

  ! Example default values
  ! Agreement (to numerical precision) of direct and FFT methods is expected if origin_interval=1
  nt              = 1000  ! max time for correlation function
  origin_interval = 1     ! This could reasonably be increased to 10 or 20 to improve efficiency
  nstep           = 2**20 ! number of steps, about a million for example
  nequil          = 10000 ! number of equilibration timesteps
  delta           = 0.01  ! timestep for simulation
  temperature     = 1.0   ! temperature for simulation

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in corfun'
  END IF

  n0 = nt / origin_interval + 1
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'number of steps in run = ',    nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'equilibration steps = ',       nequil
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'max correlation time nt = ',   nt
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'origin interval = ',           origin_interval
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'number of time origins n0 = ', n0
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'time step delta = ',           delta
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'temperature = ',               temperature

  ALLOCATE ( v(nstep), v0(n0), t0(n0) )
  ALLOCATE ( c(0:nt), c_fft(0:nt), n(0:nt) )
  fft_len = 2*nstep ! actual length of data
  ALLOCATE ( fft_inp(0:fft_len-1), fft_out(0:fft_len-1) )

  ! The memory function model is defined here
  ! Values used by Baczewski and Bond in their example are
  ! (m,kappa)
  ! (1.0,1.0)  (underdamped)
  ! (0.5,2.0)  (critically damped)
  ! (0.25,4.0) (overdamped)
  m     = 1.0
  kappa = 1.0
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.5))' ) 'm ',     m
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.5))' ) 'kappa ', kappa

  ! Coefficients used in algorithm
  x = delta*kappa
  e = EXP(-x) ! theta in B&B paper
  IF ( x > 0.0001 ) THEN
     b = 1.0 - EXP(-2.0*x)
     d = 1.0 - EXP(-x)
  ELSE  ! Taylor expansions for low x
     b = x * ( b1 + x * ( b2 + x * ( b3 + x * b4 )) )
     d = x * ( d1 + x * ( d2 + x * ( d3 + x * d4 )) )
  END IF
  b      = SQRT ( b )
  b      = b * SQRT ( kappa/2.0 ) ! alpha in B&B paper
  stddev = SQRT(2.0*temperature)

  ! Data generation
  CALL CPU_TIME ( cpu_1 )

  ! Initial values
  vt = 0.0
  s  = 0.0

  DO t = -nequil, nstep ! include an equilibration period
     vt = vt + 0.5*delta*s
     zeta = random_normal ( 0.0, stddev )
     s    = e*s - d*m*vt + b*SQRT(m)*zeta
     vt   = vt + 0.5*delta*s
     IF ( t > 0 ) v(t) = vt ! store velocities
  END DO

  CALL CPU_TIME ( cpu_2 )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'CPU time to generate data = ', cpu_2-cpu_1

  ! Data analysis (direct method)

  c(:) = 0.0
  n(:) = 0.0
  mk   = 0 ! storage location of time origin
  full = .FALSE.

  DO t = 1, nstep ! Main loop correlating data

     IF ( MOD(t-1,origin_interval) == 0 ) THEN
        mk = mk + 1
        IF ( mk > n0 ) THEN
           full = .TRUE.
           mk = mk - n0 ! overwrite older values
        END IF
        t0(mk) = t    ! store time origins
        v0(mk) = v(t) ! store velocity at time origins
     END IF

     IF ( full ) THEN
        nk = n0 ! correlate with all stored time origins
     ELSE
        nk = mk ! correlate with those stored so far
     END IF

     DO k = 1, nk ! Loop over time origins
        dt = t - t0(k)
        IF ( dt >= 0 .AND. dt <= nt ) THEN
           c(dt) = c(dt) + v(t) * v0(k) ! increment correlation function
           n(dt) = n(dt) + 1.0          ! increment normalizing factor
        END IF
     END DO ! End loop over time origins

  END DO ! End main loop correlating data
  c(:) = c(:) / n(:) ! Normalise by number of increments

  CALL CPU_TIME ( cpu_3 )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'CPU time for direct method = ', cpu_3-cpu_2

  ! Data analysis (FFT method)

  ! Prepare data for FFT
  fft_inp            = CMPLX(0.0,0.0) ! fill input array with zeros
  fft_inp(0:nstep-1) = CMPLX(v)       ! put data into first part (real only)

  ! Forward FFT
  fft_plan = fftw_plan_dft_1d ( fft_len, fft_inp, fft_out, FFTW_FORWARD, FFTW_ESTIMATE) ! set up plan
  CALL fftw_execute_dft ( fft_plan, fft_inp, fft_out )                                  ! execute FFT
  CALL fftw_destroy_plan ( fft_plan )                                                   ! release plan

  fft_out = fft_out * CONJG ( fft_out ) ! square modulus

  ! Reverse FFT
  fft_plan = fftw_plan_dft_1d ( fft_len, fft_out, fft_inp, FFTW_BACKWARD, FFTW_ESTIMATE) ! set up plan
  CALL fftw_execute_dft ( fft_plan, fft_out, fft_inp )                                   ! execute FFT
  CALL fftw_destroy_plan ( fft_plan )                                                    ! release plan

  fft_inp = fft_inp / REAL ( fft_len ) ! normalization factor associated with FFT itself
  n(:)    = [ ( nstep-t, t = 0, nt ) ]
  c_fft   = REAL ( fft_inp(0:nt) ) / n(:) ! normalization associated with time origins

  CALL CPU_TIME ( cpu_4 )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'CPU time for FFT method = ', cpu_4-cpu_3

  WRITE ( unit=output_unit, fmt='(a)' ) 'Output to corfun.out'

  OPEN ( newunit=unit, file='corfun.out', status='replace', iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error opening file', ioerr
     STOP 'Error in corfun'
  END IF

  DO t = 0, nt
     WRITE ( unit=unit, fmt='(i10,3f15.8)' ) t, c(t), c_fft(t), c_exact(t*delta)
  END DO

  CLOSE(unit=unit)

  DEALLOCATE ( v, v0, t0, c, c_fft, n )
  DEALLOCATE ( fft_inp, fft_out )

CONTAINS

  FUNCTION c_exact ( t ) RESULT ( c )
    IMPLICIT NONE

    REAL, INTENT(in) :: t ! time argument
    REAL             :: c ! result

    REAL :: del, omega, kp, km, cp, cm

    ! See AD Baczewski and SD Bond J Chem Phys 139 044107 (2013)
    ! In general the exact correlation function may be obtained from the inverse Laplace transform
    ! C(s) = (kT/m) * 1 / ( s + M(s) ) where both C(t) and M(t) are sums of exponentials in t
    ! M(s) = sum_p m_p*kappa_p / (s+kappa_p) p = 1..Np (note amplitude mp*kappa_p)
    ! C(s) = sum_p c_p / (s+k_p) p = 1..Np+1
    ! The coefficients c_p and k_p may be determined in terms of the m_p and kappa_p
    ! by solving an equation of order Np+1 and using partial fractions
    ! Here we just do this for the simplest case Np=1

    IF ( kappa > 4.0 * m ) THEN  ! Real roots

       del = SQRT ( 0.25*kappa**2 - m*kappa )
       kp  = 0.5*kappa + del
       km  = 0.5*kappa - del
       cp  = -km/(kp-km)
       cm  =  kp/(kp-km)
       c   = cp*EXP(-kp*t) + cm*EXP(-km*t)

    ELSE ! Complex roots

       omega = SQRT ( m*kappa - 0.25*kappa**2 )
       c     = EXP(-0.5*kappa*t) * ( COS(omega*t) + (0.5*kappa/omega)*SIN(omega*t) )

    END IF

  END FUNCTION c_exact

END PROGRAM corfun

