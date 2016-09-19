! error.f90
! Estimated error in correlated data
PROGRAM error
  !
  ! TODO MPA provide code
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE maths_module, ONLY : random_normal, init_random_seed

  IMPLICIT NONE

  ! Define underlying process by generalized Langevin equation
  ! with memory function expressed as a decaying exponential
  ! see G Ciccotti and JP Ryckaert Mol Phys 40 141 (1980)
  ! and AD Baczewski and SD Bond J Chem Phys 139 044107 (2013)

  REAL :: m     ! memory function coefficients
  REAL :: kappa ! memory function decay rates
  REAL :: theta ! auxiliary coefficient
  REAL :: alpha ! auxiliary coefficient
  REAL :: zeta  ! random number
  REAL :: s     ! GLE auxiliary variable
  REAL :: delta ! time step
  REAL :: vt    ! velocity at time t

  REAL, DIMENSION(:), ALLOCATABLE :: v     ! stored velocities (nstep)
  REAL, DIMENSION(:), ALLOCATABLE :: a     ! working copy of stored velocities (nstep)

  INTEGER :: nstep           ! number of timesteps in run
  INTEGER :: nequil          ! number of equilibration timesteps
  INTEGER :: t               ! time (equivalent to step number in file)
  INTEGER :: ioerr, dt, n, i_repeat
  REAL    :: average, variance, stddev, err, si, tcor, a_val, a_avg, a_var, a_var_1, a_err

  INTEGER, PARAMETER :: n_repeat = 50 ! number of simulation repeats for brute force empirical calculation

  NAMELIST /nml/ nstep, nequil, delta, variance, average

  ! Example default values
  nstep           = 2**16 ! number of steps, about 60,000 for example
  nequil          = 10000 ! number of equilibration timesteps
  delta           = 0.01  ! timestep for simulation
  variance        = 1.0   ! desired variance of data (equivalent to temperature)
  average         = 1.0   ! desired average value of data

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in error'
  END IF

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'number of steps in run = ', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'equilibration steps = ',    nequil
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'time step delta = ',        delta
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'desired average value = ',  average
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'desired variance = ',       variance
  stddev = SQRT(2.0*variance) ! NB stddev of random forces, not data

  ALLOCATE ( v(nstep), a(nstep) )

  ! The memory function model is defined here
  ! Values used by Baczewski and Bond in their example are
  ! (m,kappa)
  ! (1.0,1.0)  (underdamped)
  ! (0.5,2.0)  (critically damped)
  ! (0.25,4.0) (overdamped)
  m     = 0.25
  kappa = 4.0
  theta = EXP(-delta*kappa)
  alpha = SQRT( (1.0-theta**2)*kappa/2.0 )

  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'm ',     m
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'kappa ', kappa
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'theta ', theta
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'alpha ', alpha

  ! Give exactly-known results for this process
  tcor = 1.0/m        ! Correlation time is determined by memory function
  tcor = tcor / delta ! express in timesteps
  si   = 2*tcor       ! Statistical inefficiency
  err  = SQRT(si*variance/REAL(nstep))
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Exact correlation time in steps = ', tcor
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Run length / correlation time = ',   REAL(nstep) / tcor
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Exact value of SI = ',               si
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Exact error estimate = ',            err

  ! Data generation
  CALL init_random_seed

  ! For comparison, we do n_repeat independent runs and calculate
  ! the variation in run averages directly from these
  ! This is to give an idea of the distribution from which the run average is sampled
  
  a_avg = 0.0
  a_var = 0.0

  DO i_repeat = 1, n_repeat ! loop over repeats of simulation

     ! Initial values
     vt = 0.0
     s  = 0.0

     DO t = -nequil, nstep ! include an equilibration period
        vt = vt + 0.5*delta*s
        zeta = random_normal ( 0.0, stddev )
        s  = theta*s - (1.0-theta)*m*vt + alpha*SQRT(m)*zeta
        vt = vt + 0.5*delta*s
        IF ( t >= 1 ) v(t) = average + vt ! store velocities (adding average)
     END DO

     a_val = SUM(v) / REAL(nstep) ! the run average
     a_avg = a_avg + a_val
     a_var = a_var + a_val**2

  END DO ! end loop over repeats of simulation

  a_avg = a_avg / REAL(n_repeat)
  a_var = a_var / REAL(n_repeat)
  a_var = a_var - a_avg**2
  a_var = a_var * REAL(n_repeat)/REAL(n_repeat-1) ! bias correction
  err   = SQRT(a_var) ! empirical standard deviation in run average
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Empirical error estimate = ', err
  WRITE ( unit=output_unit, fmt='(a)') 'this should be in reasonable agreement with exact error estimate'

  ! Now analyse the last run, as if it were the only one that had been carried out
  
  ! Simple calculation of average and variance
  a       = v                            ! copy of data
  a_avg   = SUM ( a ) / REAL(nstep)      ! sample average
  a       = a - a_avg                    ! centre the data
  a_var_1 = SUM ( a**2 ) / REAL(nstep-1) ! bias-corrected sample variance
  a_err   = SQRT(a_var_1 / REAL(nstep))  ! error estimate neglecting any correlations
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'sample average value = ',                   a_avg
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'deviation from exact average = ',           a_avg-average
  WRITE ( unit=output_unit, fmt='(a)') 'deviation should (typically) be roughly within +/- exact error estimate'
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'sample variance = ',                        a_var_1
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'error estimate neglecting correlations = ', a_err
  WRITE ( unit=output_unit, fmt='(a)') 'This should be very over-optimistic'

  ! Method of Flyvbjerg and Petersen, J Chem Phys, 91, 461 (1989)
  n = nstep
  dt = 1
  WRITE ( unit=output_unit, fmt='(4a15)'        ) 'dt', 'n', 'error estimate', 'SI'
  WRITE ( unit=output_unit, fmt='(2i15,2f15.6)' ) dt, n, a_err, 1.0
  DO
     n = n / 2
     dt = dt*2
     IF ( n < 3 ) EXIT
     a(1:n) = ( a(1:2*n-1:2) + a(2:2*n:2) ) / 2.0 ! blocking transformation, halving the data set
     a_avg  = SUM(a(1:n))/REAL(n)      ! re-compute sample average
     a(1:n) = a(1:n) - a_avg           ! re-centre in case of dropped points (odd n)
     a_var  = SUM(a(1:n)**2)/REAL(n-1) ! bias-corrected variance of block averages
     a_err  = SQRT(a_var/REAL(n))      ! estimate of error from block average variance
     si     = dt*a_var/a_var_1         ! statistical inefficiency
     WRITE ( unit=output_unit, fmt='(2i15,3f15.6)' ) dt, n, a_err, si
  END DO
  WRITE ( unit=output_unit, fmt='(a)') 'Plateau at larger dt / smaller n should agree quite well with exact error estimate'
  WRITE ( unit=output_unit, fmt='(a)') 'Can plot against 1/dt and extrapolate to dt=infinity'

END PROGRAM error

