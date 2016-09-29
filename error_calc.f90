! error_calc.f90
! Estimated error in correlated data
PROGRAM error_calc
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE maths_module, ONLY : random_normal, init_random_seed

  IMPLICIT NONE

  ! Define underlying process by generalized Langevin equation
  ! with memory function expressed as a decaying exponential
  ! see G Ciccotti and JP Ryckaert Mol Phys 40 141 (1980)
  ! and AD Baczewski and SD Bond J Chem Phys 139 044107 (2013)

  REAL    :: m      ! memory function coefficients
  REAL    :: kappa  ! memory function decay rates
  REAL    :: theta  ! auxiliary coefficient
  REAL    :: alpha  ! auxiliary coefficient
  REAL    :: zeta   ! random number
  REAL    :: s      ! GLE auxiliary variable
  REAL    :: delta  ! time step
  REAL    :: at     ! dynamical variable at time t
  INTEGER :: nstep  ! number of timesteps in run
  INTEGER :: nequil ! number of equilibration timesteps
  INTEGER :: t      ! time (equivalent to step number in file)

  REAL, DIMENSION(:), ALLOCATABLE :: a ! stored data values (nstep)

  INTEGER :: ioerr, tblock, i_repeat, nblock, blk, stp1, stp2, trun
  REAL    :: average, variance, stddev, err, si, tcor, a_val, a_avg, a_var, a_var_1, a_err

  INTEGER, PARAMETER :: n_repeat = 50 ! number of simulation repeats for brute force empirical calculation

  NAMELIST /nml/ nstep, nequil, delta, variance, average

  ! Example default values
  nstep    = 2**16 ! number of steps, about 60,000 for example
  nequil   = 10000 ! number of equilibration timesteps
  delta    = 0.01  ! timestep for simulation
  variance = 1.0   ! desired variance of data (equivalent to temperature)
  average  = 1.0   ! desired average value of data

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in error_calc'
  END IF

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'number of steps in run = ', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'equilibration steps = ',    nequil
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'time step delta = ',        delta
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'desired average value = ',  average
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'desired variance = ',       variance
  stddev = SQRT(2.0*variance) ! NB stddev of random forces, not data

  ALLOCATE ( a(nstep) )

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

  ! For this process, the results of interest can be calculated exactly
  ! The time correlation function is known, and hence the statistical inefficiency (SI)
  ! From this, the error in the mean for a run of any length can be calculated

  tcor = 1.0 / m      ! Correlation time is determined by memory function
  tcor = tcor / delta ! express in timesteps
  si   = 2*tcor       ! Statistical inefficiency (SI)
  err  = SQRT(si*variance/REAL(nstep))
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Exact correlation time in steps = ', tcor
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Run length / correlation time = ',   REAL(nstep) / tcor
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Exact value of SI = ',               si
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Exact error estimate = ',            err

  ! Data generation
  CALL init_random_seed

  ! For comparison, we do n_repeat independent runs and estimate the error in run averages directly from these
  ! This is to give an empirical idea of the distribution from which the run average is sampled
  ! Of course, we expect the results to be similar to the exact values just calculated

  a_avg = 0.0 ! zero average accumulator
  a_var = 0.0 ! zero mean-squared accumulator

  DO i_repeat = 1, n_repeat ! loop over repeats of simulation

     ! Initial values
     at = 0.0
     s  = 0.0

     DO t = -nequil, nstep ! include an equilibration period
        at   = at + 0.5*delta*s
        zeta = random_normal ( 0.0, stddev )
        s    = theta*s - (1.0-theta)*m*at + alpha*SQRT(m)*zeta
        at   = at + 0.5*delta*s
        IF ( t >= 1 ) a(t) = average + at ! store values (adding average)
     END DO

     a_val = SUM(a) / REAL(nstep) ! the run average
     a_avg = a_avg + a_val        ! average over runs
     a_var = a_var + a_val**2     ! mean squared value over runs

  END DO ! end loop over repeats of simulation

  a_avg = a_avg / REAL(n_repeat)                  ! mean value
  a_var = a_var / REAL(n_repeat)                  ! mean-squared value
  a_var = a_var - a_avg**2                        ! mean-squared deviation
  a_var = a_var * REAL(n_repeat)/REAL(n_repeat-1) ! bias correction
  err   = SQRT(a_var)                             ! empirical standard deviation in run average
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Empirical error estimate = ', err
  WRITE ( unit=output_unit, fmt='(a)') 'This should be in reasonable agreement with exact error estimate'

  ! Now analyse the last run, as if it were the only one that had been carried out
  ! This is what usually happens; we rarely have the luxury of many independent runs

  ! Simple calculation of average and variance
  a_avg   = SUM ( a ) / REAL(nstep)      ! sample average
  a       = a - a_avg                    ! centre the data
  a_var_1 = SUM ( a**2 ) / REAL(nstep-1) ! bias-corrected sample variance
  a_err   = SQRT(a_var_1 / REAL(nstep))  ! error estimate neglecting any correlations
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Sample average value = ',                   a_avg
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Deviation from exact average = ',           a_avg-average
  WRITE ( unit=output_unit, fmt='(a)') 'Deviation should (typically) lie within +/- exact error estimate'
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Sample variance = ',                        a_var_1
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Error estimate neglecting correlations = ', a_err
  WRITE ( unit=output_unit, fmt='(a)') 'This should be very over-optimistic!'

  ! We must take account of the correlations between successive values in time
  ! The two common methods which follow are actually very similar
  ! They differ in the choice of block lengths

  ! Traditional block analysis
  ! The rationale here is that 20 (say) independent block averages should be enough to get a reasonable
  ! estimate of the desired quantities, and that there is not a lot to gain by looking at more (shorter) blocks.
  ! Some workers just assume that the run may be divided into 10 or 20 blocks, which they hope will be independent.
  ! This is exactly what we do in our example programs, just for simplicity.
  ! We cannot recommend this in general, unless there is good reason to support the assumption of independence.
  ! If the 20 blocks are not independent, then attention should be focused on fewer (longer) blocks,
  ! rather than more (shorter) ones, and a plot of squared error estimate, or statistical inefficiency,
  ! vs 1/tblock carried out to extrapolate to tblock=nstep. The loop below provides the data for that plot.

  WRITE ( unit=output_unit, fmt='(4a15)' ) 'tblock', 'nblock', 'error estimate', 'estimate of SI'
  DO nblock = 20, 4, -1 ! loop over number, and hence length, of blocks
     tblock = nstep / nblock              ! block length in steps (rounded down)
     trun   = nblock*tblock               ! run length in steps, accounting for rounding
     a_avg  = SUM(a(1:trun)) / REAL(trun) ! average of data
     a_var  = 0.0                         ! zero mean-square block average accumulator
     DO blk = 1, nblock ! loop over blocks
        stp1  = (blk-1)*tblock+1                            ! start of block
        stp2  = blk*tblock                                  ! end of block
        a_val = SUM ( a(stp1:stp2) - a_avg ) / REAL(tblock) ! block average
        a_var = a_var + a_val**2                            ! mean-square block average
     END DO ! end loop over blocks
     a_var = a_var / REAL(nblock-1)    ! bias-corrected variance of block averages
     a_err = SQRT(a_var/REAL(nblock))  ! estimate of error from block-average variance
     si    = tblock * a_var / a_var_1  ! statistical inefficiency
     WRITE ( unit=output_unit, fmt='(2i15,3f15.6)' ) tblock, nblock, a_err, si
  END DO ! end loop over number, and hence length, of blocks
  WRITE ( unit=output_unit, fmt='(a)' ) 'Plateau at large tblock (small nblock)'
  WRITE ( unit=output_unit, fmt='(a)' ) 'should agree quite well with exact error estimate'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Can plot SI or error**2 against 1/tblock'

  ! Method of Flyvbjerg and Petersen, J Chem Phys, 91, 461 (1989)
  ! Note that, in this implementation, we over-write the data array a
  ! Advantages of this method are the very neat method for reducing the data
  ! and the formal scaling analysis that accompanies the blocking transformation
  ! The basic calculation is the same as for the traditional block-averaging,
  ! but the block length, and hence tblock, change by a factor 2 each time.
  ! Advocates of Flyvbjerg and Petersen might argue that the additional data points
  ! at low nblock (high tblock) which are calculated in the traditional method
  ! do not actually provide much new (independent) information.
  ! One may attempt an extrapolation of si or a_err**2 as a function of 1/tblock
  ! but F&P suggest estimating a plateau in a plot vs number of blocking transformations
  ! which is essentially log2(tblock)
  nblock = nstep
  tblock = 1
  WRITE ( unit=output_unit, fmt='(4a15)' ) 'tblock', 'nblock', 'error estimate', 'estimate of SI'
  DO ! loop over number, and hence length, of blocks
     nblock = nblock / 2 ! halve the number of blocks, rounding down if nblock is odd
     tblock = tblock*2   ! double the block length
     IF ( nblock < 3 ) EXIT
     a(1:nblock) = ( a(1:2*nblock-1:2) + a(2:2*nblock:2) ) / 2.0 ! blocking transformation, halving the data set
     a_avg       = SUM ( a(1:nblock) ) / REAL(nblock)            ! re-compute sample average
     a(1:nblock) = a(1:nblock) - a_avg                           ! re-centre in case of dropped points (odd nblock)
     a_var       = SUM ( a(1:nblock)**2 ) / REAL(nblock-1)       ! bias-corrected variance of block averages
     a_err       = SQRT ( a_var / REAL(nblock) )                 ! estimate of error from block average variance
     si          = tblock * a_var / a_var_1                      ! statistical inefficiency
     WRITE ( unit=output_unit, fmt='(2i15,3f15.6)' ) tblock, nblock, a_err, si
  END DO ! end loop over number, and hence length, of blocks
  WRITE ( unit=output_unit, fmt='(a)' ) 'Plateau at large tblock (small nblock)'
  WRITE ( unit=output_unit, fmt='(a)' ) 'should agree quite well with exact error estimate'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Can plot SI or error**2 against 1/tblock'

END PROGRAM error_calc

