! error_calc.f90
! Estimated error in correlated data
PROGRAM error_calc

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

  ! In this program we analyse synthetic data, a time series whose properties are exactly known.
  ! Define underlying process by generalized Langevin equation (GLE)
  ! with memory function expressed as a decaying exponential
  ! See G Ciccotti and JP Ryckaert Mol Phys 40 141 (1980)
  ! and AD Baczewski and SD Bond J Chem Phys 139 044107 (2013)

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               maths_module,    ONLY : random_normal, init_random_seed

  IMPLICIT NONE

  REAL    :: m      ! GLE memory function coefficients
  REAL    :: kappa  ! GLE memory function decay rates
  REAL    :: s      ! GLE auxiliary variable
  REAL    :: delta  ! Time step
  REAL    :: at     ! Dynamical variable at time t
  INTEGER :: nstep  ! Number of timesteps in run
  INTEGER :: nequil ! Number of equilibration timesteps
  INTEGER :: t      ! Time (equivalent to step number in file)

  REAL, DIMENSION(:), ALLOCATABLE :: a ! Stored data values (nstep)

  INTEGER :: ioerr, tblock, i_repeat, nblock, blk, stp1, stp2, trun, n_repeat
  REAL    :: average, variance, stddev, err, si, tcor, x, e, b, d
  REAL    :: a_blk, a_run, a_avg, a_var, a_var_1, a_err

  ! Taylor series coefficients
  REAL, PARAMETER :: b1 = 2.0, b2 = -2.0,     b3 = 4.0/3.0, b4 = -2.0/3.0  ! b = 1.0 - EXP(-2.0*x)
  REAL, PARAMETER :: d1 = 1.0, d2 = -1.0/2.0, d3 = 1.0/6.0, d4 = -1.0/24.0 ! d = 1.0 - EXP(-x)

  NAMELIST /nml/ nstep, nequil, n_repeat, delta, variance, average

  ! Example default values
  nstep    = 2**16 ! Number of steps, about 60,000 for example
  nequil   = 10000 ! Number of equilibration timesteps
  n_repeat = 50    ! Number of simulation repeats for brute force empirical calculation
  delta    = 0.01  ! Timestep for simulation
  variance = 1.0   ! Desired variance of data (equivalent to temperature)
  average  = 1.0   ! Desired average value of data

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in error_calc'
  END IF

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps in run',        nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Equilibration steps',           nequil
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of brute-force repeats', n_repeat
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step delta',               delta
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Desired average value',         average
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Desired variance',              variance

  ALLOCATE ( a(nstep) )

  ! The memory function model is defined here
  ! Values used by Baczewski and Bond in their example are
  ! (m,kappa)
  ! (1.0,1.0)  (underdamped)
  ! (0.5,2.0)  (critically damped)
  ! (0.25,4.0) (overdamped)
  m     = 0.25
  kappa = 4.0
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'm ',     m
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'kappa ', kappa

  ! Coefficients used in algorithm
  x = delta*kappa
  e = EXP(-x) ! theta in B&B paper
  IF ( x > 0.0001 ) THEN
     b = 1.0 - EXP(-2.0*x)
     d = 1.0 - EXP(-x)
  ELSE ! Taylor expansions for low x
     b = x * ( b1 + x * ( b2 + x * ( b3 + x * b4 ) ) ) 
     d = x * ( d1 + x * ( d2 + x * ( d3 + x * d4 ) ) )
  END IF
  b      = SQRT ( b )
  b      = b * sqrt ( kappa/2.0 ) ! alpha in B&B paper  
  stddev = SQRT(2.0*variance)     ! NB stddev of random forces, not data

  ! For this process, the results of interest can be calculated exactly
  ! The time correlation function is known, and hence the statistical inefficiency (SI)
  ! From this, the error in the mean for a run of any length can be calculated

  tcor = 1.0 / m      ! Correlation time is determined by memory function
  tcor = tcor / delta ! Express this in timesteps
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

  a_avg = 0.0 ! Zero average accumulator
  a_var = 0.0 ! Zero mean-squared accumulator

  DO i_repeat = 1, n_repeat ! Loop over repeats of simulation

     ! Initial values, hopefully not critical
     at = 0.0
     s  = 0.0

     DO t = -nequil, nstep ! Loop over steps including an equilibration period

        ! Velocity Verlet type algorithm for at and auxiliary variable s
        at = at + 0.5 * delta * s
        s  = e * s - d * m * at + b * SQRT(m) * random_normal ( 0.0, stddev )
        at = at + 0.5 * delta * s

        IF ( t > 0 ) a(t) = average + at ! Store values (after equilibration, adding average)

     END DO ! End loop over steps including an equilibration period

     a_run = SUM(a) / REAL(nstep) ! The run average
     a_avg = a_avg + a_run        ! Average over runs
     a_var = a_var + a_run**2     ! Mean squared value over runs

  END DO ! End loop over repeats of simulation

  a_avg = a_avg / REAL(n_repeat)                  ! Mean value
  a_var = a_var / REAL(n_repeat)                  ! Mean-squared value
  a_var = a_var - a_avg**2                        ! Mean-squared deviation
  a_var = a_var * REAL(n_repeat)/REAL(n_repeat-1) ! Bias correction
  err   = SQRT(a_var)                             ! Empirical standard deviation in run average
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Empirical error estimate = ', err
  WRITE ( unit=output_unit, fmt='(a)') 'This should be in reasonable agreement with exact error estimate'

  ! Now analyse the last run, as if it were the only one that had been carried out
  ! This is what usually happens; we rarely have the luxury of many independent runs

  ! Simple calculation of average and variance

  a_avg   = SUM ( a ) / REAL(nstep)      ! Sample average
  a       = a - a_avg                    ! Centre the data
  a_var_1 = SUM ( a**2 ) / REAL(nstep-1) ! Bias-corrected sample variance
  a_err   = SQRT(a_var_1 / REAL(nstep))  ! Error estimate neglecting any correlations

  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Sample average value = ',                   a_avg
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Deviation from exact average = ',           a_avg - average
  WRITE ( unit=output_unit, fmt='(a)') 'Deviation should (typically) lie within +/- exact error estimate'
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Sample variance = ',                        a_var_1
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Error estimate neglecting correlations = ', a_err
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

  DO nblock = 20, 4, -1 ! Loop over number, and hence length, of blocks

     tblock = nstep / nblock              ! Block length in steps (rounded down)
     trun   = nblock*tblock               ! Run length in steps, accounting for rounding
     a_run  = SUM(a(1:trun)) / REAL(trun) ! Average of data
     a_var  = 0.0                         ! Zero mean-square block average accumulator

     DO blk = 1, nblock ! Loop over blocks
        stp1  = (blk-1)*tblock+1                            ! Start of block
        stp2  = blk*tblock                                  ! End of block
        a_blk = SUM ( a(stp1:stp2) - a_run ) / REAL(tblock) ! Block average of deviation
        a_var = a_var + a_blk**2                            ! Mean-square block average
     END DO ! End loop over blocks

     a_var = a_var / REAL(nblock-1)    ! Bias-corrected variance of block averages
     a_err = SQRT(a_var/REAL(nblock))  ! Estimate of error from block-average variance
     si    = tblock * a_var / a_var_1  ! Statistical inefficiency

     WRITE ( unit=output_unit, fmt='(2i15,3f15.6)' ) tblock, nblock, a_err, si

  END DO ! End loop over number, and hence length, of blocks

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

  DO ! Loop over number, and hence length, of blocks

     nblock = nblock / 2 ! Halve the number of blocks, rounding down if nblock is odd
     tblock = tblock*2   ! Double the block length
     IF ( nblock < 3 ) EXIT

     a(1:nblock) = ( a(1:2*nblock-1:2) + a(2:2*nblock:2) ) / 2.0 ! Blocking transformation, halving the data set
     a_avg       = SUM ( a(1:nblock) ) / REAL(nblock)            ! Re-compute sample average
     a(1:nblock) = a(1:nblock) - a_avg                           ! Re-centre in case of dropped points (odd nblock)
     a_var       = SUM ( a(1:nblock)**2 ) / REAL(nblock-1)       ! Bias-corrected variance of block averages
     a_err       = SQRT ( a_var / REAL(nblock) )                 ! Estimate of error from block average variance
     si          = tblock * a_var / a_var_1                      ! Statistical inefficiency

     WRITE ( unit=output_unit, fmt='(2i15,3f15.6)' ) tblock, nblock, a_err, si

  END DO ! End loop over number, and hence length, of blocks

  WRITE ( unit=output_unit, fmt='(a)' ) 'Plateau at large tblock (small nblock)'
  WRITE ( unit=output_unit, fmt='(a)' ) 'should agree quite well with exact error estimate'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Can plot SI or error**2 against 1/tblock or log2(tblock)'

END PROGRAM error_calc

