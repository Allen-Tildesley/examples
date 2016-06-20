! qmc_walk_sho.f90
! Quantum Monte Carlo, random walk, simple harmonic oscillator
PROGRAM qmc_walk_sho

  USE utility_module, ONLY: random_normal
  IMPLICIT NONE

  ! Program to calculate the ground state wavefunction
  ! for a particle in a harmonic potential, V=(x**2)/2,
  ! by solving the corresponding diffusion equation in imaginary time.

  ! In atomic units, mass=1, hbar=1, diffusion coefficient D = hbar**2/(2*m) = 1/2,
  ! so root-mean-squared displacement in time step ds is sqrt(2*D*ds) = sqrt(ds)
  ! Walkers are created and destroyed, depending on the potential energy relative to the trial energy et
  ! The program adjusts et to attempt to converge on a target number of walkers
  ! hence obtaining the exact ground state value et = 0.5
  ! This type of simulation is very sensitive to the initial guess for et,
  ! and the resulting time evolution is quite noisy:
  ! results are output as averages over the production period.
  ! The simulated ground state wave function may be compared with the exact result for this simple problem.

  ! The default parameters run the simulation with et = 0.5, the exact groundstate energy for this potential.
  ! You can then try, say, et = 0.6 and et = 0.4, observing the behaviour of the number of walkers in each case.

  ! Many of the parameters, and the updating scheme for et, are taken from the following paper:
  ! I Kostin, B Faber and K Schulten, Amer J Phys, 64, 633 (1996)

  REAL,    DIMENSION(:), ALLOCATABLE  :: x       ! position of each walker (n_max)
  REAL,    DIMENSION(:), ALLOCATABLE  :: v       ! potential energy of each walker (n_max)
  INTEGER, DIMENSION(:), ALLOCATABLE  :: replica ! number of replicas to make (n_max)
  LOGICAL, DIMENSION(:), ALLOCATABLE  :: alive   ! flags those walkers still alive (n_max)
  INTEGER, DIMENSION(:), ALLOCATABLE  :: bin     ! histogram bins for wavefunction (n_bin)
  REAL,    DIMENSION(:), ALLOCATABLE  :: psi     ! wavefunction (n_bin)

  REAL, PARAMETER :: pi = 4.0 * ATAN ( 1.0 ), const = (1.0/pi) ** (1.0/4.0)

  INTEGER :: n, n_max, n_target, n_add, n_bin, nlive, ibin
  INTEGER :: i, production_steps, equilibration_steps, output_interval, step, step_count
  REAL    :: x_max, de, ds, et, k, et_avg, pot_avg, pot
  REAL    :: zeta, dx
  REAL    :: factor, pexact, x_bin

  NAMELIST /params/ production_steps, equilibration_steps, output_interval, x_max, et, ds, n_bin, n_max, n_target

  WRITE(*,'(''qmc_walk_sho'')')
  WRITE(*,'(''Diffusion Monte Carlo simulation of a quantum oscillator'')')
  WRITE(*,'(''Results in atomic units'')')

  ! Set up default parameters for the simulation
  n_max               = 2000  ! max number of walkers
  n_target            = 500   ! target number of walkers
  production_steps    = 20000 ! number of steps for production
  equilibration_steps = 20000 ! number of steps for equilibration
  output_interval     = 200   ! output interval
  et                  = 0.5   ! initial estimate of ground state energy in atomic units
  ds                  = 0.1   ! timestep in imaginary time
  x_max               = 10.0  ! half length of the line, centred at x = 0.0 
  n_bin               = 200   ! number of bins covering -x_max .. x_max
  READ(*,nml=params)
  WRITE(*,'(a,t40,i15)'  ) 'Max number of walkers = ',                   n_max
  WRITE(*,'(a,t40,i15)'  ) 'Target number of walkers = ',                n_target
  WRITE(*,'(a,t40,i15)'  ) 'Number of steps for production = ',          production_steps
  WRITE(*,'(a,t40,i15)'  ) 'Number of steps for equilibration = ',       equilibration_steps
  WRITE(*,'(a,t40,i15)'  ) 'Output interval = ',                         output_interval
  WRITE(*,'(a,t40,f15.5)') 'Initial estimate of ground state energy = ', et
  WRITE(*,'(a,t40,f15.5)') 'Time step = ',                               ds
  WRITE(*,'(a,t40,f15.5)') 'Max x for wavefunction = ',                  x_max
  WRITE(*,'(a,t40,i15)'  ) 'Number of bins for wavefunction = ',         n_bin

  dx = 2.0*x_max / REAL(n_bin)  ! Bin interval for calculating the wavefunction
  IF ( n_target > n_max ) STOP 'n_target exceeds n_max at start'

  ALLOCATE ( x(n_max), v(n_max), replica(n_max), alive(n_max) )
  ALLOCATE ( psi(n_bin), bin(n_bin) )

  CALL RANDOM_SEED()

  ! Set up an initial delta distribution of walkers at origin
  n = n_target
  x(1:n) = 0.0

  ! Set up energy averages and array for calculating the wave function
  bin(1:n_bin) = 0
  step_count   = 0
  et_avg       = 0.0
  pot_avg      = 0.0

  WRITE(*,'(''           step       nwalkers             e0            <v>'')')

  DO step = 1, equilibration_steps+production_steps  ! Loop over steps

     ! Move each walker forward by ds in imaginary time
     DO i = 1, n
        x(i) = x(i) + random_normal( 0.0, SQRT(ds) ) ! Brownian dynamics step
     END DO

     ! Calculate potential energy of each walker
     v(1:n) = 0.5 * x(1:n)**2

     ! Calculate the replica number for each walker
     DO i = 1, n
        de = ( v(i) - et ) *  ds           ! scaled energy difference with the trial ground state
        k  = EXP( -de )                    ! new number of replicas
        CALL RANDOM_NUMBER( zeta )
        replica(i) = FLOOR ( k+zeta )      ! integerized number of replicas
        replica(i) = MIN ( replica(i), 3 ) ! guard against too rapid growth
     END DO

     ! Remove walkers with replica = 0
     alive(1:n)       = replica(1:n) > 0
     nlive            = COUNT( alive(1:n) )
     x(1:nlive)       = PACK ( x(1:n),       alive(1:n) )
     v(1:nlive)       = PACK ( v(1:n),       alive(1:n) )
     replica(1:nlive) = PACK ( replica(1:n), alive(1:n) )

     ! Replicate walkers with replica > 1
     n  = nlive
     DO i = 1, nlive
        n_add = replica(i) - 1
        IF ( n_add > 0 ) THEN
           IF ( n + n_add > n_max ) STOP 'n exceeds n_max in run'
           x(n+1:n+n_add) = x(i)
           v(n+1:n+n_add) = v(i)
           n = n + n_add
        END IF
     END DO

     pot = SUM(v(1:n))/REAL(n) ! potential averaged over walkers

     ! Stop if n is too large or too small
     IF ( n < 3 ) THEN
        STOP 'n is too small, increase the value of et'
     ELSE IF (n > (n_max-3) ) THEN
        STOP 'n is too large, decrease the value of et'
     END IF

     ! Periodically write the current number of walkers, trial energy, average potential

     IF( MOD(step,output_interval) == 0 ) THEN
        WRITE(*,'(i15,i15,2f15.5)') step, n, et, pot
     ENDIF

     ! Calculate the distribution of walkers on the line and energy averages

     IF ( step > equilibration_steps ) THEN

        DO i = 1, n
           ibin = INT( (x(i)+x_max)/dx ) + 1
           IF ( ibin >= 1 .AND. ibin <= n_bin ) THEN
              bin(ibin) = bin(ibin) +  1
           END IF
        END DO

        et_avg     = et_avg     + et
        pot_avg    = pot_avg    + pot
        step_count = step_count + 1

     END IF

     ! Reset trial energy following Kozstin et al (1996)
     et = pot + (1.0-REAL(n)/REAL(n_target))/ds

  ENDDO ! End loop over steps

  ! Normalize the wavefunction so that the integral over psi**2 = 1.0
  psi     = REAL ( bin )         ! un-normalized wave function
  factor  = SUM ( psi**2 ) * dx  ! integral, assuming that -x_max .. +x_max catches everything
  psi     = psi / SQRT( factor ) ! normalizing factor
  et_avg  = et_avg / REAL(step_count)
  pot_avg = pot_avg / REAL(step_count)

  WRITE(*,'(a,t40,f15.5)') 'Average trial energy = ',     et_avg
  WRITE(*,'(a,t40,f15.5)') 'Average potential energy = ', pot_avg

  ! Print ground state wave function
  WRITE(*, '(''              x         psi(x)   psi_exact(x)'')')

  DO i = 1, n_bin
     x_bin  = -x_max + ( REAL(i - 1) + 0.5 ) * dx
     pexact = const * EXP(- 0.5 * x_bin * x_bin)
     WRITE(*,'(3f15.5)') x_bin, psi(i), pexact
  END DO

  DEALLOCATE ( x, v, replica, alive, bin, psi )

END PROGRAM qmc_walk_sho

