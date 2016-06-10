! qmc_walk_sho.f90
! Quantum Monte Carlo, random walk, simple harmonic oscillator
PROGRAM qmc_walk_sho
  !
  ! Program to calculate the ground state wavefunction for a particle in a harmonic potential, V=(x^2)/2,
  ! by solving the corresponding diffusion equation in imaginary time.
  !
  ! Run the simulation with et = 0.5, the exact groundstate energy for this potential.
  ! Then try et = 0.6 and et = 0.4 observing the behaviour of the number of walkers in each case.
  ! This type of simulation is very senstive to the intial guess for et.
  ! Finally compare the simulated ground state wave function for the exact result for this simple problem.


  USE utility_module, ONLY: random_normal
  IMPLICIT NONE

  REAL,    DIMENSION(:), ALLOCATABLE  :: x          ! position of each walker (n_max)
  INTEGER, DIMENSION(:), ALLOCATABLE  :: replica    ! number of replicas to make (n_max)
  INTEGER, DIMENSION(:), ALLOCATABLE  :: bin        ! histogram bins for wavefunction (nbin)
  REAL,    DIMENSION(:), ALLOCATABLE  :: psi        ! wavefunction (nbin)
  REAL, PARAMETER                     :: diff = 0.5 ! diffusion coefficient in atomic units

  REAL, PARAMETER :: pi = 4.0 * ATAN ( 1.0 ), const = (1.0 / pi) ** (1.0/4.0)

  INTEGER         :: n, n_max, n_target, nbin, nk
  INTEGER         :: i, j, maxstep, step, count
  INTEGER         :: nlive, ibin, bincount
  REAL            :: lx, de, ds, var, et, resid, lx2, k
  REAL            :: zeta, dx
  REAL            :: factor, pexact, xbin

  NAMELIST /params/ maxstep, lx, et, ds, dx, n_max, n_target

  WRITE(*,'(''qmc_walk_sho'')')
  WRITE(*,'(''diffusion Monte Carlo simulation of a quantum oscillator'')')
  WRITE(*,'(''Results in atomic units'')')

  !set up parameters for the simulation

  n_max = 1000 ! max number of walkers
  n_target = 100 ! target number of walkers
  maxstep = 40000                 ! number of steps
  lx      = 8.0                   ! length of the line, centred at x = 0.0
  et      = 0.5                   ! initial estimate of ground state energy in atomic units for spring contant k=1
  ds      = 0.001                 ! timestep in imaginary time
  dx      = 0.2                   ! bin interval for calculating the wavefunction
  READ(*,nml=params)
  var     = SQRT(2.0 * diff * ds) ! the variance for gaussian random number
  lx2     = lx / 2.0              ! line goes from -lx2 to lx2
  n = n_target
  IF ( n > n_max ) STOP 'n exceeds n_max at start'
  nbin     = NINT(lx /dx)

  ALLOCATE ( x(n_max), replica(n_max) )
  ALLOCATE ( psi(nbin), bin(nbin) )

  CALL RANDOM_SEED()

  ! set up an initial delta distribution of walkers at origin
  x = 0.0

  ! set up the array for calculating the wave function

  bincount    = 0
  bin = 0


  ! loop over steps

  WRITE(*,'(''    step   nwalkers        e0'')')

  DO step = 1, maxstep

     ! Calculate the replica number for each walker

     DO i = 1, n

        x(i) = x(i) + random_normal( 0.0, var ) ! move walker forward by ds in imaginary time
        de   = (0.5 * x(i) * x(i) - et) *  ds   ! calculate the scaled energy difference with the trial ground state

        IF ( de  < 0.0 ) THEN                         ! k greater than 1
           k = EXP( -de )
           nk = FLOOR( k )
           resid = k - REAL( nk )
           CALL RANDOM_NUMBER( zeta )
           IF ( resid >= zeta ) THEN
              replica(i) = MIN ( nk + 1, 3 ) ! following Kozstin et al (1996)
           ELSE
              replica(i) = MIN ( nk, 3 ) ! following Kozstin et al (1996)
           ENDIF
        ELSE                                    ! k less than 1
           k = EXP( -de )
           resid = 1.0 - k
           CALL RANDOM_NUMBER( zeta )
           IF( resid >= zeta ) THEN
              replica(i) = 0
           ELSE
              replica(i) = 1
           ENDIF
        ENDIF

     ENDDO

     ! remove walkers with replica=0
     nlive            = COUNT( replica(1:n) > 0 )
     x(1:nlive)       = PACK ( x(1:n),       replica(1:n) > 0 )
     replica(1:nlive) = PACK ( replica(1:n), replica(1:n) > 0 )

     ! expand the array x according to the replica number

     n  = nlive

     DO i = 1, nlive

        IF ( replica(i) > 1 ) THEN ! only if more replicas are needed
           DO j = 1, replica(i) - 1 
              n = n + 1
              IF ( n > n_max ) STOP 'n exceeds n_max in run'
              x(n) = x(i)
           END DO
        END IF

     END DO

     ! Stop if n is too large or too small

     IF( n < 3 ) THEN
        STOP 'n is too small, increase the value of et'
     ELSEIF(n > (n_max-3) ) THEN
        STOP 'n is too large, decrease the value of et'
     ENDIF

     ! Periodically write the current number of walkers

     IF( MOD(step,50) == 0 ) THEN
        WRITE(*,'(i7,i9,f16.6)') step, n, et
     ENDIF

     ! Periodically calculate the distribution of walkers on the line

     IF( MOD(step, 5 ) == 0 ) THEN

        DO i = 1, n

           ibin = INT( (x(i)+lx2)/dx ) + 1
           IF ( ibin >= 0 .AND. ibin <= nbin ) THEN
              bin(ibin) = bin(ibin) +  1
           END IF
        ENDDO
        bincount = bincount + n
     ENDIF

     ! reset trial energy
     et = 0.5 * SUM ( x(1:n)**2)/REAL(n) + 0.001*(1.0-REAL(n)/REAL(n_target)) ! following Kozstin et al (1996)
     
  ENDDO ! end over steps

  ! normalize the wavefunction so that the integral over psi**2 = 1.0
  psi    =  bin * lx  / REAL( bincount) / dx
  factor = SUM ( psi**2 )
  factor = SQRT( factor * dx )
  psi    = psi / factor

  ! print ground state wave function
  WRITE(*,'(/)')
  WRITE(*, '(''     x            psi(x)    psi_exact(x)'')')
  
  DO i = 1, nbin

     xbin    = -lx2 + ( REAL(i - 1) + 0.5 ) * dx
     pexact  = const * EXP(- 0.5 * xbin * xbin)
     WRITE(*,'(3f12.6)') xbin, psi(i), pexact

  END DO

  DEALLOCATE ( x, replica, bin, psi )

END PROGRAM qmc_walk_sho

