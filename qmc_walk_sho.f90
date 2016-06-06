! qmc_walk_sho.f90
! Quantum Monte Carlo, random walk, simple harmonic oscillator
PROGRAM qmc_walk_sho
  !
  ! Program to calculate the ground state wavefunction for a particle in a harmonic potential, V=(x^2)/2, by
  ! solving the corresponding diffusion equation in imaginary time.
  !
  ! Run the simultion with et = 0.5, the exact groundstate energy for this potential.
  ! Then try et = 0.6 and et = 0.4 observing the behaviour of the number of walkers in each case.
  ! This type of simulation is very senstive to the intial guiess for et.
  ! Finally compare the simulated ground state wave function for the exact result for this simple problem.


  USE utility_module

  IMPLICIT NONE

  INTEGER, PARAMETER                  :: nw = 100         ! initial number of walkers on the line
  INTEGER, PARAMETER                  :: nw10 = 10*nw     ! maximum number of walkers allowed
  INTEGER, PARAMETER                  :: nb = 300         ! Maximum number of bin entries for wavefunction
  REAL                                :: x(nw10)
  INTEGER                             :: replica(nw10)
  INTEGER                             :: bin(nb)
  REAL                                :: psi(nb)

  INTEGER                             :: ntemp, nl, nk
  INTEGER                             :: i, j, maxstep, step
  INTEGER                             :: nc, nzero, nupper, nlive, ibin, nbin, bincount
  REAL                                :: lx, de, ds, diff, var, et, xi, resid, lx2, k
  REAL                                :: random_normal, zeta, rep1, dx
  REAL                                :: sum, factor, const, pexact, xbin


  WRITE(*,'(''quantum_mc_diffusion'')')
  WRITE(*,'(''diffusion Monte Carlo simultion of a quantum oscillator'')')
  WRITE(*,'(''Results in atomic units'')')

  !set up parameters for the simulation

  maxstep = 30000                 ! number of steps
  lx      = 8.0                   ! length of the line, centred at x = 0.0
  diff    = 0.5                   ! diffusion coefficient in atomic units
  et      = 0.5                   ! initial estimate of ground state energy in atomic units for spring contant k=1
  ds      = 0.001                 ! timestep in imaginary time
  dx      = 0.2                   ! bin interval for calculating the wavefunction
  var     = sqrt(2.0 * diff * ds) !the variance for gaussian random number
  lx2     = lx / 2.0              ! line goes from -lx2 to lx2
  nc      = nw                    ! current number of walkers

  ! set up an intial random distribution of walkers on a line

  CALL RANDOM_SEED() ! Inialise random number generator
  DO i = 1, nc

     CALL RANDOM_NUMBER( zeta )
     x(i) = -lx2 + lx * zeta

  ENDDO

  ! set up the array for calculating the wave function

  nbin        = nint(lx /dx)
  bincount    = 0

  DO i = 1, nbin

     bin(i) = 0

  ENDDO


  ! loop over steps

  step = 0
  write(*,'(''    step   nwalkers        e0'')')

  DO WHILE (step < maxstep )

     step = step + 1

     ! Calculate the replica number for each walker

     DO i = 1, nc

        xi    = x(i) + random_normal( 0.0, var )    ! move walker forward by ds in imaginary time
        x(i)  = xi                                  ! update the position array
        de    = (0.5 * xi * xi - et) *  ds          ! calculate the scaled energy difference with the trial ground state

        IF( de  < 0.0) THEN                         ! k greater than 1
           k = exp( - de )
           nk = floor( k )
           CALL RANDOM_NUMBER( zeta )
           resid = k - real( nk )
           IF ( resid >= zeta ) THEN
              replica(i) = nk + 1
           ELSE
              replica(i) = nk
           ENDIF
        ELSE                                    ! k less than 1
           k = exp( - de )
           resid = 1.0 - k
           CALL RANDOM_NUMBER( zeta )
           IF( resid >= zeta ) THEN
              replica(i) = 0
           ELSE
              replica(i) = 1
           ENDIF
        ENDIF

     ENDDO

     !sort the first nc members of the array replica in descending; organise the elements of array x in the same order

     CALL bubble_sort( replica, x, nc )

     ! count the number of end zeros

     nzero  = 0
     DO i = 1, nc

        IF( replica(i) == 0 ) THEN
           nzero = nzero + 1
        ENDIF

     ENDDO

     ! expand the x array with new replicas

     nlive   = nc - nzero                    ! current number of live particles
     nupper  = nlive                         ! upper limit of live array
     DO i = 1, nlive

        IF( replica(i) > 1 ) THEN
           rep1 = replica(i) - 1
           DO j = 1, rep1

              x(nupper+j) = x(i)

           ENDDO
           nupper = nupper + rep1
        ENDIF
     ENDDO
     nc = nupper ! update number of walkers for next step

     ! Stop if nc is too large or too small

     IF( nc < 3 ) THEN
        STOP 'nc is to small, increase the value of et'
     ELSEIF(nc > (nw10-3) ) THEN
        STOP 'nc is too large, decrease the value of et'
     ENDIF

     ! Reset replica number to one

     DO i = 1, nc

        replica(i) = 1

     ENDDO

     ! Periodically write the current number of walkers

     IF( MOD(STEP,50) == 0 ) THEN
        write(*,'(i7,i9,f16.6)') step, nc, et
     ENDIF

     ! Periodically calculate the distribution of walkers on the line

     IF( MOD(STEP, 5 ) == 0 ) THEN

        DO i = 1, nc

           xi          = x(i)
           ibin        = int( (xi+lx2)/dx ) + 1
           bin(ibin)   = bin(ibin) +  1

        ENDDO
        bincount = bincount + nc
     ENDIF

  ENDDO ! end over steps

  ! normalize the wavefunction so that the integral over psi**2 = 1.0 

  sum = 0.0

  DO i = 1, nbin

     psi(i)  =  bin(i) * lx  / REAL( bincount) / dx
     sum = sum + psi(i) * psi(i)

  ENDDO

  factor = sqrt( sum * dx )

  ! print ground state wave function
  WRITE(*,'(/)')
  WRITE(*, '(''     x            psi(x)    psi_exact(x)'')')
  const = (1.0 / 3.14159) ** 0.25
  DO i = 1, nbin

     xbin    = -lx2 + real(i - 1) * dx + dx/2.0
     psi(i)  = psi(i) / factor
     pexact  = const * exp(- 0.5 * xbin * xbin)
     WRITE(*,'(3f12.6)') xbin, psi(i), pexact

  ENDDO

CONTAINS

  SUBROUTINE bubble_sort ( ix, fx, num )

    ! sorts an array ix in descending order and moves the corresponding values of an associated array fx

    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(INOUT)           :: ix    ! list
    REAL,    DIMENSION(:), INTENT(INOUT)           :: fx    ! variables associated with list
    INTEGER    num                                          ! size of list
    INTEGER    tempi, j, k
    REAL       tempf
    LOGICAL    sorted

    sorted  = .FALSE.
    k       = 0
    DO WHILE (.NOT. sorted)

       sorted  = .TRUE.
       k       = k + 1

       DO j = 1, num - k

          IF( ix(j) < ix(j+1) ) THEN
             tempi   = ix(j)
             tempf   = fx(j)
             ix(j)   = ix(j+1)
             fx(j)   = fx(j+1)
             ix(j+1) = tempi
             fx(j+1) = tempf
             sorted  = .FALSE.
          ENDIF

       ENDDO

    ENDDO

    RETURN
  END SUBROUTINE bubble_sort

END PROGRAM qmc_walk_sho
