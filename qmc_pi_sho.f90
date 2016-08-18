! qmc_pi_sho.f90
! Quantum Monte Carlo, path-integral, harmonic oscillator
PROGRAM qmc_pi_sho

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE averages_module, ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE utility_module,  ONLY : metropolis

  IMPLICIT NONE

  ! Program to calculate the average total energy E at temperature T
  ! for a particle in a harmonic potential, V=(x**2)/2,
  ! by simulating the discretized path integral ring polymer of P beads

  ! In atomic units, classical oscillation freqency omega=1, hbar=1, mass=1
  ! so T is equivalent to kT/hbar*omega and E is equivalent to E/hbar*omega
  ! Results are output as averages over the production period.
  ! The value of <E> may be compared with the exact result for given P for this simple problem
  ! as well as the exact quantum result for P=infinity.

  ! For this simple illustration we only use crude single-particle Metropolis moves
  ! It is possible to devise smarter sampling schemes for the ring polymer

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  REAL, DIMENSION(:), ALLOCATABLE  :: x ! position of each bead (p)

  INTEGER :: p
  REAL    :: temperature, dx_max, xi, zeta, beta
  REAL    :: pot_cl, pot_cl_old, pot_cl_new, pot_qu, pot_qu_old, pot_qu_new
  REAL    :: energy, kin, k_spring, delta, move_ratio, e_qu
  INTEGER :: i, ip1, im1, nstep, production_blocks, equilibration_blocks, blk, stp
  INTEGER :: moves, ioerr

  NAMELIST /nml/ p, temperature, nstep, production_blocks, equilibration_blocks, dx_max

  WRITE ( unit=output_unit, fmt='(a)' ) 'qmc_pi_sho'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Path Integral Monte Carlo simulation of a quantum oscillator'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Results in atomic units'
  CALL time_stamp ( output_unit )

  ! Set up default parameters for the simulation
  p                    = 8     ! number of beads
  temperature          = 0.2   ! temperature
  nstep                = 50000 ! number of steps per block
  production_blocks    = 20    ! number of steps for production
  equilibration_blocks = 10    ! number of steps for equilibration
  dx_max               = 1.0   ! maximum Monte Carlo displacement
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)' ) 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)' ) 'End of file'
     STOP 'Error in qmc_pi_sho'
  END IF
  IF ( p < 2 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'p must be > 1 in this program ', p
     STOP 'Error in qmc_pi_sho'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of beads, P = ',                 p
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature = ',                        temperature
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks for production = ',    production_blocks
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks for equilibration = ', equilibration_blocks
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block = ',          nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Max displacement = ',                   dx_max
  beta     = 1.0 / temperature
  k_spring = REAL(p) * temperature**2

  ALLOCATE ( x(p) )
  CALL RANDOM_SEED()

  ! set up initial positions at origin
  x = 0.0

  ! Calculate initial values
  pot_cl = 0.5 * SUM ( x**2 ) / REAL (p )                ! Classical potential energy
  pot_qu = 0.5 * k_spring * SUM ( ( x-CSHIFT(x,1) )**2 ) ! Quantum potential energy
  kin    = 0.5 * p * temperature                         ! Kinetic energy
  energy = kin + pot_cl - pot_qu
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial classical potential energy', pot_cl
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial quantum potential energy',   pot_qu
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial energy',                     energy

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Pot-classical', 'Pot-quantum', 'Energy' ] )

  DO blk = 1, production_blocks + equilibration_blocks ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0
        DO i = 1, p ! Begin loop over beads

           ! Identify neighbours
           ip1 = i+1
           IF ( ip1 > p ) ip1 = 1
           im1 = i-1
           IF ( im1 < 1 ) im1 = p

           CALL RANDOM_NUMBER ( zeta ) ! one uniform random number in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

           xi         = x(i)
           pot_cl_old = 0.5 * xi**2 / REAL(p)
           pot_qu_old = 0.5 * k_spring * ( (xi-x(im1))**2 + (xi-x(ip1))**2 )
           xi         = xi + zeta * dx_max   ! trial move to new position
           pot_cl_new = 0.5 * xi**2 / REAL(p)
           pot_qu_new = 0.5 * k_spring * ( (xi-x(im1))**2 + (xi-x(ip1))**2 )

           delta = ( pot_cl_new + pot_qu_new - pot_cl_old - pot_qu_old ) / temperature
           IF ( metropolis ( delta ) ) THEN ! accept Metropolis test
              pot_cl    = pot_cl + pot_cl_new - pot_cl_old ! update classical potential energy
              pot_qu    = pot_qu + pot_qu_new - pot_qu_old ! update quantum potential energy
              x(i) = xi                                    ! update position
              moves  = moves + 1                           ! increment move counter
           END IF ! reject Metropolis test

        END DO ! End loop over beads

        ! Calculate all variables for this step
        move_ratio = REAL(moves) / REAL(p)
        energy = kin + pot_cl - pot_qu
        IF ( blk > equilibration_blocks ) CALL blk_add ( [move_ratio,pot_cl,pot_qu,energy] )

     END DO ! end loop over steps
     IF ( blk > equilibration_blocks ) CALL blk_end ( blk, output_unit )

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  e_qu = e_pi_sho ( p, beta )
  WRITE ( unit=output_unit, fmt='(a,i0.0,a,t40,f15.5)' ) 'Exact P=',p,' energy', e_qu
  e_qu = 0.5 / TANH(0.5*beta)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)'        ) 'Exact P=infinity energy', e_qu

  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final classical potential energy', pot_cl
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final quantum potential energy',   pot_qu
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final energy',                     energy
  pot_cl = 0.5 * SUM ( x**2 ) / REAL (p )                ! Classical potential energy
  pot_qu = 0.5 * k_spring * SUM ( ( x-CSHIFT(x,1) )**2 ) ! Quantum potential energy
  energy = kin + pot_cl - pot_qu
  WRITE ( unit=output_unit, fmt='(a)'           ) 'Final check'
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final classical potential energy', pot_cl
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final quantum potential energy',   pot_qu
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final energy',                     energy

CONTAINS

  FUNCTION e_pi_sho ( p, beta ) RESULT ( e )
    INTEGER, INTENT(in) :: p
    REAL,    INTENT(in) :: beta
    REAL                :: e

    REAL :: t, s, alpha, q1, q2, q1p, q2p

    ! Exact formulae given by
    ! KS Schweizer, RM Stratt, D Chandler, and PG Wolynes, J Chem Phys, 75, 1347 (1981)
    ! M Takahashi and M Imada, J Phys Soc Japan, 53, 3765 (1984)

    ! For any given, not-too-high, P it is probably best to express the results as a ratio 
    ! of polynomials in alpha, with integer coefficients, most conveniently in partial-fraction form.
    ! We give these up to P=8, and they are easy to obtain using a computer algebra package.

    ! Otherwise, we use the floating-point formula, but this might become
    ! inaccurate for certain values of the parameters

    IF ( p < 1 ) THEN
       WRITE ( unit=error_unit, fmt='(a,i15)') 'Error in value of p ', p
       STOP 'Error in e_pi_sho'
    END IF

    t = 1 / beta
    s = ( REAL(p)*t ) ** 2

    SELECT CASE ( p )

    CASE (1)
       e = t

    CASE (2)
       e = 1.0 + 1.0 / polynomial ( [1,4], s )
       e = e * t

    CASE (3)
       e = 1.0 + 2.0 / polynomial ( [1,3], s )
       e = e * t

    CASE (4)
       e = 1.0 + 1.0 / polynomial ( [1,4], s )
       e = e + 2.0 / polynomial ( [1,2], s )
       e = e * t

    CASE (5)
       e = 1.0 + polynomial ( [4,10], s ) / polynomial ( [1,5,5], s )
       e = e * t

    CASE (6)
       e = 1.0 + 1.0 / polynomial ( [1,4], s )
       e = e + 2.0 / polynomial ( [1,1], s )
       e = e + 2.0 / polynomial ( [1,3], s )
       e = e * t

    CASE (7)
       e = 1.0 + polynomial ( [6,28,28], s ) / polynomial ( [1,7,14,7], s )
       e = e * t

    CASE (8)
       e = 1.0 + 1.0 / polynomial ( [1,4], s )
       e = e + 2.0 / polynomial ( [1,2], s )
       e = e + polynomial ( [4,8], s ) / polynomial ( [1,4,2], s )
       e = e * t

    CASE default
       alpha = 0.5 * beta / REAL(p)
       q1 = SQRT(1.0+alpha**2) + alpha
       q2 = SQRT(1.0+alpha**2) - alpha
       q1p = q1 ** p
       q2p = q2 ** p
       e = (q1p+q2p) / ( (q1p-q2p) * (q1+q2) )
    END SELECT

  END FUNCTION e_pi_sho

  FUNCTION polynomial ( c, x ) RESULT ( f )
    INTEGER, DIMENSION(:), INTENT(in) :: c ! coefficients of x**0, x**1, x**2 etc
    REAL,                  INTENT(in) :: x ! argument
    REAL                              :: f ! result

    INTEGER :: i

    f = 0.0
    DO i = SIZE(c), 1, -1
       f = f * x + REAL ( c(i) )
    END DO

  END FUNCTION polynomial
  
END PROGRAM qmc_pi_sho

