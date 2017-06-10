! qmc_pi_sho.f90
! Quantum Monte Carlo, path-integral, harmonic oscillator
PROGRAM qmc_pi_sho

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

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               averages_module, ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               maths_module,    ONLY : metropolis

  IMPLICIT NONE

  REAL, DIMENSION(:), ALLOCATABLE  :: x ! Position of each bead (p)

  ! Most important variables
  INTEGER :: p           ! Number of beads in ring polymer
  REAL    :: temperature ! Specified temperature
  REAL    :: dx_max      ! Maximum Monte Carlo displacement
  REAL    :: pot_cl      ! Classical potential energy
  REAL    :: pot_qu      ! Quantum potential energy

  REAL    :: xi, zeta, beta
  REAL    :: pot_cl_old, pot_cl_new, pot_qu_old, pot_qu_new
  REAL    :: k_spring, delta, e_qu, m_ratio
  INTEGER :: i, ip1, im1, nstep, nblock, nequil, blk, stp
  INTEGER :: moves, ioerr

  NAMELIST /nml/ p, temperature, nstep, nblock, nequil, dx_max

  WRITE ( unit=output_unit, fmt='(a)' ) 'qmc_pi_sho'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Path Integral Monte Carlo simulation of a quantum oscillator'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Results in atomic units'

  CALL RANDOM_SEED() ! Initialize random number generator

  ! Set sensible default run parameters for testing
  p           = 8
  temperature = 0.2
  nstep       = 50000
  nblock      = 20
  nequil      = 10
  dx_max      = 1.0

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
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

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of beads, P = ',                 p
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature = ',                        temperature
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks for production = ',    nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks for equilibration = ', nequil
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block = ',          nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Max displacement = ',                   dx_max
  beta     = 1.0 / temperature
  k_spring = REAL(p) * temperature**2

  ALLOCATE ( x(p) )
  x = 0.0 ! Set up initial positions at origin

  ! Calculate initial values
  pot_cl = 0.5 * SUM ( x**2 ) / REAL (p )                ! Classical potential energy
  pot_qu = 0.5 * k_spring * SUM ( ( x-CSHIFT(x,1) )**2 ) ! Quantum potential energy

  ! Initialize arrays for averaging and write column headings
  m_ratio = 0.0
  CALL run_begin ( calc_variables() )

  DO blk = 1-nequil, nblock ! Begin loop over blocks (including equilibration)

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0
        DO i = 1, p ! Begin loop over beads

           ! Identify neighbours
           ip1 = i+1
           IF ( ip1 > p ) ip1 = 1
           im1 = i-1
           IF ( im1 < 1 ) im1 = p

           CALL RANDOM_NUMBER ( zeta ) ! One uniform random number in range (0,1)
           zeta = 2.0*zeta - 1.0       ! Now in range (-1,+1)

           xi         = x(i)
           pot_cl_old = 0.5 * xi**2 / REAL(p)
           pot_qu_old = 0.5 * k_spring * ( (xi-x(im1))**2 + (xi-x(ip1))**2 )
           xi         = xi + zeta * dx_max   ! Trial move to new position
           pot_cl_new = 0.5 * xi**2 / REAL(p)
           pot_qu_new = 0.5 * k_spring * ( (xi-x(im1))**2 + (xi-x(ip1))**2 )

           delta = ( pot_cl_new + pot_qu_new - pot_cl_old - pot_qu_old ) / temperature
           IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
              pot_cl = pot_cl + pot_cl_new - pot_cl_old ! Update classical potential energy
              pot_qu = pot_qu + pot_qu_new - pot_qu_old ! Update quantum potential energy
              x(i)   = xi                               ! Update position
              moves  = moves + 1                        ! Increment move counter
           END IF ! Reject Metropolis test

        END DO ! End loop over beads

        m_ratio = REAL(moves) / REAL(p)

        IF ( blk > 0 ) THEN ! Test for production phase
           ! Calculate and accumulate variables for this step
           CALL blk_add ( calc_variables() )
        END IF ! End test for production phase

     END DO ! End loop over steps

     IF ( blk > 0 ) CALL blk_end ( blk ) ! Output block averages

  END DO ! End loop over blocks (including equilibration)

  CALL run_end ( calc_variables() ) ! Output run averages

  e_qu = e_pi_sho ( p, beta )
  WRITE ( unit=output_unit, fmt='(a,i0.0,a,t40,f15.6)' ) 'Exact P=', p, ' energy', e_qu
  e_qu = 0.5 / TANH(0.5*beta)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'        ) 'Exact P=infinity energy', e_qu

CONTAINS

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE averages_module, ONLY : variable_type
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(4) :: variables ! The 4 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: m_r, pe_cl, pe_qu, energy
    REAL                :: kin

    ! Preliminary calculations
    kin = 0.5 * p * temperature  ! Kinetic energy for P-bead system

    ! Move ratio
    m_r = variable_type ( nam = 'Move:ratio', val = m_ratio, instant = .FALSE. )

    ! Classical potential energy
    pe_cl = variable_type ( nam = 'PE:classical', val = pot_cl )

    ! Quantum potential energy
    pe_qu = variable_type ( nam = 'PE:quantum', val = pot_qu )

    ! Energy
    energy = variable_type ( nam = 'Energy', val = kin+pot_cl-pot_qu )

    ! Collect together for averaging
    variables = [ m_r, pe_cl, pe_qu, energy ]

  END FUNCTION calc_variables

  FUNCTION e_pi_sho ( p, beta ) RESULT ( e )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: p
    REAL,    INTENT(in) :: beta
    REAL                :: e

    REAL :: t, s, alpha, q1, q2, q1p, q2p

    ! Exact formulae given by
    ! KS Schweizer, RM Stratt, D Chandler, and PG Wolynes, J Chem Phys, 75, 1347 (1981)
    ! M Takahashi and M Imada, J Phys Soc Japan, 53, 3765 (1984)

    ! For not-too-high P, we may express the results as a ratio of polynomials in alpha, 
    ! with integer coefficients, most conveniently in partial-fraction form.
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
    IMPLICIT NONE
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

