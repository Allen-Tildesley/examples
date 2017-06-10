! corfun.f90
! Time correlation function, directly and by FFT
PROGRAM corfun

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

  ! Define underlying process by generalized Langevin equation
  ! with memory function expressed as a decaying exponential function
  ! See G Ciccotti and JP Ryckaert Mol Phys 40 141 (1980)
  ! and AD Baczewski and SD Bond J Chem Phys 139 044107 (2013)

  ! We assume that the program is linked with the FFTW library (version 3)
  ! but it could easily be adapted to use a different library
  ! With older FFT routines it was necessary to transform a number of points = exact power of 2
  ! which, for us, would mean nstep = 2**j (for integer j) and fft_len=2*nstep
  ! This is no longer true, so we do not restrict nstep in this way,
  ! but the transform should be more efficient if nstep happens to take such values
  ! Advantage can be taken of the fact that the data is real, but for clarity
  ! we just use the complex FFT with imaginary parts of the data set to zero

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE, INTRINSIC :: iso_c_binding
  USE               maths_module,    ONLY : random_normal

  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

  REAL :: m     ! Memory function coefficients
  REAL :: kappa ! Memory function decay rates
  REAL :: s     ! GLE auxiliary variables
  REAL :: delta ! Time step
  REAL :: vt    ! Velocity at time t

  REAL,    DIMENSION(:), ALLOCATABLE :: v     ! Stored velocities (nstep)
  REAL,    DIMENSION(:), ALLOCATABLE :: v0    ! Stored velocity origins (n0)
  INTEGER, DIMENSION(:), ALLOCATABLE :: t0    ! Times of origins (n0)
  REAL,    DIMENSION(:), ALLOCATABLE :: c     ! Velocity correlation function (direct method) (0:nt)
  REAL,    DIMENSION(:), ALLOCATABLE :: c_fft ! Velocity correlation function (FFT method) (0:nt)
  REAL,    DIMENSION(:), ALLOCATABLE :: n     ! Normalizing function (0:nt)

  INTEGER :: nt              ! Number of timesteps to correlate
  INTEGER :: nstep           ! Number of timesteps in run
  INTEGER :: nequil          ! Number of equilibration timesteps
  INTEGER :: n0              ! Number of time origins to store
  INTEGER :: origin_interval ! Interval for time origins
  INTEGER :: dt              ! Time difference
  INTEGER :: t               ! Time (equivalent to step number in file)

  LOGICAL :: full
  INTEGER :: k, mk, nk
  INTEGER :: unit, ioerr
  REAL    :: temperature, stddev, x, e, b, d
  REAL    :: cpu_1, cpu_2, cpu_3, cpu_4

  INTEGER(C_INT)                                       :: fft_len  ! the number of points for FFT
  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: fft_inp  ! data to be transformed (0:fft_len-1)
  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: fft_out  ! data to be transformed (0:fft_len-1)
  TYPE(C_PTR)                                          :: fft_plan ! plan needed for FFTW

  ! Taylor series coefficients
  REAL, PARAMETER :: b1 = 2.0, b2 = -2.0,     b3 = 4.0/3.0, b4 = -2.0/3.0  ! b = 1.0 - EXP(-2.0*x)
  REAL, PARAMETER :: d1 = 1.0, d2 = -1.0/2.0, d3 = 1.0/6.0, d4 = -1.0/24.0 ! d = 1.0 - EXP(-x)

  NAMELIST /nml/ nt, origin_interval, nstep, nequil, delta, temperature

  WRITE( unit=output_unit, fmt='(a)' ) 'corfun'
  WRITE( unit=output_unit, fmt='(a)' ) 'Illustrates methods for calculating time correlation functions'
  WRITE( unit=output_unit, fmt='(a)' ) 'using synthetic data from a generalized Langevin equation'

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Example default values
  ! Agreement (to numerical precision) of direct and FFT methods is expected if origin_interval=1
  nt              = 1000  ! Max time for correlation function
  origin_interval = 1     ! This could reasonably be increased to 10 or 20 to improve efficiency
  nstep           = 2**20 ! Number of steps, about a million for example
  nequil          = 10000 ! Number of equilibration timesteps
  delta           = 0.01  ! Timestep for simulation
  temperature     = 1.0   ! Temperature for simulation

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in corfun'
  END IF

  n0 = nt / origin_interval + 1
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps in run = ',    nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Equilibration steps = ',       nequil
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Max correlation time nt = ',   nt
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Origin interval = ',           origin_interval
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of time origins n0 = ', n0
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step delta = ',           delta
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature = ',               temperature

  ALLOCATE ( v(nstep), v0(n0), t0(n0) )
  ALLOCATE ( c(0:nt), c_fft(0:nt), n(0:nt) )
  fft_len = 2*nstep ! Actual length of FFT data
  ALLOCATE ( fft_inp(0:fft_len-1), fft_out(0:fft_len-1) )

  ! The memory function model is defined here
  ! Values used by Baczewski and Bond in their example are
  ! (m,kappa)
  ! (1.0,1.0)  (underdamped)
  ! (0.5,2.0)  (critically damped)
  ! (0.25,4.0) (overdamped)
  m     = 1.0
  kappa = 1.0
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.6))' ) 'm = ',     m
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.6))' ) 'kappa = ', kappa

  ! Coefficients used in algorithm
  x = delta*kappa
  e = EXP(-x) ! theta in B&B paper
  IF ( x > 0.0001 ) THEN
     b = 1.0 - EXP(-2.0*x)
     d = 1.0 - EXP(-x)
  ELSE  ! Taylor expansions for low x
     b = x * ( b1 + x * ( b2 + x * ( b3 + x * b4 ) ) )
     d = x * ( d1 + x * ( d2 + x * ( d3 + x * d4 ) ) )
  END IF
  b      = SQRT ( b )
  b      = b * SQRT ( kappa/2.0 ) ! alpha in B&B paper
  stddev = SQRT(2.0*temperature)

  ! Data generation
  CALL CPU_TIME ( cpu_1 )

  ! Initial values
  vt = 0.0
  s  = 0.0

  DO t = -nequil, nstep ! Loop over steps including an equilibration period

     ! Velocity Verlet type algorithm for vt and auxiliary variable s
     vt = vt + 0.5 * delta * s
     s  = e * s - d * m * vt + b * SQRT(m) * random_normal ( 0.0, stddev )
     vt = vt + 0.5*delta*s

     IF ( t > 0 ) v(t) = vt ! Store velocities, after equilibration

  END DO ! End loop over steps including an equilibration period

  CALL CPU_TIME ( cpu_2 )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'CPU time to generate data = ', cpu_2-cpu_1

  ! Data analysis (direct method)

  c(:) = 0.0
  n(:) = 0.0
  mk   = 0 ! Storage location of time origin
  full = .FALSE.

  DO t = 1, nstep ! Main loop correlating data

     IF ( MOD(t-1,origin_interval) == 0 ) THEN
        mk = mk + 1
        IF ( mk > n0 ) THEN
           full = .TRUE.
           mk = mk - n0 ! Overwrite older values
        END IF
        t0(mk) = t    ! Store time origins
        v0(mk) = v(t) ! Store velocity at time origins
     END IF

     IF ( full ) THEN
        nk = n0 ! Correlate with all stored time origins
     ELSE
        nk = mk ! Correlate with those stored so far
     END IF

     DO k = 1, nk ! Loop over time origins

        dt = t - t0(k)

        IF ( dt >= 0 .AND. dt <= nt ) THEN ! Check that dt is in range
           c(dt) = c(dt) + v(t) * v0(k) ! Increment correlation function
           n(dt) = n(dt) + 1.0          ! Increment normalizing factor
        END IF ! End check that dt is in range

     END DO ! End loop over time origins

  END DO ! End main loop correlating data

  IF ( ANY ( n(:) < 0.5 ) ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Normalization array error'
     STOP 'Error in corfun'
  END IF

  c(:) = c(:) / n(:) ! Normalise by number of increments

  CALL CPU_TIME ( cpu_3 )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'CPU time for direct method = ', cpu_3-cpu_2

  ! Data analysis (FFT method)

  ! Prepare data for FFT
  fft_inp            = CMPLX(0.0) ! Fill input array with zeros
  fft_inp(0:nstep-1) = CMPLX(v)   ! Put data into first part (real only)

  ! Forward FFT
  fft_plan = fftw_plan_dft_1d ( fft_len, fft_inp, fft_out, FFTW_FORWARD, FFTW_ESTIMATE) ! Set up plan
  CALL fftw_execute_dft ( fft_plan, fft_inp, fft_out )                                  ! Execute FFT
  CALL fftw_destroy_plan ( fft_plan )                                                   ! Release plan

  fft_out = fft_out * CONJG ( fft_out ) ! Square modulus

  ! Reverse FFT
  fft_plan = fftw_plan_dft_1d ( fft_len, fft_out, fft_inp, FFTW_BACKWARD, FFTW_ESTIMATE) ! Set up plan
  CALL fftw_execute_dft ( fft_plan, fft_out, fft_inp )                                   ! Execute FFT
  CALL fftw_destroy_plan ( fft_plan )                                                    ! Release plan

  fft_inp = fft_inp / REAL ( fft_len ) ! Normalization factor associated with FFT itself

  n(:) = [ ( nstep-t, t = 0, nt ) ] ! Normalization factors associated with number of time origins
  IF ( ANY ( n(:) < 0.5 ) ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Normalization array error'
     STOP 'Error in corfun'
  END IF
  c_fft   = REAL ( fft_inp(0:nt) ) / n(:) ! Apply normalization associated with number of time origins

  CALL CPU_TIME ( cpu_4 )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'CPU time for FFT method = ', cpu_4-cpu_3

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
    REAL             :: c ! Returns analytically known correlation function, for
    REAL, INTENT(in) :: t ! given time

    ! See AD Baczewski and SD Bond J Chem Phys 139 044107 (2013)
    ! In general the exact correlation function may be obtained from the inverse Laplace transform
    ! C(s) = (kT/m) * 1 / ( s + M(s) ) where both C(t) and M(t) are sums of exponentials in t
    ! M(s) = sum_p m_p*kappa_p / (s+kappa_p) p = 1..Np (note amplitude mp*kappa_p)
    ! C(s) = sum_p c_p / (s+k_p) p = 1..Np+1
    ! The coefficients c_p and k_p may be determined in terms of the m_p and kappa_p
    ! by solving an equation of order Np+1 and using partial fractions
    ! Here we just do this for the simplest case Np=1
    ! Agreement with the simulated function is only expected to within statistical error

    REAL :: del, omega, kp, km, cp, cm

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

