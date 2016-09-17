! corfun.f90
! Time correlation function, directly and by FFT
PROGRAM corfun
  !
  ! TODO MPA provide code

  ! Define underlying process by generalized Langevin equation
  ! with memory function expressed as sum of exponentials
  ! see G Ciccotti and JP Ryckaert Mol Phys 40 141 (1980)
  ! and AD Baczewski and SD Bond J Chem Phys 139 044107 (2013)

  integer                         :: nk    ! number of memory function terms
  real, dimension(:), allocatable :: c     ! memory function coefficients (nk)
  real, dimension(:), allocatable :: tau   ! memory function decay times (nk)
  real, dimension(:), allocatable :: theta ! auxiliary coefficient (nk)
  real, dimension(:), allocatable :: alpha ! auxiliary coefficient (nk)
  real, dimension(:), allocatable :: zeta  ! random numbers
  real, dimension(:), allocatable :: s     ! GLE auxiliary variables
  real                            :: delta ! time step
  real                            :: vt    ! velocity at time t
  real, dimension(:), allocatable :: v     ! stored velocities (0:nstep)
  REAL, DIMENSION(:), ALLOCATABLE :: vacf  ! velocity correlation function (0:nt)

  INTEGER :: nt              ! number of timesteps to correlate
  INTEGER :: n0              ! number of time origins to store
  INTEGER :: origin_interval ! interval for time origins
  INTEGER :: dt              ! time difference
  INTEGER :: t               ! time (equivalent to step number in file)

  NAMELIST /nml/ nt, origin_interval, nstep

  ! Example default values
  nt              = 1000
  origin_interval = 10 

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in corfun'
  END IF

  n0 = nt / origin_interval + 1

  nstep = 100000
  ALLOCATE ( v(0:nstep), v0(n0), t0(n0) )
  ALLOCATE ( vacf(0:nt), norm(0:nt) )

  ! Decide on model and allocate arrays
  nk = 1
  allocate ( c(nk), tau(nk), theta(nk), alpha(nk), zeta, s(nk) )
  ! Values used by Baczewski and Bond are c=tau=1 (underdamped) 0.5 (critically damped) and 0.25 (overdamped)
  c(1) = 1.0
  tau(1) = 1.0
  
  delta = 0.01
  temperature = 1.0
  stddev = sqrt(2.0*temperature)
  
  theta = exp(-delta/tau)
  alpha = sqrt((1.0-theta**2)/(2.0*tau))

  ! Initial values
  vt = 0.0
  s = 0.0
  v(0) = vt

  do t = 1, nstep
     vt = vt + 0.5*delta*sum(s)
     call random_normals ( 0.0, stddev, zeta )
     s = theta*s - (1.0-theta)*c*vt + alpha*sqrt(c)*zeta
     vt = vt + 0.5*delta*sum(s)
     v(t) = vt ! store velocities
  end do

  vacf(:) = 0.0
  norm(:) = 0.0
  t       = 0 ! time of each velocity
  mk      = 0 ! storage location of time origin
  full    = .FALSE.

  DO t = 0, nstep ! main loop correlating data

     IF ( MOD(t,origin_interval) == 0 ) THEN
        mk = mk + 1
        IF ( mk > n0 ) THEN
           full = .TRUE.
           mk = mk - n0 ! overwrite older values
        END IF
        t0(mk) = t    ! store time origin
        v0(mk) = v(t) ! store velocity origins
     END IF

     IF ( full ) THEN
        nk = n0 ! correlate with all stored time origins
     ELSE
        nk = mk ! correlate with those stored so far
     END IF

     DO k = 1, nk
        dt = t - t0(k)
        IF ( dt >= 0 .AND. dt <= nt ) THEN
           vacf(dt) = vacf(dt) + v(t) * v0(k) ! increment correlation function
           norm(dt) = norm(dt) + 1.0          ! increment normalizing factor
        END IF
     END DO

  END DO  ! End main loop correlating data
  vacf(:) = vacf(:) / norm(:) ! normalise correlation function

  DO t = 0, nt
     WRITE ( unit=output_unit, fmt='(i10,f15.8)' ) t, vacf(t)
  END DO

contains

  function vacf_exact ( t ) result ( vacf )
    real, intent(in) :: t
    real             :: vacf

    ! in general the exact correlation function may be obtained from the inverse Laplace transform
    ! vacf(s) = temperature / ( s + memory(s) ) where
    ! memory(s) = sum_k c_k / (s+tau_k) from the Laplace transform of a decaying exponential
    ! and this is expressed as a 
  end function vacf_exact
  
END PROGRAM corfun

