! diffusion.f90
! Calculates vacf and msd
PROGRAM diffusion
  USE iso_fortran_env, ONLY : iostat_end
  IMPLICIT NONE

  INTEGER :: n ! number of atoms
  INTEGER :: nt ! number of timesteps to correlate
  INTEGER :: n0 ! number of time origins to store
  INTEGER :: origin_interval ! interval for time origins
  INTEGER :: dt ! time difference
  INTEGER :: t  ! time (equivalent to record number in file)

  REAL,    DIMENSION(:,:),   ALLOCATABLE :: r, r1 ! positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: v     ! velocities (3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: r0    ! stored position origins (3,n,n0)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: v0    ! stored velocity origins (3,n,n0)
  INTEGER, DIMENSION(:),     ALLOCATABLE :: t0    ! times of origins
  REAL,    DIMENSION(:),     ALLOCATABLE :: vacf  ! velocity correlation function (0:nt)
  REAL,    DIMENSION(:),     ALLOCATABLE :: msd   ! mean-squared displacement (0:nt)
  REAL,    DIMENSION(:),     ALLOCATABLE :: norm  ! normalization factors (0:nt)

  INTEGER :: trj_unit, iostat
  LOGICAL :: full
  INTEGER :: k, mk, nk
  REAL    :: box

  NAMELIST /diffusion_parameters/ nt, origin_interval

  ! Sample defaults
  nt = 1000
  origin_interval = 10 
  READ(*,nml=diffusion_parameters)
  n0 = nt / origin_interval + 1

  OPEN(newunit=trj_unit,file='diffusion.trjinp',status='old', &
       & form = 'unformatted', access = 'sequential',action='read',iostat=iostat)
  IF ( iostat /= 0 ) STOP 'Error opening trajectory file'

  READ(trj_unit) n, box

  ALLOCATE ( r(3,n), r1(3,n), v(3,n), r0(3,n,n0), v0(3,n,n0), t0(n0) )
  ALLOCATE ( msd(0:nt), vacf(0:nt), norm(0:nt) )

  msd(:)  = 0.0
  vacf(:) = 0.0
  norm(:) = 0.0
  t       = 0 ! time of each velocity
  mk      = 0 ! storage location of time origin
  full    = .FALSE.

  DO ! Single sweep through data until end

     READ ( trj_unit, iostat=iostat ) r, v
     IF ( iostat == iostat_end ) EXIT
     IF ( iostat /= 0          ) STOP 'Error reading trajectory file'

     IF ( t > 0 ) CALL unfold ( r1, r )

     IF ( MOD(t,origin_interval) == 0 ) THEN
        mk = mk + 1
        IF ( mk > n0 ) THEN
           full = .TRUE.
           mk = mk - n0 ! overwrite older values
        END IF
        t0(mk)     = t      ! store time origin
        r0(:,:,mk) = r(:,:) ! store position origins
        v0(:,:,mk) = v(:,:) ! store velocity origins
     END IF

     IF ( full ) THEN
        nk = n0 ! correlate with all stored time origins
     ELSE
        nk = mk ! correlate with those stored so far
     END IF

     DO k = 1, nk
        dt = t - t0(k)
        IF ( dt >= 0 .AND. dt <= nt ) THEN
           msd(dt)  = msd(dt)  + SUM ( ( r(:,:) - r0(:,:,k) ) ** 2 ) ! increment msd
           vacf(dt) = vacf(dt) + SUM (   v(:,:) * v0(:,:,k)        ) ! increment correlation function
           norm(dt) = norm(dt) + 1.0                                 ! increment normalizing factor
        END IF
     END DO

     r1 = r     ! ready for next step
     t  = t + 1 ! number of next step

  END DO  ! End main loop, reading and correlating data

  CLOSE ( unit = trj_unit )

  msd(:)  = msd(:)  / norm(:) / REAL ( n ) ! normalise mean-squared displacement
  vacf(:) = vacf(:) / norm(:) / REAL ( n ) ! normalise correlation function

  DO t = 0, nt
     WRITE ( *, '(i10,2f15.8)') t, vacf(t), msd(t)
  END DO

CONTAINS

  SUBROUTINE unfold ( r1, r )
    IMPLICIT NONE

    REAL, DIMENSION(:,:), INTENT(in)    :: r1 ! previous configuration
    REAL, DIMENSION(:,:), INTENT(inout) :: r  ! current configuration

    ! Removes effects of periodic boundaries on particle trajectories
    ! Assumes that no particle actually moves more than half a box length
    ! between successive configurations

    r = r - r1                      ! convert r to displacements relative to r1
    r = r - ANINT ( r / box ) * box ! apply periodic boundaries to displacements
    r = r + r1                      ! convert r back to absolute coordinates

  END SUBROUTINE unfold

END PROGRAM diffusion
