! diffusion.f90
! Calculates vacf and msd
PROGRAM diffusion

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms

  IMPLICIT NONE

  ! Reads a trajectory from a sequence of configuration files
  ! For illustration and simplicity, we adopt a scheme of formatted files of the same kind
  ! as those that are saved at the end of each block of our MD simulation examples
  ! We assume that the initial configuration of a run has been copied to cnf.000
  ! and subsequent configurations are called cnf.001 cnf.002 etc., up to (at most) cnf.999
  ! Obviously, in a practical application, a binary trajectory file would fulfil this role.

  ! Cubic periodic boundary conditions are assumed
  ! r and box are assumed to be in the same units (e.g. LJ sigma)
  ! box is assumed to be unchanged throughout

  ! Values of basic parameters nt, origin_interval are read from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Results are written to a file 'diffusion.out' with diagnostics to standard output

  INTEGER :: n               ! number of atoms
  INTEGER :: nt              ! number of timesteps to correlate
  INTEGER :: n0              ! number of time origins to store
  INTEGER :: origin_interval ! interval for time origins
  INTEGER :: dt              ! time difference
  INTEGER :: t               ! time (equivalent to number of file)

  REAL,    DIMENSION(:,:),   ALLOCATABLE :: r, r1 ! positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: v     ! velocities (3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: r0    ! stored position origins (3,n,n0)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: v0    ! stored velocity origins (3,n,n0)
  INTEGER, DIMENSION(:),     ALLOCATABLE :: t0    ! times of origins
  REAL,    DIMENSION(:),     ALLOCATABLE :: vacf  ! velocity correlation function (0:nt)
  REAL,    DIMENSION(:),     ALLOCATABLE :: rvcf  ! displacement-velocity cross correlation function (0:nt)
  REAL,    DIMENSION(:),     ALLOCATABLE :: msd   ! mean-squared displacement (0:nt)
  REAL,    DIMENSION(:),     ALLOCATABLE :: norm  ! normalization factors (0:nt)

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3)            :: sav_tag
  INTEGER                     :: unit, ioerr
  LOGICAL                     :: full, exists
  INTEGER                     :: k, mk, nk
  REAL                        :: box

  NAMELIST /nml/ nt, origin_interval

  ! Example default values
  nt              = 500
  origin_interval = 10 

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in diffusion'
  END IF

  n0 = nt / origin_interval + 1
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'max correlation time nt = ',   nt
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'origin interval = ',           origin_interval
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'number of time origins n0 = ', n0

  WRITE ( sav_tag, fmt='(i3.3)' ) 0 ! use initial configuration to get n
  INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! check the file exists
  IF ( .NOT. exists ) THEN
     WRITE ( unit=error_unit, fmt='(a,a)') 'File does not exist: ', cnf_prefix//sav_tag
     STOP 'Error in diffusion'
  END IF
  CALL read_cnf_atoms ( cnf_prefix//sav_tag, n, box ) ! read n and box

  ALLOCATE ( r(3,n), r1(3,n), v(3,n), r0(3,n,n0), v0(3,n,n0), t0(n0) )
  ALLOCATE ( msd(0:nt), rvcf(0:nt), vacf(0:nt), norm(0:nt) )

  msd(:)  = 0.0
  rvcf(:) = 0.0
  vacf(:) = 0.0
  norm(:) = 0.0
  t       = 0 ! time of each velocity
  mk      = 0 ! storage location of time origin
  full    = .FALSE.

  DO ! Single sweep through data until end

     IF ( t >= 1000 ) EXIT ! our naming scheme only goes up to cnf.999

     WRITE ( sav_tag, fmt='(i3.3)' ) t                      ! number of configuration
     INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! check the file exists
     IF ( .NOT. exists ) EXIT                               ! we have come to the end of the data

     CALL read_cnf_atoms ( cnf_prefix//sav_tag, n, box, r, v ) ! read configuration
     WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Processing ', t

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
           msd(dt)  = msd(dt)  + SUM ( ( r(:,:) - r0(:,:,k) ) ** 2        ) ! increment msd
           rvcf(dt) = rvcf(dt) + SUM ( ( r(:,:) - r0(:,:,k) ) * v0(:,:,k) ) ! increment cross correlation function
           vacf(dt) = vacf(dt) + SUM (   v(:,:) * v0(:,:,k)               ) ! increment autocorrelation function
           norm(dt) = norm(dt) + 1.0                                        ! increment normalizing factor
        END IF
     END DO

     r1 = r     ! ready for next step
     t  = t + 1 ! number of next step

  END DO  ! End main loop, reading and correlating data

  msd(:)  = msd(:)  / norm(:) / REAL ( n ) ! normalise mean-squared displacement
  rvcf(:) = rvcf(:) / norm(:) / REAL ( n ) ! normalise cross-correlation function
  vacf(:) = vacf(:) / norm(:) / REAL ( n ) ! normalise autocorrelation function

  WRITE ( unit=output_unit, fmt='(a)' ) 'Output to diffusion.out'

  OPEN ( newunit=unit, file='diffusion.out', status='replace', iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error opening file', ioerr
     STOP 'Error in diffusion'
  END IF

  DO t = 0, nt
     WRITE ( unit=unit, fmt='(i10,3f15.8)' ) t, vacf(t), rvcf(t), msd(t)
  END DO
  CLOSE(unit=unit)

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
