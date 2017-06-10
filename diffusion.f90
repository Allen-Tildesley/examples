! diffusion.f90
! Calculates vacf and msd
PROGRAM diffusion

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

  ! Reads a trajectory from a sequence of configuration files
  ! Calculates velocity autocorrelation function, mean square displacement,
  ! and cross-correlation between initial velocity and displacement
  ! Results are written to a file 'diffusion.out' with diagnostics to standard output

  ! For illustration and simplicity, we adopt a scheme of formatted files of the same kind
  ! as those that are saved at the end of each block of our MD simulation examples
  ! We assume that the initial configuration of a run has been copied to cnf.000
  ! and subsequent configurations are called cnf.001 cnf.002 etc., up to (at most) cnf.999
  ! Obviously, in a practical application, a binary trajectory file would fulfil this role.

  ! Cubic periodic boundary conditions are assumed
  ! r and box are assumed to be in the same units (e.g. LJ sigma)
  ! box is assumed to be unchanged throughout

  ! Note that we never apply periodic boundary conditions to the atomic positions
  ! We unfold the trajectory by applying PBCs to the displacements between successive configurations
  ! This assumes that the atoms never move more than box/2 during that interval

  ! Values of basic parameters are read from standard input using a namelist nml
  ! Although a default value of delta=0.05 is supplied, it is really only a place-holder
  ! for the correct user-supplied value (time interval between configurations)

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms

  IMPLICIT NONE

  INTEGER :: n               ! Number of atoms
  INTEGER :: nt              ! Number of timesteps to correlate
  INTEGER :: n0              ! Number of time origins to store
  INTEGER :: origin_interval ! Interval for time origins
  INTEGER :: dt              ! Time difference (in timesteps)
  INTEGER :: t               ! Time (in timesteps, equivalent to number of file)
  REAL    :: delta           ! Time interval (simulation units: only used in output file)

  REAL,    DIMENSION(:,:),   ALLOCATABLE :: r     ! Positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: r_old ! Previous positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: v     ! Velocities (3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: r0    ! Stored position origins (3,n,n0)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: v0    ! Stored velocity origins (3,n,n0)
  INTEGER, DIMENSION(:),     ALLOCATABLE :: t0    ! Times of origins
  REAL,    DIMENSION(:),     ALLOCATABLE :: vacf  ! Velocity correlation function (0:nt)
  REAL,    DIMENSION(:),     ALLOCATABLE :: rvcf  ! Displacement-velocity cross correlation function (0:nt)
  REAL,    DIMENSION(:),     ALLOCATABLE :: msd   ! Mean-squared displacement (0:nt)
  REAL,    DIMENSION(:),     ALLOCATABLE :: norm  ! Normalization factors (0:nt)

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3)            :: sav_tag
  INTEGER                     :: k, mk, nk, unit, ioerr
  LOGICAL                     :: full, exists
  REAL                        :: box

  NAMELIST /nml/ nt, origin_interval, delta

  ! Example default values
  nt              = 500  ! Max correlation time (as a multiple of interval between configurations)
  origin_interval = 10   ! We only take time origins at these intervals, for efficiency
  delta           = 0.05 ! This should be set to the actual time interval between configurations

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in diffusion'
  END IF

  n0 = nt / origin_interval + 1 ! Enough origins to span max correlation time

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Max correlation time nt = ',       nt
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Origin interval = ',               origin_interval
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of time origins n0 = ',     n0
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time interval between configs = ', delta

  sav_tag = '000' ! Use initial configuration to get N
  INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! Check the file exists
  IF ( .NOT. exists ) THEN
     WRITE ( unit=error_unit, fmt='(a,a)') 'File does not exist: ', cnf_prefix//sav_tag
     STOP 'Error in diffusion'
  END IF
  CALL read_cnf_atoms ( cnf_prefix//sav_tag, n, box ) ! Read n and box

  ALLOCATE ( r(3,n), r_old(3,n), v(3,n), r0(3,n,n0), v0(3,n,n0), t0(n0) )
  ALLOCATE ( msd(0:nt), rvcf(0:nt), vacf(0:nt), norm(0:nt) )

  msd(:)  = 0.0
  rvcf(:) = 0.0
  vacf(:) = 0.0
  norm(:) = 0.0
  t       = 0 ! Time of each velocity
  mk      = 0 ! Storage location of time origin
  full    = .FALSE.

  DO ! Loop reading and correlating data

     IF ( t > 999 ) EXIT ! Our naming scheme only goes up to cnf.999

     WRITE ( sav_tag, fmt='(i3.3)' ) t                      ! Number of this configuration
     INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! Check the file exists
     IF ( .NOT. exists ) EXIT                               ! We have come to the end of the data

     CALL read_cnf_atoms ( cnf_prefix//sav_tag, n, box, r, v ) ! Read configuration
     WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Processing ', t

     IF ( t > 0 ) CALL unfold ( r_old, r )

     IF ( MOD ( t, origin_interval ) == 0 ) THEN ! Test to store as time origin

        mk = mk + 1 ! Increment origin counter

        IF ( mk > n0 ) THEN ! Test for over-running origin arrays
           full = .TRUE.
           mk = mk - n0 ! Wrap-around to overwrite older origins
        END IF ! End test for over-running origin arrays

        t0(mk)     = t      ! Store time origin
        r0(:,:,mk) = r(:,:) ! Store position origins
        v0(:,:,mk) = v(:,:) ! Store velocity origins

     END IF ! End test to store as time origin

     IF ( full ) THEN
        nk = n0 ! Correlate with all stored time origins
     ELSE
        nk = mk ! Correlate with those stored so far
     END IF

     DO k = 1, nk ! Loop over time origins

        dt = t - t0(k) ! Time difference between current configuration and time origin

        IF ( dt < 0 ) THEN ! should never happen
           WRITE ( unit=error_unit, fmt='(a,i5)') 'dt error ', dt
           STOP 'Error in diffusion'
        END IF

        IF ( dt <= nt ) THEN ! Check that dt is within desired range

           ! Sum over all N atoms and over 3 Cartesian components
           msd(dt)  = msd(dt)  + SUM ( ( r(:,:) - r0(:,:,k) ) ** 2        ) ! Increment msd
           rvcf(dt) = rvcf(dt) + SUM ( ( r(:,:) - r0(:,:,k) ) * v0(:,:,k) ) ! Increment cross correlation function
           vacf(dt) = vacf(dt) + SUM (   v(:,:) * v0(:,:,k)               ) ! Increment autocorrelation function

           norm(dt) = norm(dt) + 1.0 ! Increment normalizing factor

        END IF ! End check that dt is within desired range

     END DO ! End loop over time origins

     r_old = r     ! Ready to unfold next step
     t     = t + 1 ! Number of next step

  END DO ! End loop reading and correlating data

  ! Normalize by N as well as time-origin normalizing factors
  msd(:)  = msd(:)  / norm(:) / REAL ( n ) ! 3D mean-squared displacement
  rvcf(:) = rvcf(:) / norm(:) / REAL ( n ) ! 3D cross-correlation function
  vacf(:) = vacf(:) / norm(:) / REAL ( n ) ! 3D autocorrelation function

  WRITE ( unit=output_unit, fmt='(a)' ) 'Output to diffusion.out'

  OPEN ( newunit=unit, file='diffusion.out', status='replace', iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error opening file', ioerr
     STOP 'Error in diffusion'
  END IF

  DO t = 0, nt
     WRITE ( unit=unit, fmt='(f15.6,3f15.8)' ) t*delta, vacf(t), rvcf(t), msd(t) ! Include delta in time
  END DO

  CLOSE(unit=unit)

CONTAINS

  SUBROUTINE unfold ( r_old, r )
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(in)    :: r_old ! Previous configuration
    REAL, DIMENSION(:,:), INTENT(inout) :: r     ! Current configuration

    ! Removes effects of periodic boundaries on particle trajectories
    ! Assumes that no particle actually moves more than half a box length
    ! between successive configurations

    r = r - r_old                   ! Convert r to displacements relative to r_old
    r = r - ANINT ( r / box ) * box ! Apply periodic boundaries to displacements
    r = r + r_old                   ! Convert r back to absolute coordinates

  END SUBROUTINE unfold

END PROGRAM diffusion
