! dpd.f90
! Dissipative particle dynamics
PROGRAM dpd

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : lowercase
  USE dpd_module,       ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, lowe, shardlow, r, v, f, n, potential_type

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts dissipative particle dynamics using Shardlow or Lowe-Andersen algorithm
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in and we assume mass=1 throughout
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! The range parameter (cutoff distance) is taken as unity

  ! The model is defined in dpd_module
  ! The typical DPD model described by Groot and Warren, J Chem Phys 107, 4423 (1997)
  ! has temperature kT=1, density rho=3, noise level sigma=3, gamma=sigma**2/(2*kT)=4.5
  ! and force strength parameter a=25 (more generally 75*kT/rho).
  ! They recommend timestep=0.04, but we recommend a somewhat smaller value.
  ! They also give the approximate pressure as rho*kT + alpha*a*rho**2 where alpha=0.101

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: density     ! Density
  REAL :: a           ! Force strength parameter
  REAL :: dt          ! Time step
  REAL :: gamma       ! Thermalization rate (inverse time)
  REAL :: temperature ! Temperature (specified)

  ! Quantities to be averaged
  REAL :: p  ! Pressure
  REAL :: en ! Internal energy per atom
  REAL :: tk ! Kinetic temperature
  REAL :: tc ! Configurational temperature

  ! Composite interaction = pot & vir & lap variables
  TYPE(potential_type) :: total

  REAL, PARAMETER :: alpha = 0.101 ! constant in approximate equation of state

  INTEGER :: blk, stp, nstep, nblock, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number
  CHARACTER(len=10)           :: method

  ! Define a procedure pointer with an interface like that of lowe
  PROCEDURE(lowe), POINTER :: thermalize => NULL()

  NAMELIST /nml/ nblock, nstep, dt, temperature, a, gamma, method

  WRITE ( unit=output_unit, fmt='(a)' ) 'dpd'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Dissipative particle dynamics, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 and cutoff=1 throughout'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  dt          = 0.02
  temperature = 1.0
  a           = 75.0 ! actually a*rho/kT: to be multiplied by kT/rho later
  gamma       = 4.5
  method      = 'Lowe'

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in dpd'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',              nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',     nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Time step',                     dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Specified temperature',         temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Force strength a*rho/kT',       a
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Friction / thermal rate gamma', gamma

  IF ( INDEX( lowercase(method), 'shardlow' ) /= 0 ) THEN
     thermalize => shardlow
     WRITE ( unit=output_unit, fmt='(a)' ) 'Shardlow integration method'
     WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'DPD sigma parameter', SQRT ( 2.0 * gamma * temperature )
  ELSE IF ( INDEX( lowercase(method), 'lowe' ) /= 0 ) THEN
     thermalize => lowe
     WRITE ( unit=output_unit, fmt='(a)' ) 'Lowe thermalization method'
     IF ( gamma*dt > 1.0 ) THEN
        WRITE ( unit=error_unit, fmt='(a,f15.5)') 'gamma*dt too large', gamma*dt
        STOP 'Error in dpd'
     END IF
  ELSE
     WRITE ( unit=error_unit, fmt='(a,a)' ) 'Unrecognized thermalization method ', method
     STOP 'Error in dpd'
  END IF

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  a       = a * temperature / density ! Scale force strength accordingly
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density',          density
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Force strength a', a

  CALL allocate_arrays ( box )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial energy etc
  CALL force ( box, a, total )
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'E/N', 'T (kin)', 'T (con)', 'P' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! Shardlow or Lowe-Andersen step
        CALL thermalize ( box, temperature, gamma*dt )

        ! Velocity Verlet step
        CALL kick_propagator ( dt/2.0 ) ! Kick half-step
        CALL drift_propagator ( dt )    ! Drift step 
        CALL force ( box, a, total )    ! Force evaluation
        CALL kick_propagator ( dt/2.0 ) ! Kick half-step

        ! Calculate all variables for this step
        CALL calculate ( )
        CALL blk_add ( [en,tk,tc,p] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! Save configuration

  END DO ! End loop over blocks      

  CALL run_end ( output_unit )

  ! Final energy etc
  CALL force ( box, a, total )
  CALL calculate ( 'Final values' )
  CALL time_stamp ( output_unit )

  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Approx EOS pressure', density*temperature+alpha*a*density**2

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE kick_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    v(:,:) = v(:,:) + t * f(:,:)

  END SUBROUTINE kick_propagator

  SUBROUTINE drift_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt)

    r(:,:) = r(:,:) + t * v(:,:) / box ! Positions in box=1 units
    r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  END SUBROUTINE drift_propagator

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates the properties of interest from total values
    ! and optionally writes them out (e.g. at the start and end of the run)
    ! The DPD potential is short ranged, zero at, and beyond, r_cut
    ! so issues of shifted potentials and long-range corrections do not arise

    REAL :: kin

    kin = 0.5*SUM(v**2)
    tk  = 2.0 * kin / REAL ( 3*(n-1) )

    en = ( total%pot + kin ) / REAL ( n ) ! total%pot is the total PE

    p = density * temperature + total%vir / box**3 ! total%vir is the total virial

    tc = SUM(f**2) / total%lap ! total%lap is the total Laplacian

    IF ( PRESENT ( string ) ) THEN ! Output required
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N',     en
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (kin)', tk
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (con)', tc
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P',       p
    END IF

  END SUBROUTINE calculate

END PROGRAM dpd

