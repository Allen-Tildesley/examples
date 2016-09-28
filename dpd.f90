! dpd.f90
! Dissipative particle dynamics
PROGRAM dpd

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : lowercase
  USE dpd_module,       ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, make_ij, lowe, shardlow, r, v, f, n

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts dissipative dynamics using Shardlow or Lowe-Andersen algorithm
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in and we assume mass=1 throughout
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! The range parameter, and cutoff distance, is taken as unity

  ! The model is defined in dpd_module
  ! The typical DPD model described by Groot and Warren, J Chem Phys 107, 4423 (1997)
  ! has temperature kT=1, density rho=3, noise level sigma=3, gamma=sigma**2/2kT=4.5
  ! and force strength parameter alpha=25. They recommend timestep=0.04.
  ! They also give the approximate excess pressure as 0.101*alpha*rho**2

  ! Most important variables
  REAL :: box         ! box length
  REAL :: density     ! density
  REAL :: alpha       ! force strength parameter
  REAL :: dt          ! time step
  REAL :: gamma       ! thermalization rate (inverse time)
  REAL :: pot         ! total potential energy
  REAL :: kin         ! total kinetic energy
  REAL :: vir         ! total virial
  REAL :: lap         ! Laplacian
  REAL :: temperature ! temperature (specified)
  REAL :: pres_virial ! virial pressure (to be averaged)
  REAL :: temp_kinet  ! kinetic temperature (to be averaged)
  REAL :: temp_config ! configurational temperature (to be averaged)
  REAL :: energy      ! total energy per atom (to be averaged)

  INTEGER :: blk, stp, nstep, nblock, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number
  CHARACTER(len=10)           :: method

  ! Define a procedure pointer with an interface like that of lowe
  PROCEDURE(lowe), POINTER :: thermalize => NULL()

  NAMELIST /nml/ nblock, nstep, dt, temperature, alpha

  WRITE ( unit=output_unit, fmt='(a)' ) 'dpd'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Dissipative particle dynamics, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 and cutoff=1 throughout'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  dt          = 0.02
  temperature = 1.0
  alpha       = 25.0
  gamma       = 4.5
  method      = 'lowe'

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
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Friction / thermal rate gamma', gamma
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Specified temperature',         temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Force strength alpha',          alpha
  IF ( INDEX( lowercase(method), 'shardlow' ) /= 0 ) THEN
     thermalize => shardlow
     WRITE ( unit=output_unit, fmt='(a)' ) 'Shardlow integration method'
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
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density

  CALL allocate_arrays

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL make_ij ( box ) ! construct initial list of pairs within range
  CALL force ( box, alpha, pot, vir, lap )
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Energy', 'Temp-kinet', 'Temp-config', 'Virial Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! Shardlow or Lowe-Andersen step
        CALL thermalize ( box, temperature, gamma*dt )

        ! Velocity Verlet step
        v(:,:) = v(:,:) + 0.5 * dt * f(:,:)      ! Kick half-step
        r(:,:) = r(:,:) + dt * v(:,:) / box      ! Drift step (positions in box=1 units)
        r(:,:) = r(:,:) - ANINT ( r(:,:) )       ! Periodic boundaries
        CALL make_ij ( box )                     ! Construct list of pairs within range
        CALL force ( box, alpha, pot, vir, lap ) ! Force evaluation
        v(:,:) = v(:,:) + 0.5 * dt * f(:,:)      ! Kick half-step

        CALL calculate ( )
        
        CALL blk_add ( [energy,temp_kinet,temp_config,pres_virial] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! save configuration

  END DO ! End loop over blocks      

  CALL run_end ( output_unit )

  CALL force ( box, alpha, pot, vir, lap )
  CALL calculate ( 'Final values' )
  CALL time_stamp ( output_unit )

  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Approx EOS pressure', density*temperature+0.101*alpha*density**2

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    kin         = 0.5*SUM(v**2)
    energy      = ( pot + kin ) / REAL ( n )
    temp_kinet  = 2.0 * kin / REAL ( 3*(n-1) )
    temp_config = SUM(f**2) / lap
    pres_virial = density * temperature + vir / box**3

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Total energy',    energy
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temp-kinet',      temp_kinet
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temp-config',     temp_config
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Virial pressure', pres_virial
    END IF

  END SUBROUTINE calculate
  
END PROGRAM dpd

