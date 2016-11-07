! diffusion_test.f90
! Generate test data for diffusion.f90
PROGRAM diffusion_test

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : write_cnf_atoms
  USE maths_module,     ONLY : random_normals

  IMPLICIT NONE

  ! Generates random configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts Brownian dynamics without atomic interactions
  ! Resulting vacf should be exponential with supplied decay rate gamma
  ! Long-time mean square displacement should be diffusive with D = kT/gamma

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Input configuration, output configuration, all calculations, and results 
  ! are given in simulation units defined by the model
  ! We assume mass m=1 throughout

  ! Most important variables
  INTEGER :: n           ! Number of atoms
  REAL    :: box         ! Box length
  REAL    :: dt          ! Time step
  REAL    :: temperature ! Temperature (specified)
  REAL    :: gamma       ! Friction coefficient

  REAL, DIMENSION(:,:), ALLOCATABLE :: r    ! Positions (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: v    ! Velocities (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: zeta ! Random numbers (3,n)

  INTEGER :: blk, stp, nstep, nblock, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ n, nblock, nstep, dt, gamma, temperature, box

  WRITE ( unit=output_unit, fmt='(a)' ) 'diffusion_test'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Brownian dynamics without interactions, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass m=1 throughout'

  ! Set sensible default run parameters for testing
  n           = 100
  nblock      = 999
  nstep       = 25
  dt          = 0.002
  gamma       = 1.0
  temperature = 1.0
  box         = 1.0

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in diffusion_test'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of atoms',           n
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Friction coefficient',      gamma
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Ideal diffusion coefft',    temperature / gamma

  ALLOCATE ( r(3,n), v(3,n), zeta(3,n) )

  ! Set random positions
  CALL RANDOM_SEED ()
  CALL RANDOM_NUMBER ( r )
  r = r - 0.5 ! Now in range (-1/2,1/2)
  r = r * box ! Now in range (-box/2,box/2)

  ! Set random velocities
  CALL random_normals ( 0.0, SQRT(temperature), v )

  sav_tag = '000' ! Initial configuration output file
  CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r, v ) ! Save configuration

  DO blk = 1, nblock ! Begin loop over blocks

     DO stp = 1, nstep ! Begin loop over steps

        CALL a_propagator ( dt/2.0 ) ! A drift half-step
        CALL o_propagator ( dt )     ! O random velocities and friction step
        CALL a_propagator ( dt/2.0 ) ! A drift half-step

     END DO ! End loop over steps

     r(:,:) = r(:,:) - box * ANINT ( r(:,:) / box )             ! Periodic boundaries
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk           ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r, v ) ! Save configuration

  END DO ! End loop over blocks

  DEALLOCATE ( r, v, zeta )

CONTAINS

  SUBROUTINE a_propagator ( t ) ! A propagator (drift)
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    r(:,:) = r(:,:) + t * v(:,:) ! Positions

  END SUBROUTINE a_propagator

  SUBROUTINE o_propagator ( t ) ! O propagator (friction and random contributions)
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt)

    REAL            :: x, c
    REAL, PARAMETER :: c1 = 2.0, c2 = -2.0, c3 = 4.0/3.0, c4 = -2.0/3.0

    x = gamma * t
    IF ( x > 0.0001 ) THEN
       c = 1-EXP(-2.0*x)
    ELSE
       c = x * ( c1 + x * ( c2 + x * ( c3 + x * c4 )) ) ! Taylor expansion for low x
    END IF
    c = SQRT ( c )

    CALL random_normals ( 0.0, SQRT(temperature), zeta ) ! Random momenta
    v = EXP(-x) * v + c * zeta

  END SUBROUTINE o_propagator

END PROGRAM diffusion_test
