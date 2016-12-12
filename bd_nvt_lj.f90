! bd_nvt_lj.f90
! Brownian dynamics, NVT ensemble
PROGRAM bd_nvt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE maths_module,     ONLY : random_normals
  USE md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, r, v, f, n, potential_type

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using BAOAB algorithm of BJ Leimkuhler and C Matthews
  ! Appl. Math. Res. eXpress 2013, 34–56 (2013); J. Chem. Phys. 138, 174102 (2013)
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1
  ! We assume mass m=1 throughout

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: dt          ! Time step
  REAL :: r_cut       ! Potential cutoff distance
  REAL :: temperature ! Temperature (specified)
  REAL :: gamma       ! Friction coefficient

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  ! Composite interaction = pot & cut & vir & lap & ovr variables
  TYPE(potential_type) :: total

  INTEGER :: blk, stp, nstep, nblock, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, r_cut, dt, gamma, temperature

  WRITE ( unit=output_unit, fmt='(a)' ) 'bd_nvt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Brownian dynamics, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass m=1 throughout'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 20000
  r_cut       = 2.5
  dt          = 0.005
  gamma       = 1.0
  temperature = 1.0

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in bd_nvt_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Friction coefficient',      gamma
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Ideal diffusion coefft',    temperature / gamma

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n) / box**3
  CALL allocate_arrays ( box, r_cut )
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial forces, potential, etc plus overlap check
  CALL force ( box, r_cut, total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in bd_nvt_lj'
  END IF
  CALL calculate ( 'Initial values' )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( variables )
  
  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL b_propagator ( dt/2.0 ) ! B kick half-step
        CALL a_propagator ( dt/2.0 ) ! A drift half-step
        CALL o_propagator ( dt )     ! O random velocities and friction step
        CALL a_propagator ( dt/2.0 ) ! A drift half-step

        CALL force ( box, r_cut, total )
        IF ( total%ovr ) THEN
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Error in bd_nvt_lj'
        END IF

        CALL b_propagator ( dt/2.0 ) ! B kick half-step

        ! Calculate and accumulate variables for this step
        CALL calculate ( )
        CALL blk_add ( variables )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                           ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ! Output run averages

  ! Final forces, potential, etc plus overlap check
  CALL force ( box, r_cut, total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in bd_nvt_lj'
  END IF
  CALL calculate ( 'Final values' )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE a_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    ! Implements the A propagator (drift)
    
    r(:,:) = r(:,:) + t * v(:,:) / box ! Positions in box=1 units
    r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries (box=1 units)

  END SUBROUTINE a_propagator

  SUBROUTINE b_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    ! Implements the B propagator (kick)
    
    v(:,:) = v(:,:) + t * f(:,:)

  END SUBROUTINE b_propagator

  SUBROUTINE o_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt)

    ! Implements the O propagator (friction and random contributions)
    
    REAL                 :: x, c
    REAL, DIMENSION(3,n) :: zeta

    REAL, PARAMETER :: c1 = 2.0, c2 = -2.0, c3 = 4.0/3.0, c4 = -2.0/3.0 ! Taylor series coefficients

    x = gamma * t
    IF ( x > 0.0001 ) THEN ! Use formula
       c = 1-EXP(-2.0*x)
    ELSE ! Use Taylor expansion for low x
       c = x * ( c1 + x * ( c2 + x * ( c3 + x * c4 ) ) ) 
    END IF
    c = SQRT ( c )

    CALL random_normals ( 0.0, SQRT(temperature), zeta ) ! Random momenta
    v = EXP(-x) * v + c * zeta

  END SUBROUTINE o_propagator

  SUBROUTINE calculate ( string )
    USE lrc_module,      ONLY : potential_lrc, pressure_lrc
    USE averages_module, ONLY : write_variables, msd
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: e_s, p_s, c_s, e_f, p_f, c_f, t_k, t_c
    REAL                :: kin, fsq, vol, rho

    ! Preliminary calculations
    kin = 0.5*SUM(v**2) ! Total kinetic energy
    fsq = SUM(f**2)     ! Total squared force
    vol = box**3        ! Volume
    rho = REAL(n) / vol ! Density

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Internal energy (cut-and-shifted ) per atom
    ! Total KE plus cut-and-shifted PE divided by N
    e_s = variable_type ( nam = 'E/N cut&shifted', val = (kin+total%pot)/REAL(n) )

    ! Pressure (cut-and-shifted)
    ! Ideal gas part plus total virial divided by V
    p_s = variable_type ( nam = 'P cut&shifted', val = rho*temperature + total%vir/vol )

    ! Internal energy (full, including LRC) per atom
    ! LRC plus total KE plus cut (but not shifted) PE, divided by N
    e_f = variable_type ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + (kin+total%cut)/REAL(n) )

    ! Pressure (full, including LRC)
    ! LRC + ideal gas part plus total virial divided by V plus LRC
    p_f = variable_type ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*temperature + total%vir/vol )

    ! Kinetic temperature
    ! Momentum is not conserved, hence 3N degrees of freedom
    t_k = variable_type ( nam = 'T kinetic', val = 2.0*kin/REAL(3*n) )

    ! Configurational temperature
    ! Total squared force divided by total Laplacian
    t_c = variable_type ( nam = 'T config', val = fsq/total%lap )

    ! Heat capacity (cut-and-shifted)
    ! Total energy divided by temperature and sqrt(N) to make result intensive
    c_s = variable_type ( nam = 'Cv/N cut&shifted', val = (kin+total%pot)/(temperature*SQRT(REAL(n))), method = msd )

    ! Heat capacity (full)
    ! Total energy divided by temperature and sqrt(N) to make result intensive; LRC does not contribute
    c_f = variable_type ( nam = 'Cv/N full', val = (kin+total%cut)/(temperature*SQRT(REAL(n))), method = msd )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ e_s, p_s, e_f, p_f, t_k, t_c, c_s, c_f ]
    
    IF ( PRESENT ( string ) ) THEN ! Output required
       WRITE ( unit=output_unit, fmt='(a)' ) string
       CALL write_variables ( variables(1:6) ) ! Don't write MSD variables
    END IF

  END SUBROUTINE calculate

END PROGRAM bd_nvt_lj
