! bd_nvt_lj.f90
! Brownian dynamics, NVT ensemble
PROGRAM bd_nvt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
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
  REAL :: density     ! Density
  REAL :: dt          ! Time step
  REAL :: r_cut       ! Potential cutoff distance
  REAL :: temperature ! Temperature (specified)
  REAL :: gamma       ! Friction coefficient

  ! Quantities to be averaged
  REAL :: e_s ! Internal energy (cut-and-shifted ) per atom
  REAL :: p_s ! Pressure (cut-and-shifted)
  REAL :: e_f ! Internal energy (full, including LRC) per atom
  REAL :: p_f ! Pressure (full, including LRC)
  REAL :: t_k ! Kinetic temperature
  REAL :: t_c ! Configurational temperature

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
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 5000
  r_cut       = 2.5
  dt          = 0.002
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
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Friction coefficient',      gamma
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Ideal diffusion coefft',    temperature / gamma

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density
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
  CALL run_begin ( [ CHARACTER(len=15) :: 'E (cut&shift)', 'P (cut&shift)', &
       &            'E (full)', 'P (full)', 'T (kin)', 'T (con)' ] )
  
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

        CALL calculate ( )                       ! Calculate instantaneous values
        CALL blk_add ( [e_s,p_s,e_f,p_f,t_k,t_c] ) ! Accumulate block averages

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )                              ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit ) ! Output run averages

  CALL force ( box, r_cut, total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in bd_nvt_lj'
  END IF
  CALL calculate ( 'Final values' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE a_propagator ( t ) ! A propagator (drift)
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    r(:,:) = r(:,:) + t * v(:,:) / box ! Positions in box=1 units
    r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries (box=1 units)

  END SUBROUTINE a_propagator

  SUBROUTINE b_propagator ( t ) ! B propagator (kick)
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    v(:,:) = v(:,:) + t * f(:,:)

  END SUBROUTINE b_propagator

  SUBROUTINE o_propagator ( t ) ! O propagator (friction and random contributions)
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt)

    REAL                 :: x, c
    REAL, DIMENSION(3,n) :: zeta

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

  SUBROUTINE calculate ( string )
    USE md_module, ONLY : potential_lrc, pressure_lrc
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates variables of interest and (optionally) writes them out

    REAL :: kin

    kin = 0.5*SUM(v**2)            ! Total KE
    t_k = 2.0 * kin / REAL ( 3*n ) ! Momentum is not conserved, hence 3N degrees of freedom

    e_s = ( total%pot + kin ) / REAL ( n )       ! total%pot is the total cut-and-shifted PE
    e_f = ( total%cut + kin ) / REAL ( n )       ! total%cut is the total cut (but not shifted) PE
    e_f = e_f + potential_lrc ( density, r_cut ) ! Add LRC

    p_s = density * temperature + total%vir / box**3 ! total%vir is the total virial
    p_f = p_s + pressure_lrc ( density, r_cut )      ! Add LRC

    t_c = SUM(f**2) / total%lap ! total%lap is the total Laplacian

    IF ( PRESENT ( string ) ) THEN ! output required
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E (cut&shift)', e_s
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (cut&shift)', p_s
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E (full)',      e_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (full)',      p_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (kin)',       t_k
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (con)',       t_c
    END IF

  END SUBROUTINE calculate

END PROGRAM bd_nvt_lj
