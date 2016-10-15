! bd_nvt_lj.f90
! Brownian dynamics, NVT ensemble
PROGRAM bd_nvt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : random_normals
  USE md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, r, v, f, n, energy_lrc

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using BAOAB algorithm of BJ Leimkuhler and C Matthews
  ! Appl. Math. Res. eXpress 2013, 34–56 (2013); J. Chem. Phys. 138, 174102 (2013)
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in, and we assume mass=1 throughout
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  ! Most important variables
  REAL :: box         ! box length
  REAL :: density     ! density
  REAL :: dt          ! time step
  REAL :: r_cut       ! potential cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: pot_sh      ! total shifted potential energy
  REAL :: kin         ! total kinetic energy
  REAL :: vir         ! total virial
  REAL :: lap         ! total Laplacian
  REAL :: temperature ! temperature (specified)
  REAL :: pres_virial ! virial pressure (to be averaged)
  REAL :: temp_kinet  ! kinetic temperature (to be averaged)
  REAL :: temp_config ! configurational temperature (to be averaged)
  REAL :: energy      ! total energy per atom (to be averaged)
  REAL :: energy_sh   ! total shifted energy per atom (to be averaged)
  REAL :: gamma       ! friction coefficient

  INTEGER :: blk, stp, nstep, nblock, ioerr
  REAL    :: pot_lrc, vir_lrc

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, r_cut, dt, gamma, temperature

  WRITE ( unit=output_unit, fmt='(a)' ) 'bd_nvt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Brownian dynamics, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  r_cut       = 2.5
  dt          = 0.005
  gamma       = 1.0
  temperature = 0.7

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in bd_nvt_lj'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Friction coefficient',      gamma
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Ideal diffusion coefft',    temperature / gamma

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density

  CALL allocate_arrays ( box, r_cut )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL force ( box, r_cut, pot, pot_sh, vir, lap )
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Energy', 'Shifted Energy', 'Temp-kinet', 'Temp-config', 'Virial Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! BAOAB algorithm
        CALL b_propagator ( dt/2.0 )                      ! B kick half-step
        CALL a_propagator ( dt/2.0 )                      ! A drift half-step
        CALL o_propagator ( dt )                          ! O random velocities and friction step
        CALL a_propagator ( dt/2.0 )                      ! A drift half-step
        CALL force ( box, r_cut, pot, pot_sh, vir, lap )  ! Force evaluation
        CALL b_propagator ( dt/2.0 )                      ! B kick half-step

        CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )

        ! Calculate all variables for this step
        CALL calculate ( )
        CALL blk_add ( [energy,energy_sh,temp_kinet,temp_config,pres_virial] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL force ( box, r_cut, pot, pot_sh, vir, lap )
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  CALL calculate ( 'Final values' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE a_propagator ( t ) ! A propagator (drift)
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! time over which to propagate (typically dt/2)
    r(:,:) = r(:,:) + t * v(:,:) / box ! positions in box=1 units
    r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! periodic boundaries (box=1 units)
  END SUBROUTINE a_propagator

  SUBROUTINE b_propagator ( t ) ! B propagator (kick)
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! time over which to propagate (typically dt/2)
    v(:,:) = v(:,:) + t * f(:,:)
  END SUBROUTINE b_propagator

  SUBROUTINE o_propagator ( t ) ! O propagator (friction and random contributions)
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! time over which to propagate (typically dt)

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

    CALL random_normals ( 0.0, SQRT(temperature), zeta ) ! random momenta
    v = EXP(-x) * v + c * zeta

  END SUBROUTINE o_propagator

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    kin         = 0.5*SUM(v**2)
    energy      = ( pot + pot_lrc + kin ) / REAL ( n )
    energy_sh   = ( pot_sh + kin ) / REAL ( n )
    temp_kinet  = 2.0 * kin / REAL ( 3*(n-1) )
    temp_config = SUM ( f**2 ) / lap
    pres_virial = density * temperature + ( vir + vir_lrc ) / box**3

    IF ( PRESENT ( string ) ) THEN ! output required
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Energy',          energy
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Shifted energy',  energy_sh
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temp-kinet',      temp_kinet
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temp-config',     temp_config
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Virial pressure', pres_virial
    END IF

  END SUBROUTINE calculate

END PROGRAM bd_nvt_lj
