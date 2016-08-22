! bd_nvt_lj.f90
! Brownian dynamics, NVT ensemble, Lennard-Jones atoms
PROGRAM bd_nvt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,   ONLY : random_normals
  USE md_lj_module,     ONLY : allocate_arrays, deallocate_arrays, force, r, v, f, n, energy_lrc

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using BAOAB algorithm of BJ Leimkuhler and C Matthews
  ! Appl. Math. Res. eXpress 2013, 34–56 (2013); J. Chem. Phys. 138, 174102 (2013)
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1, mass = 1

  ! Most important variables
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: dt          ! time step
  REAL :: r_cut       ! potential cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: pot_sh      ! total shifted potential energy
  REAL :: kin         ! total kinetic energy
  REAL :: vir         ! total virial
  REAL :: pressure    ! pressure (LJ sigma=1 units, to be averaged)
  REAL :: temperature ! temperature (LJ sigma=1 units, specified)
  REAL :: energy      ! total energy per atom (LJ sigma=1 units, to be averaged)
  REAL :: energy_sh   ! total shifted energy per atom (LJ sigma=1 units, to be averaged)
  REAL :: gamma       ! friction coefficient

  INTEGER :: blk, stp, nstep, nblock, ioerr
  REAL    :: pot_lrc, vir_lrc, c0, c1

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, r_cut, dt, gamma, temperature

  WRITE ( unit=output_unit, fmt='(a)' ) 'bd_nvt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Brownian dynamics, constant-NVT, Lennard-Jones'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Results in units epsilon = sigma = 1'
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  r_cut       = 2.5
  dt          = 0.005
  gamma       = 0.1
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

  ! Dimensionless coefficients for random momentum term
  gamma = gamma * dt ! incorporate timestep into friction
  c0 = EXP(-gamma)
  IF ( gamma > 1.e-5 ) THEN
     c1 = 1-c0**2
  ELSE
     c1 = gamma * ( 2.0 - gamma * ( 2.0 - 4.0 * gamma / 3.0 ) ) ! Taylor expansion for low gamma
  END IF
  c1 = SQRT ( c1 )
  WRITE ( unit=output_unit, fmt='(a,t40,es15.4)' ) 'Algorithm coefficient c0', c0
  WRITE ( unit=output_unit, fmt='(a,t40,es15.4)' ) 'Algorithm coefficient c1', c1

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box (in sigma units)', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Reduced density', density

  CALL allocate_arrays ( box, r_cut )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL force ( box, r_cut, pot, pot_sh, vir )
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  pot         = pot + pot_lrc
  vir         = vir + vir_lrc
  kin         = 0.5*SUM(v**2)
  energy      = ( pot + kin ) / REAL ( n )
  energy_sh   = ( pot_sh + kin ) / REAL ( n )
  temperature = 2.0 * kin / REAL ( 3*(n-1) )
  pressure    = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial total energy (sigma units)',   energy
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial shifted energy (sigma units)', energy_sh
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial temperature (sigma units)',    temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial pressure (sigma units)',       pressure

  CALL run_begin ( [ CHARACTER(len=15) :: 'Energy', 'Shifted Energy', 'Temperature', 'Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! BAOAB algorithm
        v(:,:) = v(:,:) + 0.5 * dt * f(:,:)               ! Kick half-step (B)
        r(:,:) = r(:,:) + 0.5 * dt * v(:,:) / box         ! Drift half-step (A) (positions in box=1 units)
        CALL random_normals ( 0.0, SQRT(temperature), f ) ! Use f to hold random momenta
        v(:,:) = c0 * v(:,:) + c1 * f                     ! Friction and random effects (O)
        r(:,:) = r(:,:) + 0.5 * dt * v(:,:) / box         ! Drift half-step (A) (positions in box=1 units)
        r(:,:) = r(:,:) - ANINT ( r(:,:) )                ! Periodic boundaries (box=1 units)
        CALL force ( box, r_cut, pot, pot_sh, vir )       ! Force evaluation
        v(:,:) = v(:,:) + 0.5 * dt * f(:,:)               ! Kick half-step (B)

        CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
        pot         = pot + pot_lrc
        vir         = vir + vir_lrc
        kin         = 0.5*SUM(v**2)
        energy      = ( pot + kin ) / REAL ( n )
        energy_sh   = ( pot_sh + kin ) / REAL ( n )
        temperature = 2.0 * kin / REAL ( 3*(n-1) )
        pressure    = density * temperature + vir / box**3

        ! Calculate all variables for this step
        CALL blk_add ( [energy,energy_sh,temperature,pressure] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL force ( box, r_cut, pot, pot_sh, vir )
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  pot         = pot + pot_lrc
  vir         = vir + vir_lrc
  kin         = 0.5*SUM(v**2)
  energy      = ( pot + kin ) / REAL ( n )
  energy_sh   = ( pot_sh + kin ) / REAL ( n )
  temperature = 2.0 * kin / REAL ( 3*(n-1) )
  pressure    = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final total energy (sigma units)',   energy
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final shifted energy (sigma units)', energy_sh
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final temperature (sigma units)',    temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final pressure (sigma units)',       pressure
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )

  CALL deallocate_arrays

END PROGRAM bd_nvt_lj

