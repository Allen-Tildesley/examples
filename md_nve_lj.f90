! md_nve_lj.f90
! Molecular dynamics, NVE ensemble, Lennard-Jones atoms
PROGRAM md_nve_lj
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit
  USE utility_module, ONLY : read_cnf_atoms, write_cnf_atoms, time_stamp, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add
  USE md_lj_module,   ONLY : initialize, finalize, force, r, v, f, n, energy_lrc
  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using velocity Verlet algorithm
  ! Uses no special neighbour lists

  ! Box is taken to be of unit length during the dynamics
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1, mass = 1

  ! Most important variables
  REAL :: sigma       ! atomic diameter (in units where box=1)
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: dt          ! time step
  REAL :: r_cut       ! potential cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: pot_sh      ! total shifted potential energy
  REAL :: kin         ! total kinetic energy
  REAL :: vir         ! total virial
  REAL :: pressure    ! pressure (LJ sigma=1 units, to be averaged)
  REAL :: temperature ! temperature (LJ sigma=1 units, to be averaged)
  REAL :: energy      ! total energy per atom (LJ sigma=1 units, to be averaged)
  REAL :: energy_sh   ! total shifted energy per atom (LJ sigma=1 units, to be averaged)

  INTEGER :: blk, stp, nstep, nblock
  REAL    :: pot_lrc, vir_lrc

  CHARACTER(len=13), PARAMETER :: cnf_prefix = 'md_nve_lj.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /run_parameters/ nblock, nstep, r_cut, dt

  WRITE(*,'(''md_nve_lj'')')
  WRITE(*,'(''Molecular dynamics, constant-NVE, Lennard-Jones'')')
  WRITE(*,'(''Results in units epsilon = sigma = 1'')')
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  r_cut       = 2.5
  dt          = 0.005

  READ(*,nml=run_parameters)
  WRITE(*,'(''Number of blocks'',         t40,i15)'  ) nblock
  WRITE(*,'(''Number of steps per block'',t40,i15)'  ) nstep
  WRITE(*,'(''Potential cutoff distance'',t40,f15.5)') r_cut
  WRITE(*,'(''Time step'',                t40,f15.5)') dt

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box )
  WRITE(*,'(''Number of particles'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density

  ! Convert run and potential parameters to box units
  sigma  = sigma / box
  r_cut  = r_cut / box
  dt     = dt / box
  WRITE(*,'(''sigma  (in box units)'',t40,f15.5)') sigma
  WRITE(*,'(''r_cut  (in box units)'',t40,f15.5)') r_cut
  WRITE(*,'(''dt     (in box units)'',t40,f15.5)') dt
  IF ( r_cut > 0.5  ) STOP 'r_cut too large '

  CALL initialize ( r_cut )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v )

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL force ( sigma, r_cut, pot, pot_sh, vir )
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot         = pot + pot_lrc
  vir         = vir + vir_lrc
  kin         = 0.5*SUM(v**2)
  energy      = ( pot + kin ) / REAL ( n )
  energy_sh   = ( pot_sh + kin ) / REAL ( n )
  temperature = 2.0 * kin / REAL ( 3*(n-1) )
  pressure    = density * temperature + vir / box**3
  WRITE(*,'(''Initial total energy (sigma units)'',  t40,f15.5)') energy
  WRITE(*,'(''Initial shifted energy (sigma units)'',t40,f15.5)') energy_sh
  WRITE(*,'(''Initial temperature (sigma units)'',   t40,f15.5)') temperature
  WRITE(*,'(''Initial pressure (sigma units)'',      t40,f15.5)') pressure

  CALL run_begin ( [ CHARACTER(len=15) :: 'Energy', 'Shifted Energy', 'Temperature', 'Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! Velocity Verlet algorithm
        v(:,:) = v(:,:) + 0.5 * dt * f(:,:)           ! Kick half-step
        r(:,:) = r(:,:) + dt * v(:,:)                 ! Drift step
        r(:,:) = r(:,:) - anint ( r(:,:) )            ! Periodic boundaries
        CALL force ( sigma, r_cut, pot, pot_sh, vir ) ! Force evaluation
        v(:,:) = v(:,:) + 0.5 * dt * f(:,:)           ! Kick half-step

        CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
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

  CALL force ( sigma, r_cut, pot, pot_sh, vir )
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot         = pot + pot_lrc
  vir         = vir + vir_lrc
  kin         = 0.5*SUM(v**2)
  energy      = ( pot + kin ) / REAL ( n )
  energy_sh   = ( pot_sh + kin ) / REAL ( n )
  temperature = 2.0 * kin / REAL ( 3*(n-1) )
  pressure    = density * temperature + vir / box**3
  WRITE(*,'(''Final total energy (sigma units)'',  t40,f15.5)') energy
  WRITE(*,'(''Final shifted energy (sigma units)'',t40,f15.5)') energy_sh
  WRITE(*,'(''Final temperature (sigma units)'',   t40,f15.5)') temperature
  WRITE(*,'(''Final pressure (sigma units)'',      t40,f15.5)') pressure
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )

  CALL finalize

END PROGRAM md_nve_lj

