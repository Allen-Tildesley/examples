! md_nvt_lj_le.f90
! MD, NVT ensemble, LJ atoms, Lees-Edwards boundaries
PROGRAM md_nvt_lj_le

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE md_lj_le_module,  ONLY : initialize, finalize, force, r, v, f, n, energy_lrc

  IMPLICIT NONE

  ! MPA TODO: double check this program!!!
  
  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions, with Lees-Edwards shear
  ! Conducts molecular dynamics, SLLOD algorithm, with isokinetic thermostat
  ! Refs: Pan et al J Chem Phys 122 094114 (2005)
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Box is taken to be of unit length during the dynamics
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1, mass = 1
  ! property      program units            Lennard-Jones units
  ! length        box                      sigma
  ! mass          m                        m
  ! energy        epsilon                  epsilon
  ! time          box*sqrt(m/epsilon)      sigma*sqrt(m/epsilon)
  ! velocity      sqrt(epsilon/m)          sqrt(epsilon/m)
  ! pressure      epsilon/box**3           epsilon/sigma**3
  ! strain        dimensionless            dimensionless
  ! strain_rate   sqrt(epsilon/m)/box      sqrt(epsilon/m)/sigma

  ! Most important variables
  REAL :: sigma       ! atomic diameter (in units where box=1)
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: dt          ! time step
  REAL :: strain_rate ! strain_rate (velocity gradient) dv_x/dr_y
  REAL :: strain      ! strain (integrated velocity gradient) dr_x/dr_y
  REAL :: r_cut       ! potential cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: pot_sh      ! total shifted potential energy
  REAL :: kin         ! total kinetic energy
  REAL :: vir         ! total virial
  REAL :: pressure    ! pressure (LJ sigma=1 units, to be averaged)
  REAL :: temperature ! temperature (LJ sigma=1 units, to be averaged)
  REAL :: energy      ! total energy per atom (LJ sigma=1 units, to be averaged)
  REAL :: energy_sh   ! total shifted energy per atom (LJ sigma=1 units, to be averaged)

  INTEGER :: blk, stp, nstep, nblock, ioerr
  REAL    :: pot_lrc, vir_lrc
  REAL    :: c1, c2, alpha, beta, e, h, d_strain, dt_factor, prefactor

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, r_cut, dt, strain_rate

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_nvt_lj_le'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVT, Lennard-Jones, Lees-Edwards'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Results in units epsilon = sigma = 1'
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  r_cut       = 2.5
  dt          = 0.005
  strain_rate = 0.01

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_nvt_lj_le'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Strain rate',               strain_rate

  ! For simplicity we assume that input configuration has strain = 0
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box (in sigma units)', box
  sigma = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Reduced density', density

  ! Convert run and potential parameters to box units
  sigma       = sigma / box
  r_cut       = r_cut / box
  dt          = dt / box
  strain_rate = strain_rate * box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'sigma  (in box units)',      sigma
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'r_cut  (in box units)',      r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'dt     (in box units)',      dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'strain rate (in box units)', strain_rate
  IF ( r_cut > 0.5  ) THEN
     WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_cut too large ', r_cut
     STOP 'Error in md_nvt_lj_le'
  END IF

  CALL initialize ( r_cut )

  ! For simplicity we assume that input configuration has strain = 0
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v
  strain = 0.0

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(1,:) = r(1,:) - ANINT ( r(2,:) ) * strain ! Extra correction
  r(:,:) = r(:,:) - ANINT ( r(:,:) )          ! Periodic boundaries

  CALL force ( sigma, r_cut, strain, pot, pot_sh, vir )
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
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

        ! Isokinetic SLLOD algorithm (Pan et al)

        ! Operator A for half a step
        d_strain = 0.5 * dt * strain_rate
        r(1,:)   = r(1,:) + d_strain * r(2,:)  ! Extra strain term
        r(:,:)   = r(:,:) + 0.5 * dt * v(:,:)  ! Drift half-step
        strain   = strain + d_strain           ! Advance boundaries

        r(1,:) = r(1,:) - ANINT ( r(2,:) ) * strain   ! Extra correction
        r(:,:) = r(:,:) - ANINT ( r(:,:) )            ! Periodic boundaries

        ! Operator B1 for half a step
        d_strain = 0.5*dt * strain_rate
        c1       = d_strain * SUM ( v(1,:)*v(2,:) ) / SUM ( v(:,:)**2 )
        c2       = ( d_strain**2 ) * SUM ( v(2,:)**2 ) / SUM ( v(:,:)**2 )
        v(1,:)   = v(1,:) - d_strain*v(2,:)
        v(:,:)   = v(:,:) / SQRT ( 1.0 - 2.0*c1 + c2 )
        
        CALL force ( sigma, r_cut, strain, pot, pot_sh, vir ) ! Force evaluation

        ! Operator B2 for a full step
        alpha     = SUM ( f(:,:)*v(:,:) ) / SUM ( v(:,:)**2 )
        beta      = SQRT ( SUM ( f(:,:)**2 ) / SUM ( v(:,:)**2 ) )
        h         = ( alpha + beta ) / ( alpha - beta )
        e         = EXP ( -beta * dt )
        dt_factor = ( 1 + h - e - h / e ) / ( ( 1 - h ) * beta )
        prefactor = ( 1 - h ) / ( e - h / e )
        v(:,:)    = prefactor * ( v(:,:) + dt_factor * f(:,:) )

        ! Operator B1 for half a step
        d_strain = 0.5*dt * strain_rate
        c1       = d_strain * SUM ( v(1,:)*v(2,:) ) / SUM ( v(:,:)**2 )
        c2       = (d_strain**2 ) * SUM ( v(2,:)**2 ) / SUM ( v(:,:)**2 )
        v(1,:)   = v(1,:) - d_strain*v(2,:)
        v(:,:)   = v(:,:) / SQRT ( 1.0 - 2.0*c1 + c2 )

        ! Operator A for half a step
        d_strain = 0.5 * dt * strain_rate
        r(1,:)   = r(1,:) + d_strain * r(2,:)  ! Extra strain term
        r(:,:)   = r(:,:) + 0.5 * dt * v(:,:)  ! Drift half-step
        strain   = strain + d_strain           ! Advance boundaries

        r(1,:) = r(1,:) - ANINT ( r(2,:) ) * strain   ! Extra correction
        r(:,:) = r(:,:) - ANINT ( r(:,:) )            ! Periodic boundaries

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

  CALL force ( sigma, r_cut, strain, pot, pot_sh, vir )
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
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

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )
  CALL time_stamp ( output_unit )

  CALL finalize

END PROGRAM md_nvt_lj_le

