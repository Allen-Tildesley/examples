! md_lj_mts.f90
! Molecular dynamics, NVE, multiple timesteps
PROGRAM md_lj_mts

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE md_module,        ONLY : model_description, allocate_arrays, deallocate_arrays, &
       &                       force, r, v, f, n

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using velocity Verlet algorithm
  ! Uses no special neighbour lists, for clarity

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! This program just illustrates the idea of splitting the non-bonded interactions
  ! using a criterion based on distance, for use in a MTS scheme
  ! This would hardly ever be efficient for a simple potential of the LJ kind alone

  ! This program uses mass = 1 throughout
  ! Unlike most other example programs, positions are not divided by box length at all
  ! Input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  ! Most important variables
  REAL :: box         ! box length
  REAL :: density     ! density
  REAL :: dt          ! time step (smallest)
  REAL :: kin         ! total kinetic energy
  REAL :: pres_virial ! virial pressure (to be averaged)
  REAL :: temp_kinet  ! kinetic temperature (to be averaged)
  REAL :: energy      ! total energy per atom (to be averaged)
  REAL :: lambda      ! healing length for switch function

  INTEGER, PARAMETER        :: k_max = 3   ! number of shells
  REAL,    DIMENSION(k_max) :: r_cut       ! cutoff distance for each shell
  REAL,    DIMENSION(k_max) :: pot         ! total potential energy for each shell
  REAL,    DIMENSION(k_max) :: vir         ! total virial for each shell
  INTEGER, DIMENSION(k_max) :: n_mts       ! successive ratios of number of steps for each shell

  INTEGER :: blk, stp1, stp2, stp3, nstep, nblock, k, ioerr
  REAL    :: pairs

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

  NAMELIST /nml/ nblock, nstep, r_cut, lambda, dt, n_mts

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_lj_mts'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE ensemble, multiple time steps'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL model_description ( output_unit )
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock = 10
  nstep  = 1000
  r_cut  = [ 2.4, 3.5, 4.0 ]
  n_mts  = [ 1, 4, 2 ]
  dt     = 0.002
  lambda = 0.1

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_lj_mts'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'      ) 'Number of blocks',           nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'      ) 'Number of steps per block',  nstep
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.5))' ) 'Potential cutoff distances', r_cut

  DO k = 1, k_max
     IF ( k == 1 ) THEN
        pairs = r_cut(k)**3
     ELSE
        pairs = r_cut(k)**3 - r_cut(k-1)**3
        IF ( r_cut(k)-r_cut(k-1) < lambda ) THEN
           WRITE ( unit=error_unit, fmt='(a,3f15.5)' ) 'r_cut interval error', r_cut(k-1), r_cut(k), lambda
           STOP 'Error in md_lj_mts'
        END IF
     END IF
     pairs = REAL(n*(n-1)/2) * (4.0/3.0)*pi * pairs / box**3
     WRITE ( unit=output_unit, fmt='(a,i1,t40,i15)' ) 'Estimated pairs in shell ', k, NINT ( pairs )
  END DO

  WRITE ( unit=output_unit, fmt='(a,t40,*(i15))'  ) 'Multiple step ratios', n_mts(:)
  IF ( n_mts(1) /= 1 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)' ) 'n_mts(1) must be 1', n_mts(1)
     STOP 'Error in md_lj_mts'
  END IF
  IF ( ANY ( n_mts <= 0 ) ) THEN
     WRITE ( unit=error_unit, fmt='(a,*(i15))' ) 'n_mts values must be positive', n_mts
     STOP 'Error in md_lj_mts'
  END IF
  DO k = 1, k_max
     WRITE ( unit=output_unit, fmt='(a,i1,t40,f15.5)' ) 'Time step for shell ', k, PRODUCT(n_mts(1:k))*dt
  END DO

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box ** 3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density
  IF ( r_cut(k_max) > box/2.0  ) THEN
     WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_cut(k_max) too large ', r_cut(k_max)
     STOP 'Error in md_lj_mts'
  END IF

  CALL allocate_arrays ( r_cut )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call to get r and v

  r(:,:) = r(:,:) - ANINT ( r(:,:) / box ) * box ! Periodic boundaries

  ! Calculate forces and pot, vir contributions for each shell
  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, pot(k), vir(k) )
  END DO
  kin         = 0.5*SUM(v**2)
  energy      = ( SUM(pot) + kin ) / REAL ( n )
  temp_kinet  = 2.0 * kin / REAL ( 3*(n-1) )
  pres_virial = density * temp_kinet + SUM(vir) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial total energy',    energy
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial temp-kinet',      temp_kinet
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial virial pressure', pres_virial

  CALL run_begin ( [ CHARACTER(len=15) :: 'Energy', 'Temp-Kinet', 'Virial Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     ! The following set of nested loops is specific to k_max=3

     DO stp3 = 1, nstep ! Begin loop over steps

        ! Outer shell 3: a single step of size n_mts(3) * n_mts(2) * dt
        v(:,:) = v(:,:) + 0.5 * n_mts(3) * n_mts(2) * dt * f(:,:,3) ! Kick half-step (outer shell)

        DO stp2 = 1, n_mts(3) ! Middle shell 2: n_mts(3) steps of size n_mts(2) * dt

           v(:,:) = v(:,:) + 0.5 * n_mts(2) * dt * f(:,:,2) ! Kick half-step (middle shell)

           DO stp1 = 1, n_mts(2) ! Inner shell 1: n_mts(3) * n_mts(2) steps of size dt

              v(:,:) = v(:,:) + 0.5 * dt * f(:,:,1)                ! Kick half-step (inner shell)
              r(:,:) = r(:,:) + dt * v(:,:)                        ! Drift step
              r(:,:) = r(:,:) - ANINT ( r(:,:)/box ) * box         ! Periodic boundaries
              CALL force ( box, r_cut, lambda, 1, pot(1), vir(1) ) ! Force evaluation (inner shell)
              v(:,:) = v(:,:) + 0.5 * dt * f(:,:,1)                ! Kick half-step (inner shell)

           END DO ! End inner shell 1

           CALL force ( box, r_cut, lambda, 2, pot(2), vir(2) ) ! Force evaluation (middle shell)
           v(:,:) = v(:,:) + 0.5 * n_mts(2) * dt * f(:,:,2)     ! Kick half-step (middle shell)

        END DO ! End middle shell 2

        CALL force ( box, r_cut, lambda, 3, pot(3), vir(3) )        ! Force evaluation (outer shell)
        v(:,:) = v(:,:) + 0.5 * n_mts(3) * n_mts(2) * dt * f(:,:,3) ! Kick half-step (outer shell)
        ! End outer shell 3

        kin         = 0.5*SUM(v**2)
        energy      = ( SUM(pot) + kin ) / REAL ( n )
        temp_kinet  = 2.0 * kin / REAL ( 3*(n-1) )
        pres_virial = density * temp_kinet + SUM(vir) / box**3

        ! Calculate all variables for this step
        CALL blk_add ( [energy,temp_kinet,pres_virial] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk           ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, pot(k), vir(k) )
  END DO
  kin         = 0.5*SUM(v**2)
  energy      = ( SUM(pot) + kin ) / REAL ( n )
  temp_kinet  = 2.0 * kin / REAL ( 3*(n-1) )
  pres_virial = density * temp_kinet + SUM(vir) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final total energy',    energy
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final temp-kinet',      temp_kinet
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final virial pressure', pres_virial
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r, v )

  CALL deallocate_arrays

END PROGRAM md_lj_mts

