! md_chain.f90
! Molecular dynamics, NVE ensemble, WCA Lennard-Jones chain
PROGRAM md_chain
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit
  USE utility_module,  ONLY : read_cnf_atoms, write_cnf_atoms, time_stamp, lowercase, &
       &                      run_begin, run_end, blk_begin, blk_end, blk_add
  USE md_chain_module, ONLY : allocate_arrays, deallocate_arrays, check_constraints, force, &
       &                      milcshake_a, milcshake_b, rattle_a, rattle_b, r, v, n
  IMPLICIT NONE

  ! Takes in a configuration of atoms in a linear chain (positions, velocities)
  ! NO periodic boundary conditions, no box
  ! Conducts molecular dynamics with constraints (RATTLE or MILCSHAKE)
  ! Uses no special neighbour lists

  ! However, input configuration, output configuration,
  ! all calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1, mass = 1

  ! Most important variables
  REAL :: dt          ! time step
  REAL :: bond        ! bond length
  REAL :: pot         ! total potential energy
  REAL :: kin         ! total kinetic energy
  REAL :: temperature ! temperature (LJ sigma=1 units, to be averaged)
  REAL :: energy      ! total energy per atom (LJ sigma=1 units, to be averaged)
  REAL :: wc          ! constraint virial (not used in this example)

  INTEGER :: blk, stp, nstep, nblock

  CHARACTER(len=12), PARAMETER :: cnf_prefix = 'md_chain.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number
  CHARACTER(len=10)            :: constraints
  PROCEDURE(rattle_a), pointer :: move_a => null()
  PROCEDURE(rattle_b), pointer :: move_b => null()

  NAMELIST /params/ nblock, nstep, dt, constraints

  WRITE(*,'(''md_chain'')')
  WRITE(*,'(''Molecular dynamics, constant-NVE, repulsive Lennard-Jones chain'')')
  WRITE(*,'(''Results in units epsilon = sigma = 1'')')
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  dt          = 0.002
  constraints = 'rattle'

  READ(*,nml=params)
  WRITE(*,'(''Number of blocks'',         t40,i15)'  ) nblock
  WRITE(*,'(''Number of steps per block'',t40,i15)'  ) nstep
  WRITE(*,'(''Time step'',                t40,f15.5)') dt
  IF ( INDEX( lowercase(constraints), 'rattle' ) /= 0 ) THEN
     move_a => rattle_a
     move_b => rattle_b
     WRITE(*,'(a)') 'RATTLE constraint method'
  ELSE IF ( INDEX( lowercase(constraints), 'milcshake' ) /= 0 ) THEN
     move_a => milcshake_a
     move_b => milcshake_b
     WRITE(*,'(a)') 'MILCSHAKE constraint method'
  ELSE
     WRITE(*,'(a,t40,a)') 'Unrecognized constraint method', constraints
     STOP
  END IF
  
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond )
  WRITE(*,'(''Number of particles'', t40,i15)'          ) n
  WRITE(*,'(''Bond length (in sigma units)'',t40,f15.5)') bond

  CALL allocate_arrays

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r, v )
  CALL check_constraints ( bond )

  CALL force ( pot )
  kin         = 0.5*SUM(v**2)
  energy      = ( pot + kin ) / REAL ( n )
  temperature = 2.0 * kin / REAL ( 2*(n-1) ) ! NB degrees of freedom = 3(n-1) - (n-1)
  WRITE(*,'(''Initial total energy (sigma units)'',  t40,f15.5)') energy
  WRITE(*,'(''Initial temperature (sigma units)'',   t40,f15.5)') temperature

  CALL run_begin ( [ CHARACTER(len=15) :: 'Energy', 'Temperature' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! Velocity Verlet algorithm
        CALL move_a ( dt, bond )     ! RATTLE/MILCSHAKE part A
        CALL force ( pot )           ! Force evaluation
        CALL move_b ( dt, bond, wc ) ! RATTLE/MILCSHAKE part B

        kin         = 0.5*SUM(v**2)
        energy      = ( pot + kin ) / REAL ( n )
        temperature = 2.0 * kin / REAL ( 2*(n-1) ) ! NB degrees of freedom = 3(n-1) - (n-1)

        ! Calculate all variables for this step
        CALL blk_add ( [energy,temperature] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL force ( pot )
  kin         = 0.5*SUM(v**2)
  energy      = ( pot + kin ) / REAL ( n )
  temperature = 2.0 * kin / REAL ( 2*(n-1) ) ! NB degrees of freedom = 3(n-1) - (n-1)
  WRITE(*,'(''Final total energy (sigma units)'',  t40,f15.5)') energy
  WRITE(*,'(''Final temperature (sigma units)'',   t40,f15.5)') temperature
  CALL check_constraints ( bond )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r, v )

  CALL deallocate_arrays

END PROGRAM md_chain

