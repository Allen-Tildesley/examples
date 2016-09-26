! md_chain_nve_lj.f90
! Molecular dynamics, NVE ensemble, chain molecule
PROGRAM md_chain_nve_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : lowercase
  USE md_module,        ONLY : model_description, allocate_arrays, deallocate_arrays, &
       &                       force, r, v, n, &
       &                       worst_bond, milcshake_a, milcshake_b, rattle_a, rattle_b

  IMPLICIT NONE

  ! Takes in a configuration of atoms in a linear chain (positions, velocities)
  ! NO periodic boundary conditions, no box
  ! Conducts molecular dynamics with constraints (RATTLE or MILCSHAKE)
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Input configuration, output configuration, all calculations, and all results 
  ! are given in mass = 1 units, and in simulation units defined by the model 
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  ! Most important variables
  REAL :: dt         ! time step
  REAL :: bond       ! bond length
  REAL :: pot        ! total potential energy
  REAL :: kin        ! total kinetic energy
  REAL :: temp_kinet ! kinetic temperature (LJ sigma=1 units, to be averaged)
  REAL :: energy     ! total energy per atom (LJ sigma=1 units, to be averaged)
  REAL :: wc         ! constraint virial (not used in this example)

  INTEGER :: blk, stp, nstep, nblock, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number
  CHARACTER(len=10)           :: constraints

  ! Define procedure pointers with interfaces like those of rattle_a and rattle_b
  PROCEDURE(rattle_a), pointer :: move_a => null()
  PROCEDURE(rattle_b), pointer :: move_b => null()

  NAMELIST /nml/ nblock, nstep, dt, constraints

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_chain'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE ensemble, chain molecule'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL model_description ( output_unit )
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  dt          = 0.002
  constraints = 'rattle'

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_chain'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Time step',                 dt
  IF ( INDEX( lowercase(constraints), 'rattle' ) /= 0 ) THEN
     move_a => rattle_a
     move_b => rattle_b
     WRITE ( unit=output_unit, fmt='(a)' ) 'RATTLE constraint method'
  ELSE IF ( INDEX( lowercase(constraints), 'milcshake' ) /= 0 ) THEN
     move_a => milcshake_a
     move_b => milcshake_b
     WRITE ( unit=output_unit, fmt='(a)' ) 'MILCSHAKE constraint method'
  ELSE
     WRITE ( unit=error_unit, fmt='(a,a)' ) 'Unrecognized constraint method ', constraints
     STOP 'Error in md_chain'
  END IF
  
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond ) ! First call is just to get n and bond
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',          n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Bond length (in sigma units)', bond

  CALL allocate_arrays

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r, v ) ! Second call gets r and v
  WRITE ( unit=output_unit, fmt='(a,t40,es15.5)' ) 'Worst bond length deviation = ', worst_bond ( bond )

  CALL force ( pot )
  call calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Energy', 'Temp-Kinet' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! Velocity Verlet algorithm
        CALL move_a ( dt, bond )     ! RATTLE/MILCSHAKE part A
        CALL force ( pot )           ! Force evaluation
        CALL move_b ( dt, bond, wc ) ! RATTLE/MILCSHAKE part B

        ! Calculate all variables for this step
        call calculate ( )
        CALL blk_add ( [energy,temp_kinet] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL force ( pot )
  call calculate ( 'Final values' )
  WRITE ( unit=output_unit, fmt='(a,t40,es15.5)' ) 'Worst bond length deviation = ', worst_bond ( bond )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r, v )

  CALL deallocate_arrays

CONTAINS

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    kin        = 0.5*SUM(v**2)
    energy     = ( pot + kin ) / REAL ( n )
    temp_kinet = 2.0 * kin / REAL ( 2*(n-1) ) ! NB degrees of freedom = 3(n-1) - (n-1)

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)'  ) 'Total energy', energy
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)'  ) 'Temp-kinet',   temp_kinet
    END IF

  END SUBROUTINE calculate


END PROGRAM md_chain_nve_lj

