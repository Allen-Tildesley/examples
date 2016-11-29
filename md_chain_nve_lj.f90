! md_chain_nve_lj.f90
! Molecular dynamics, NVE ensemble, chain molecule
PROGRAM md_chain_nve_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE maths_module,     ONLY : lowercase
  USE md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       zero_cm, force, worst_bond, r, v, n, potential_type, &
       &                       milcshake_a, milcshake_b, rattle_a, rattle_b

  IMPLICIT NONE

  ! Takes in a configuration of atoms in a linear chain (positions, velocities)
  ! NO periodic boundary conditions, no box
  ! Conducts molecular dynamics with constraints (RATTLE or MILC-SHAKE)
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Input configuration, output configuration, all calculations, and all results 
  ! are given in mass = 1 units, and in simulation units defined by the model 
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  ! Most important variables
  REAL :: dt   ! Time step
  REAL :: bond ! Bond length
  REAL :: wc   ! Constraint virial (not used in this example)

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  ! Composite interaction = pot & ovr variables
  TYPE(potential_type) :: total

  INTEGER :: blk, stp, nstep, nblock, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number
  CHARACTER(len=30)           :: constraints

  ! Define procedure pointers with interfaces like those of rattle_a and rattle_b
  PROCEDURE(rattle_a), POINTER :: move_a => NULL()
  PROCEDURE(rattle_b), POINTER :: move_b => NULL()

  NAMELIST /nml/ nblock, nstep, dt, constraints

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_chain'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE ensemble, chain molecule'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Uses a potential shifted to vanish smoothly at cutoff, no LRC'
  WRITE ( unit=output_unit, fmt='(a)' ) 'No periodic boundaries'
  CALL introduction

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  dt          = 0.002
  constraints = 'rattle'

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_chain'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step',                 dt
  IF ( INDEX( lowercase(constraints), 'ratt' ) /= 0 ) THEN
     move_a => rattle_a
     move_b => rattle_b
     WRITE ( unit=output_unit, fmt='(a)' ) 'RATTLE constraint method'
  ELSE IF ( INDEX( lowercase(constraints), 'milc' ) /= 0 ) THEN
     move_a => milcshake_a
     move_b => milcshake_b
     WRITE ( unit=output_unit, fmt='(a)' ) 'MILCSHAKE constraint method'
  ELSE
     WRITE ( unit=error_unit, fmt='(a,a)' ) 'Unrecognized constraint method ', constraints
     STOP 'Error in md_chain'
  END IF

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond ) ! First call is just to get n and bond
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',          n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Bond length (in sigma units)', bond
  CALL allocate_arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r, v ) ! Second call gets r and v
  CALL zero_cm ! Set centre-of-mass position and velocity to zero
  WRITE ( unit=output_unit, fmt='(a,t40,es15.5)' ) 'Worst bond length deviation = ', worst_bond ( bond )

  ! Calculate initial values and forces
  CALL force ( total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in md_chain_nve_lj'
  END IF
  CALL calculate ( 'Initial values' )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( variables )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL move_a ( dt, bond )     ! RATTLE/MILCSHAKE part A

        CALL force ( total )         ! Force evaluation
        IF ( total%ovr ) THEN
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Error in md_chain_nve_lj'
        END IF

        CALL move_b ( dt, bond, wc ) ! RATTLE/MILCSHAKE part B

        ! Calculate and accumulate variables for this step
        CALL calculate ( )
        CALL blk_add ( variables )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                        ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ! Output run averages

  CALL force ( total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in md_chain_nve_lj'
  END IF
  CALL calculate ( 'Final values' )
  WRITE ( unit=output_unit, fmt='(a,t40,es15.5)' ) 'Worst bond length deviation = ', worst_bond ( bond )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r, v ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE calculate ( string )
    USE averages_module, ONLY : write_variables
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    ! In this example we simulate using a specified potential (e.g. WCA LJ)
    ! which goes to zero smoothly at the cutoff to highlight energy conservation
    ! so long-range corrections do not arise

    TYPE(variable_type) :: e_f, t_k
    REAL                :: kin

    ! Preliminary calculations
    kin  = 0.5*SUM(v**2)

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Internal energy per atom
    ! Total KE plus total LJ nonbonded energy divided by N
    e_f = variable_type ( nam = 'E/N', val = (kin+total%pot)/REAL(n) )

    ! Kinetic temperature
    ! Remove 6 degrees of freedom for conserved linear and angular momentum
    ! and also (n-1) degrees of freedom for the bond constraints
    t_k = variable_type ( nam = 'T kinetic', val = 2.0*kin/REAL(3*n-(n-1)-6) )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ e_f, t_k ]

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       CALL write_variables ( variables )
    END IF

  END SUBROUTINE calculate


END PROGRAM md_chain_nve_lj

