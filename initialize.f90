! initialize.f90
! Sets up initial configuration for MD or MC
PROGRAM initialize

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module,  ONLY : write_cnf_atoms, write_cnf_mols
  USE maths_module,      ONLY : lowercase
  USE initialize_module, ONLY : allocate_arrays, deallocate_arrays, &
       &                        initialize_positions_lattice, initialize_orientations_lattice, &
       &                        initialize_positions_random,  initialize_orientations_random, &
       &                        initialize_chain_lattice, initialize_chain_random, initialize_chain_velocities, &
       &                        initialize_velocities, initialize_angular_velocities, &
       &                        n, r, e, v, w

  IMPLICIT NONE

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Most important variables
  INTEGER :: nc                  ! Number of fcc unit cells in each coordinate direction; n = 4*nc**3
  REAL    :: temperature         ! Specified temperature
  REAL    :: inertia             ! Specified moment of inertia (for linear and nonlinear molecules)
  REAL    :: density             ! Specified number density
  REAL    :: bond                ! Specified bond length (for chain molecules)
  LOGICAL :: velocities          ! User option requiring velocities for MD (otherwise just positions)
  LOGICAL :: random_positions    ! User option for random positions (otherwise on a lattice)
  LOGICAL :: random_orientations ! User option for random orientations (otherwise on a lattice)

  CHARACTER(len=10)  :: molecules       ! Character string, input, used to specify molecule_option
  INTEGER            :: molecule_option ! User option for atoms, linear, nonlinear or chain molecule
  INTEGER, PARAMETER :: atoms = 0, linear = 1, nonlinear = 2, chain = 3

  INTEGER            :: ioerr
  REAL               :: box ! Deduced from N and density (for atoms, linear, and nonlinear molecules)

  CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp' ! Will be used as an input file by later simulations

  NAMELIST /nml/ nc, n, temperature, inertia, density, bond, &
       &         velocities, molecules, random_positions, random_orientations

  WRITE ( unit=output_unit, fmt='(a)' ) 'initialize'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Sets up initial configuration file for various simulations'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Options for molecules are "atoms", "chain", "linear", "nonlinear"'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass m=1 throughout'

  CALL RANDOM_SEED()

  ! Set default parameters
  n                   = 0       ! nc takes precedence unless n is explicitly specified
  nc                  = 4       ! Default is N = 4*(4**3) = 256 on a fcc lattice, a small system
  temperature         = 1.0     ! Should lie in the liquid region for density > 0.7 or so
  inertia             = 1.0     ! Only relevant for linear and nonlinear molecules
  density             = 0.75    ! Should lie in the liquid region for temperature > 0.9 or so
  bond                = 1.0     ! Only relevant for chains
  velocities          = .FALSE. ! By default, produce positions only, for MC simulations
  molecules           = 'atoms' ! Options are 'atoms', 'chain', 'linear', 'nonlinear'
  random_positions    = .FALSE. ! By default, arrange atoms on a lattice
  random_orientations = .FALSE. ! By default, use predetermined molecular orientations

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in initialize'
  END IF

  IF ( n <= 0 ) THEN ! Test for unspecified N

     IF ( nc <= 0 ) THEN ! Test for unphysical nc

        WRITE ( unit=error_unit, fmt='(a,2i15)') 'nc must be positive', nc
        STOP 'Error in initialize'

     ELSE ! Deduce N from nc

        WRITE ( unit=output_unit, fmt='(a,t40,i15)') 'nc = ', nc
        n = 4*nc**3
        WRITE ( unit=output_unit, fmt='(a,t40,i15)') 'n = ', n

     END IF

  ELSE ! N has been specified

     WRITE ( unit=output_unit, fmt='(a,t40,i15)') 'n = ', n

  END IF ! End test for unspecified N

  ! Use molecules string to deduce molecule_option
  molecules = lowercase(molecules)
  
  IF ( INDEX ( molecules, 'chain' ) /= 0 ) THEN

     molecule_option = chain
     WRITE ( unit=output_unit, fmt='(a)' ) 'Chain of atoms, no periodic boundaries'

  ELSE IF ( INDEX ( molecules, 'nonlinear' ) /= 0 ) THEN

     molecule_option = nonlinear
     WRITE ( unit=output_unit, fmt='(a)' ) 'Nonlinear molecules, periodic boundaries'

  ELSE IF ( INDEX ( molecules, 'linear') /= 0 ) THEN

     molecule_option = linear
     WRITE ( unit=output_unit, fmt='(a)' ) 'Linear molecules, periodic boundaries'

  ELSE IF ( INDEX ( molecules, 'atoms') /= 0 ) THEN

     molecule_option = atoms
     WRITE ( unit=output_unit, fmt='(a)' ) 'Atoms, periodic boundaries'

  ELSE

     WRITE ( unit=error_unit, fmt='(a,a)') 'Unrecognized molecules option: ', molecules
     STOP 'Error in initialize'

  END IF

  ! Allocate arrays appropriately according to molecule_option
  CALL allocate_arrays ( quaternions = ( molecule_option == nonlinear ) )

  ! Initialize positions and optionally velocities
  
  SELECT CASE ( molecule_option)

  CASE ( chain )

     IF ( random_positions ) THEN
        CALL initialize_chain_random ! Random with unit bond length
     ELSE
        CALL initialize_chain_lattice ! On a lattice with unit bond length
     END IF

     IF ( velocities ) CALL initialize_chain_velocities ( temperature )

  CASE default

     IF ( random_positions ) THEN
        CALL initialize_positions_random ! Random within unit box
     ELSE
        CALL initialize_positions_lattice ! On a lattice within unit box
     END IF

     IF ( velocities ) CALL initialize_velocities ( temperature )

  END SELECT

  ! For linear and nonlinear molecules, initialize orientations and optionally angular velocities
  
  SELECT CASE ( molecule_option )

  CASE ( linear, nonlinear )

     IF ( random_orientations ) THEN
        CALL initialize_orientations_random
     ELSE
        CALL initialize_orientations_lattice
     END IF

     IF ( velocities ) CALL initialize_angular_velocities ( temperature, inertia )

  END SELECT

  ! Write out configuration
  
  WRITE ( unit=output_unit, fmt='(a,a)' ) 'Writing configuration to filename ', filename

  SELECT CASE ( molecule_option )

  CASE ( atoms )

     ! Write out coordinates in same units as box
     box = ( REAL(n) / density ) ** ( 1.0/3.0 )
     WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density',    density
     WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box length', box
     IF ( velocities ) THEN
        CALL write_cnf_atoms ( filename, n, box, box*r, v )
     ELSE
        CALL write_cnf_atoms ( filename, n, box, box*r )
     END IF

  CASE ( chain )

     ! We do not use periodic boundaries for this system
     ! Instead, use "box" variable to store bond length
     WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Bond length', bond
     IF ( velocities ) THEN
        CALL write_cnf_atoms ( filename, n, bond, bond*r, v )
     ELSE
        CALL write_cnf_atoms ( filename, n, bond, bond*r )
     END IF

  CASE ( linear, nonlinear )

     ! Write out coordinates in same units as box
     box = ( REAL(n) / density ) ** ( 1.0/3.0 )
     WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density',    density
     WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box length', box
     IF ( velocities ) THEN
        CALL write_cnf_mols ( filename, n, box, box*r, e, v, w )
     ELSE
        CALL write_cnf_mols ( filename, n, box, box*r, e )
     END IF

  END SELECT

  CALL deallocate_arrays

END PROGRAM initialize
