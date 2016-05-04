! initialize.f90
! Sets up initial configuration for MD or MC
PROGRAM initialize
  USE utility_module,    ONLY : write_cnf_atoms, write_cnf_mols, lowercase
  USE initialize_module, ONLY : allocate_arrays, deallocate_arrays, &
       &                        initialize_positions_lattice, initialize_orientations_lattice, &
       &                        initialize_positions_random,  initialize_orientations_random, &
       &                        initialize_chain_lattice, initialize_chain_random, initialize_chain_velocities, &
       &                        initialize_velocities, initialize_angular_velocities, &
       &                        n, r, e, v, w
  IMPLICIT NONE

  INTEGER            :: nc
  REAL               :: temperature, inertia, density, box, bond
  LOGICAL            :: velocities, random_positions, random_orientations
  INTEGER            :: molecule_option
  CHARACTER(len=10)  :: molecules
  INTEGER, PARAMETER :: atoms = 0, linear = 1, nonlinear = 2, chain = 3

  NAMELIST /params/ nc, n, temperature, inertia, density, bond, &
       &            velocities, molecules, random_positions, random_orientations

  ! Default values
  ! We must specify either n or nc but not both
  n           = -1
  nc          = -1
  temperature = 1.0
  inertia     = 1.0
  density     = 0.5
  bond        = 1.0
  velocities  = .FALSE.
  molecules   = 'atoms'
  random_positions = .FALSE.
  random_orientations = .FALSE.

  READ(*,nml=params)

  IF ( n < 0 ) THEN
     IF ( nc < 0 ) THEN
        STOP 'Must specify either n or nc'
     ELSE
        WRITE(*,'(''nc = '',t40,i5)') nc
        n = 4*nc**3 ! fcc lattice
        WRITE(*,'(''n = '',t40,i5)') n
     END IF
  ELSE
     IF ( nc < 0 ) THEN
        WRITE(*,'(''n = '',t40,i5)') n
     ELSE
        STOP 'Must specify either n or nc, not both'
     END IF
  END IF

  IF ( INDEX(lowercase(molecules), 'chain') /= 0 ) THEN
     molecule_option = chain
     WRITE(*,'(a)') 'Chain of atoms, no periodic boundaries'
  ELSE IF  ( INDEX(lowercase(molecules), 'nonlinear') /= 0 ) THEN
     molecule_option = nonlinear
     WRITE(*,'(a)') 'Nonlinear molecules, periodic boundaries'
  ELSE IF  ( INDEX(lowercase(molecules), 'linear') /= 0 ) THEN
     molecule_option = linear
     WRITE(*,'(a)') 'Linear molecules, periodic boundaries'
  ELSE IF  ( INDEX(lowercase(molecules), 'atoms') /= 0 ) THEN
     molecule_option = atoms
     WRITE(*,'(a)') 'Atoms, periodic boundaries'
  ELSE
     STOP 'Unrecognized molecules option'
  END IF

  CALL allocate_arrays ( quaternions = ( molecule_option == nonlinear ) )

  CALL RANDOM_SEED()

  SELECT CASE ( molecule_option)

  CASE ( chain )

     IF ( random_positions ) THEN ! chain bonds are randomly oriented, avoiding overlaps
        CALL initialize_chain_random ! unit bond length
     ELSE ! chain is a close-packed cube of atoms surrounded by vacuum
        CALL initialize_chain_lattice ! unit bond length
     END IF
     IF ( velocities ) THEN
        CALL initialize_chain_velocities ( temperature )
     END IF

  CASE default

     IF ( random_positions ) THEN ! coordinates chosen randomly
        CALL initialize_positions_random ! unit box
     ELSE ! close packed lattice
        CALL initialize_positions_lattice ! unit box
     END IF
     IF ( velocities ) THEN
        CALL initialize_velocities ( temperature )
     END IF

  END SELECT

  SELECT CASE ( molecule_option )

  CASE ( linear, nonlinear )
     IF ( random_orientations ) THEN
        CALL initialize_orientations_random
     ELSE
        CALL initialize_orientations_lattice
     END IF

     IF ( velocities ) THEN
        CALL initialize_angular_velocities ( temperature, inertia )
     END IF

  END SELECT

  SELECT CASE ( molecule_option )

  CASE ( atoms )

     ! Write out coordinates in same units as box
     box = ( REAL(n) / density ) ** ( 1.0/3.0 )
     IF ( velocities ) THEN
        CALL write_cnf_atoms ( 'md.cnf', n, box, box*r, v )
     ELSE
        CALL write_cnf_atoms ( 'mc.cnf', n, box, box*r )
     END IF

  CASE ( chain )

     ! We do not use periodic boundaries for this system
     ! Instead, use "box" variable to store bond length
     IF ( velocities ) THEN
        CALL write_cnf_atoms ( 'md.cnf', n, bond, bond*r, v )
     ELSE
        CALL write_cnf_atoms ( 'mc.cnf', n, bond, bond*r )
     END IF

  CASE ( linear, nonlinear )

     ! Write out coordinates in same units as box
     box = ( REAL(n) / density ) ** ( 1.0/3.0 )
     IF ( velocities ) THEN
        CALL write_cnf_mols ( 'md.cnf', n, box, box*r, e, v, w )
     ELSE
        CALL write_cnf_mols ( 'mc.cnf', n, box, box*r, e )
     END IF

  END SELECT

  CALL deallocate_arrays

END PROGRAM initialize
