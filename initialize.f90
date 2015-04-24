! initialize.f90
! Sets up initial configuration for MD or MC
PROGRAM initialize
  USE utility_module,    ONLY : write_cnf_atoms, write_cnf_molecules
  USE initialize_module, ONLY : initialize_positions, initialize_orientations, &
       &                        initialize_velocities, initialize_angular_velocities, &
       &                        n, r, e, v, w
  IMPLICIT NONE

  INTEGER :: nc
  REAL    :: temperature, inertia, density, box
  LOGICAL :: velocities, molecules

  NAMELIST /initialize_parameters/ nc, n, temperature, inertia, density, velocities, molecules

  ! Default values
  ! We must specify either n or nc but not both
  n           = -1
  nc          = -1
  temperature = 1.0
  inertia     = 1.0
  density     = 0.5
  velocities  = .FALSE.
  molecules   = .FALSE.

  READ(*,nml=initialize_parameters)

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
        nc = NINT ( ( REAL(n)/4.0 )**(1.0/3.0) )
        IF ( n /= 4*nc**3 ) STOP 'n not of the form  4*nc**3'
        WRITE(*,'(''nc = '',t40,i5)') nc
     ELSE
        STOP 'Must specify either n or nc, not both'
     END IF
  END IF

  CALL RANDOM_SEED()

  CALL initialize_positions ( nc ) ! coordinates initialized in unit box
  IF ( molecules ) THEN
     CALL initialize_orientations ( nc )
  END IF
  IF ( velocities ) THEN
     CALL initialize_velocities ( temperature )
     IF (molecules ) THEN
        CALL initialize_angular_velocities ( temperature, inertia )
     END IF
  END IF

  box = ( REAL(n) / density ) ** ( 1.0/3.0 )

  ! Write out coordinates in same units as box
  IF ( molecules ) THEN
     IF ( velocities ) THEN
        CALL write_cnf_molecules ( 'initialize.cnf', n, box, box*r, e, v, w )
     ELSE
        CALL write_cnf_molecules ( 'initialize.cnf', n, box, box*r, e )
     END IF
  ELSE
     IF ( velocities ) THEN
        CALL write_cnf_atoms ( 'initialize.cnf', n, box, box*r, v )
     ELSE
        CALL write_cnf_atoms ( 'initialize.cnf', n, box, box*r )
     END IF
  END IF
END PROGRAM initialize
