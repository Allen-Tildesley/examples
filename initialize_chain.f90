! initialize_chain.f90
! Sets up initial chain configuration for MD or MC
PROGRAM initialize_chain
  USE utility_module,    ONLY : write_cnf_atoms
  USE initialize_module, ONLY : initialize_chain_lattice, initialize_chain_random, initialize_chain_velocities, n, r, v
  IMPLICIT NONE

  INTEGER :: nc
  REAL    :: temperature, bond
  LOGICAL :: velocities, lattice

  NAMELIST /initialize_parameters/ nc, n, temperature, velocities, lattice, bond

  ! Default values
  ! We must specify either n or nc but not both
  n           = -1
  nc          = -1
  temperature = 1.0
  bond        = 1.0
  velocities  = .FALSE.
  lattice     = .true.

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
        WRITE(*,'(''nc = '',t40,i5)') nc
     ELSE
        STOP 'Must specify either n or nc, not both'
     END IF
  END IF

  CALL RANDOM_SEED()

  IF ( lattice ) THEN
  ! The chain takes the form of a close-packed cube of atoms surrounded by vacuum
     CALL initialize_chain_lattice ! coordinates initialized for unit bond length
  ELSE
     ! The chain bonds are randomly oriented, avoiding overlaps
     CALL initialize_chain_random ! coordinates initialized for unit bond length
  END IF
  r = r * bond ! positions scaled to give desired bond length
  
  IF ( velocities ) THEN
     CALL initialize_chain_velocities ( temperature )
  END IF

  ! We do not use periodic boundaries for this system
  ! Instead, use box variable to store bond length

  ! Write out coordinates
  IF ( velocities ) THEN
     CALL write_cnf_atoms ( 'md_chain.cnf', n, bond, r, v )
  ELSE
     CALL write_cnf_atoms ( 'mc_chain.cnf', n, bond, r )
  END IF
END PROGRAM initialize_chain
