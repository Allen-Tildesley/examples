  ! initialize.f90 (uses initialize_module.f90, utility_module.f90)
  ! Sets up initial configuration for MD or MC
PROGRAM initialize
  USE utility_module, ONLY : write_cnf_atoms, write_cnf_molecules
  USE initialize_module, ONLY : n, r, e, v, w
  IMPLICIT NONE

  NAMELIST /initialize_parameters/ nc, n

  ! We must specify either n or nc but not both
  n = -1
  nc = -1

  READ(*,nml=initialize_parameters)

  IF ( n < 0 ) THEN
     IF ( nc < 0 ) THEN
        stop 'Must specify either n or nc'
     ELSE
        n = 4*nc**3 ! fcc lattice
     END IF
  ELSE
     IF ( nc < 0 ) THEN
        nc = NINT ( ( REAL(n)/4.0 )**(1.0/3.0) )
        IF ( n /= 4*nc**3 ) STOP 'n not of the form  4*nc**3'
     ELSE
        STOP 'Must specify either n or nc, not both'
     END IF
  END IF

  CALL random_SEED()


  CALL init_fcc ( nc )
  call init_vel ( temperature )

END PROGRAM initialize
