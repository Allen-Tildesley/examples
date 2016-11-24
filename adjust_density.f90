! adjust_density.f90
! Utility program to allow user to change the density of MC or MD configuration
PROGRAM adjust_energy

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions and optionally, velocities)
  ! Cubic periodic boundary conditions

  ! Input configuration and output configuration are assumed to be in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1
  ! There is nothing here specific to Lennard-Jones

  ! Most important variables
  REAL :: box     ! Box length
  REAL :: density ! Desired density
  REAL :: scale   ! Scaling factor

  LOGICAL :: velocities
  INTEGER :: n

  REAL, DIMENSION(:,:), ALLOCATABLE :: r, v ! Positions and velocities (3,n)

  CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp'

  WRITE ( unit=output_unit, fmt='(a)'   ) 'adjust_density'
  WRITE ( unit=output_unit, fmt='(a)'   ) 'Allows user to fix density of MC or MD configuration'
  WRITE ( unit=output_unit, fmt='(a,a)' ) 'Input and output configuration in ', filename

  WRITE ( unit=output_unit, fmt='(a)', advance='no' ) 'File includes velocities? (.t./.f.) > '
  READ ( unit=input_unit, fmt=* ) velocities

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( filename, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n) / box**3
  ALLOCATE ( r(3,n), v(3,n) )
  IF ( velocities ) THEN
     CALL read_cnf_atoms ( filename, n, box, r, v ) ! Second call gets r and v
  ELSE
     CALL read_cnf_atoms ( filename, n, box, r ) ! Second call gets r
  END IF

  ! Read desired density and scale coordinates
  WRITE ( unit=output_unit, fmt='(a)', advance='no' ) 'Enter desired density > '
  READ ( unit=input_unit, fmt=* ) density
  scale  = ( (REAL(n)/density)**(1.0/3.0) ) / box
  box    = box * scale
  r(:,:) = r(:,:) * scale

  IF ( velocities ) THEN
     CALL write_cnf_atoms ( filename, n, box, r, v ) ! Write out final configuration
  ELSE
     CALL write_cnf_atoms ( filename, n, box, r ) ! Write out final configuration
  END IF

  DEALLOCATE ( r, v )

END PROGRAM adjust_energy

