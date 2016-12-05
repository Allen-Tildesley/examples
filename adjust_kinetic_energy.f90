! adjust_kinetic_energy.f90
! Utility program to allow user to change total energy of MD configuration
PROGRAM adjust_kinetic_energy

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Scales the velocities by an amount sufficient to shift the energy per particle by delta

  REAL, DIMENSION(:,:), ALLOCATABLE :: r, v

  REAL    :: kin_old, kin_new, delta, dummy
  INTEGER :: n

  CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp'

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( filename, n, dummy ) ! First call is just to get n and bond/box length
  ALLOCATE ( r(3,n), v(3,n) )
  CALL read_cnf_atoms ( filename, n, dummy, r, v ) ! Second call gets r and v

  ! Read desired change in energy per atom
  READ ( unit=input_unit, fmt=* ) delta

  kin_old = 0.5*SUM(v**2)           ! Kinetic energy
  kin_new = kin_old + delta*REAL(n)
  IF ( kin_new < 0.0 ) STOP 'Impossible energy change'
  v = v * SQRT(kin_new/kin_old)

  CALL write_cnf_atoms ( filename, n, dummy, r, v ) ! Write out final configuration

  DEALLOCATE ( r, v )
  
END PROGRAM adjust_kinetic_energy

