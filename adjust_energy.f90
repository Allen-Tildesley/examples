! adjust_energy.f90
! Utility program to allow user to fix total energy of MD configuration
PROGRAM adjust_energy

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, r, v, n, potential_type

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions

  ! Positions r are divided by box length after reading in and we assume mass=1 throughout
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! There is nothing here specific to Lennard-Jones; the model is defined in md_module

  ! Most important variables
  REAL            :: box         ! Box length
  REAL, PARAMETER :: r_cut = 2.5 ! Potential cutoff distance (adjust this if necessary)

  ! Composite interaction = pot & cut & vir & lap & ovr variables
  TYPE(potential_type) :: total

  REAL, DIMENSION(3) :: vcm
  REAL               :: kin_old, kin_new, es_old, es_new

  CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp'

  WRITE ( unit=output_unit, fmt='(a)'   ) 'adjust_energy'
  WRITE ( unit=output_unit, fmt='(a)'   ) 'Allows user to fix total energy of MD configuration'
  WRITE ( unit=output_unit, fmt='(a)'   ) 'Particle mass=1 throughout'
  WRITE ( unit=output_unit, fmt='(a,a)' ) 'Input and output configuration in ', filename
  CALL introduction

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( filename, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n) / box**3
  CALL allocate_arrays ( box, r_cut )
  CALL read_cnf_atoms ( filename, n, box, r, v )            ! Second call gets r and v
  r(:,:) = r(:,:) / box                                     ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) )                        ! Periodic boundaries
  vcm(:) = SUM ( v(:,:), dim=2 ) / REAL(n)                  ! Centre-of mass velocity
  v(:,:) = v(:,:) - SPREAD ( vcm(:), dim = 2, ncopies = n ) ! Set COM velocity to zero

  ! Initial forces, potential, etc plus overlap check
  CALL force ( box, r_cut, total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in adjust_energy'
  END IF

  kin_old = 0.5*SUM(v**2)           ! Kinetic energy
  es_old = (kin_old+total%pot)/REAL(n)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'E/N cut-and-shifted', es_old
  WRITE ( unit=output_unit, fmt='(a)'           ) 'Enter desired E/N'
  READ  ( unit=input_unit,  fmt=*               ) es_new
  kin_new = es_new*REAL(n) - total%pot
  IF ( kin_new < 0.0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,t40,f15.6)' ) 'Would require negative KE', kin_new
     STOP 'Error in adjust_energy'
  END IF
  v = v * SQRT(kin_new/kin_old)

  CALL write_cnf_atoms ( filename, n, box, r*box, v ) ! Write out final configuration

  CALL deallocate_arrays

END PROGRAM adjust_energy

