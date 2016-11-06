! eos_lj.f90
! Equation of State for Lennard-Jones pair potential
PROGRAM eos_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE mc_module, ONLY : potential_lrc, pressure_lrc, pressure_delta

  USE eos_lj_module, ONLY : plj, ulj

  IMPLICIT NONE

  REAL    :: temperature, density, r_cut, u_full, u_cut, p_full, p_cut
  INTEGER :: ioerr

  NAMELIST /nml/ temperature, density, r_cut

  ! Set sensible default values
  temperature = 1.0
  density     = 0.75
  r_cut       = 2.5

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in eos_lj'
  END IF

  WRITE ( unit=output_unit, fmt='(a,f15.5)' ) 'T   = ', temperature
  WRITE ( unit=output_unit, fmt='(a,f15.5)' ) 'rho = ', density
  p_full = plj ( temperature, density )
  p_cut  = p_full - pressure_lrc ( density, r_cut ) + pressure_delta ( density, r_cut )
  u_full = ulj ( temperature, density )
  u_cut  = u_full - potential_lrc ( density, r_cut )
  WRITE ( unit=output_unit, fmt='(a)'       ) 'Full potential, including ideal gas part'
  WRITE ( unit=output_unit, fmt='(a,f15.5)' ) 'P   = ', p_full
  WRITE ( unit=output_unit, fmt='(a,f15.5)' ) 'E/N = ', u_full
  WRITE ( unit=output_unit, fmt='(a)'       ) 'Cut (but not shifted) potential, including ideal gas part'
  WRITE ( unit=output_unit, fmt='(a,f15.5)' ) 'P   = ', p_cut
  WRITE ( unit=output_unit, fmt='(a,f15.5)' ) 'E/N = ', u_cut

END PROGRAM eos_lj
