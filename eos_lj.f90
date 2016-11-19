! eos_lj.f90
! Equation of State for Lennard-Jones pair potential
PROGRAM eos_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE mc_module,     ONLY : potential_lrc, pressure_lrc, pressure_delta
  USE eos_lj_module, ONLY : a_res_full, a_res_cutshift

  ! The routines in the above module use the fitting function described and parametrized in
  ! M Thol, G Rutkai, R Span, J Vrabec, R Lustig, Int J Thermophys 36, 25 (2015)
  ! M Thol, G Rutkai, A Koester, R Lustig, R Span, J Vrabec, J Phys Chem Ref Data 45, 023101 (2016)
  ! Those authors also supply C++ codes (in the supplementary information of those papers)
  ! They are NOT responsible for this Fortran code, which was written independently by Michael P Allen
  ! A similar notation, consistent with the papers, is retained for clarity.

  ! Formulae for P, E/N etc in terms of the scaled free energy derivatives a_res(0,1) etc
  ! may be found in the abbove papers

  IMPLICIT NONE

  REAL    :: temperature, density, e, p, cv, cp
  INTEGER :: ioerr, i, j

  REAL, DIMENSION(0:2,0:2) :: a_res ! Residual free energy and scaled derivatives

  REAL, PARAMETER :: r_cut = 2.5 ! The cut-and-shifted subroutine applies only to this value

  NAMELIST /nml/ temperature, density

  ! Set sensible default values
  temperature = 1.0
  density     = 0.75

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in eos_lj'
  END IF

  WRITE ( unit=output_unit, fmt='(a,f15.5)' ) 'T   = ', temperature
  WRITE ( unit=output_unit, fmt='(a,f15.5)' ) 'rho = ', density

  WRITE ( unit=output_unit, fmt='(a)' ) 'Full Lennard-Jones potential'

  a_res = a_res_full ( temperature, density )

  DO i = 0, 2
     DO j = 0, 2
        IF ( i+j > 2 ) CYCLE ! Only interested in some of the results
        WRITE ( unit=output_unit, fmt='(a4,2i1,2x,f15.8)' ) 'Ares', i, j, a_res(i,j)
     END DO
  END DO

  p = density * temperature * ( 1.0 + a_res(0,1) )
  e = temperature * ( 1.5 + a_res(1,0) )
  cv = 1.5 - a_res(2,0)
  cp = 2.5 - a_res(2,0)+(1.0+a_res(0,1)-a_res(1,1))*(1.0+a_res(0,1)-a_res(1,1))/(1.0+2.0*a_res(0,1)+a_res(0,2)) - 1.0
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'P', p
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'E/N', e
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'Cv/NkB', cv
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'Cp/NkB', cp

  WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential cut (but not shifted) at 2.5 sigma'
  p  = p - pressure_lrc ( density, r_cut ) + pressure_delta ( density, r_cut )
  e  = e - potential_lrc ( density, r_cut )
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'P', p
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'E/N', e

  WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential cut-and-shifted at 2.5 sigma'

  a_res = a_res_cutshift ( temperature, density )

  DO i = 0, 2
     DO j = 0, 2
        IF ( i+j > 2 ) CYCLE ! Only interested in some of the results
        WRITE ( unit=output_unit, fmt='(a4,2i1,2x,f15.8)' ) 'Ares', i, j, a_res(i,j)
     END DO
  END DO

  p = density * temperature * ( 1.0 + a_res(0,1) )
  e = temperature * ( 1.5 + a_res(1,0) )
  cv = 1.5 - a_res(2,0)
  cp = 2.5 - a_res(2,0)+(1.0+a_res(0,1)-a_res(1,1))*(1.0+a_res(0,1)-a_res(1,1))/(1.0+2.0*a_res(0,1)+a_res(0,2)) - 1.0
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'P', p
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'E/N', e
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'Cv/NkB', cv
  WRITE ( unit=output_unit, fmt='(a,t20,f15.5)' ) 'Cp/NkB', cp

END PROGRAM eos_lj
