! eos_lj.f90
! Equation of State for Lennard-Jones pair potential
PROGRAM eos_lj

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             !
  ! to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     !
  ! published by Oxford University Press ("the publishers").                                       !
  !                                                                                                !
  ! LICENCE                                                                                        !
  ! Creative Commons CC0 Public Domain Dedication.                                                 !
  ! To the extent possible under law, the authors have dedicated all copyright and related         !
  ! and neighboring rights to this software to the PUBLIC domain worldwide.                        !
  ! This software is distributed without any warranty.                                             !
  ! You should have received a copy of the CC0 Public Domain Dedication along with this software.  !
  ! If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               !
  !                                                                                                !
  ! DISCLAIMER                                                                                     !
  ! The authors and publishers make no warranties about the software, and disclaim liability       !
  ! for all uses of the software, to the fullest extent permitted by applicable law.               !
  ! The authors and publishers do not recommend use of this software for any purpose.              !
  ! It is made freely available, solely to clarify points made in the text. When using or citing   !
  ! the software, you should not imply endorsement by the authors or publishers.                   !
  !------------------------------------------------------------------------------------------------!

  ! The routines in the above module use the fitting function described and parametrized in
  ! M Thol, G Rutkai, R Span, J Vrabec, R Lustig, Int J Thermophys 36, 25 (2015)
  ! M Thol, G Rutkai, A Koester, R Lustig, R Span, J Vrabec, J Phys Chem Ref Data 45, 023101 (2016)
  ! Those authors also supply C++ codes (in the supplementary information of those papers)
  ! They are NOT responsible for this Fortran code, which was written independently by Michael P Allen
  ! A similar notation, consistent with the papers, is retained for clarity.

  ! Formulae for P, E/N etc in terms of the scaled free energy derivatives a_res(0,1) etc
  ! may be found in the above papers

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               lrc_module,      ONLY : potential_lrc, pressure_lrc, pressure_delta
  USE               eos_lj_module,   ONLY : a_res_full, a_res_cutshift

  IMPLICIT NONE

  REAL    :: temperature, density, e, p, cv, cp, mu, z
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

  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature T', temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density rho',   density

  ! Results for full potential from Thol et al (2016) fitting formula
  
  WRITE ( unit=output_unit, fmt='(/,a,/)' ) 'Full Lennard-Jones potential'

  a_res = a_res_full ( temperature, density )

  DO i = 0, 2
     DO j = 0, 2
        IF ( i+j > 2 ) CYCLE ! Only interested in some of the results
        WRITE ( unit=output_unit, fmt='(a4,2i1,t40,f15.6)' ) 'Ares', i, j, a_res(i,j)
     END DO
  END DO

  p  = density * temperature * ( 1.0 + a_res(0,1) )
  e  = temperature * ( 1.5 + a_res(1,0) )
  cv = 1.5 - a_res(2,0)
  cp = 2.5 - a_res(2,0)+(1.0+a_res(0,1)-a_res(1,1))*(1.0+a_res(0,1)-a_res(1,1))/(1.0+2.0*a_res(0,1)+a_res(0,2)) - 1.0
  mu = temperature * ( LOG(density) + a_res(0,0) + a_res(0,1) )
  z  = density * EXP ( a_res(0,0) + a_res(0,1) )
  WRITE ( unit=output_unit, fmt='(/,a,t40,f15.6)' ) 'Pressure P',            p
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Energy E/N',            e
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Heat capacity Cv/NkB',  cv
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Heat capacity Cp/NkB',  cp
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Chemical potential mu', mu
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Activity z',            z

  ! Estimates for cut (but not shifted) potential by reverse-application of long-range & delta corrections
  
  WRITE ( unit=output_unit, fmt='(/,a,/)' ) 'Lennard-Jones potential cut (but not shifted) at 2.5 sigma'
  p  = p - pressure_lrc ( density, r_cut ) + pressure_delta ( density, r_cut )
  e  = e - potential_lrc ( density, r_cut )
  mu = mu - 2.0 * potential_lrc ( density, r_cut )
  z  = z * EXP ( -2.0* potential_lrc ( density, r_cut ) / temperature )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Pressure P',            p
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Energy E/N',            e
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Chemical potential mu', mu
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Activity z',            z

  ! Results for cut-and-shifted potential from Thol et al (2015) fitting formula
  
  WRITE ( unit=output_unit, fmt='(/,a,/)' ) 'Lennard-Jones potential cut-and-shifted at 2.5 sigma'

  a_res = a_res_cutshift ( temperature, density )

  DO i = 0, 2
     DO j = 0, 2
        IF ( i+j > 2 ) CYCLE ! Only interested in some of the results
        WRITE ( unit=output_unit, fmt='(a4,2i1,t40,f15.6)' ) 'Ares', i, j, a_res(i,j)
     END DO
  END DO

  p  = density * temperature * ( 1.0 + a_res(0,1) )
  e  = temperature * ( 1.5 + a_res(1,0) )
  cv = 1.5 - a_res(2,0)
  cp = 2.5 - a_res(2,0)+(1.0+a_res(0,1)-a_res(1,1))*(1.0+a_res(0,1)-a_res(1,1))/(1.0+2.0*a_res(0,1)+a_res(0,2)) - 1.0
  mu = temperature * ( LOG(density) + a_res(0,0) + a_res(0,1) )
  z  = density * EXP ( a_res(0,0) + a_res(0,1) )
  WRITE ( unit=output_unit, fmt='(/,a,t40,f15.6)' ) 'Pressure P',            p
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Energy E/N',            e
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Heat capacity Cv/NkB',  cv
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Heat capacity Cp/NkB',  cp
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Chemical potential mu', mu
  WRITE ( unit=output_unit, fmt='(  a,t40,f15.6)' ) 'Activity z',            z

END PROGRAM eos_lj
