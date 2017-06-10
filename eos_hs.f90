! eos_hs.f90
! Equation of State for hard sphere potential
PROGRAM eos_hs

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

  ! This program uses the function derived by H Hansen-Goos, J Chem Phys, 16, 164506 (2016)
  ! which is claimed to be an excellent fit to simulation data over the whole fluid density range
  ! That paper also gives references to previous approximate equations of state (such as the
  ! venerable Carnahan-Starling equation).

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  IMPLICIT NONE

  REAL            :: density, p, z, eta
  INTEGER         :: ioerr
  REAL, PARAMETER :: pi = 4.0*ATAN(1.0)

  ! The coefficients appear in Table I of Hansen-Goos (2016)
  REAL,                 PARAMETER :: a = 8.0
  REAL, DIMENSION(0:7), PARAMETER :: b = [ 9.0,     -19.0,      47.0/3.0, -2.635232, &
       &                                  -1.265575,  0.041212, 0.248245, -0.096495 ]

  NAMELIST /nml/ density

  ! Set sensible default values
  density = 0.75

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in eos_hs'
  END IF

  eta = pi * density / 6.0 ! Packing fraction
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density rho',          density
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Packing fraction eta', eta

  ! Equation (6) of Hansen-Goos (2016)
  z = a * LOG ( 1.0-eta ) / eta
  z = z + polynomial ( b, eta ) / ( 1.0 - eta ) ** 3 ! Compressibility factor P/(rho*kT)
  p = z * density                                    ! Pressure P / kT
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Pressure P',                            p
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Compressibility factor Z = P/(rho*kT)', z

CONTAINS

  FUNCTION polynomial ( c, x ) RESULT ( f )
    REAL, DIMENSION(:), INTENT(in) :: c ! coefficients of x**0, x**1, x**2 etc
    REAL,               INTENT(in) :: x ! argument
    REAL                           :: f ! result

    INTEGER :: i

    f = 0.0
    DO i = SIZE(c), 1, -1
       f = f * x + c(i)
    END DO

  END FUNCTION polynomial

END PROGRAM eos_hs
