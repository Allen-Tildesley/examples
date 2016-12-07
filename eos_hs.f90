! eos_hs.f90
! Equation of State for hard sphere potential
PROGRAM eos_hs

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  ! This program uses the function derived by H Hansen-Goos, J Chem Phys, 16, 164506 (2016)
  ! which is claimed to be an excellent fit to simulation data over the whole fluid density range
  ! That paper also gives references to previous approximate equations of state (such as the
  ! venerable Carnhahan-Starling equation).

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
  WRITE ( unit=output_unit, fmt='(a,f15.6)' ) 'rho = ', density
  WRITE ( unit=output_unit, fmt='(a,f15.6)' ) 'eta = ', eta

  ! Equation (6) of Hansen-Goos (2016)
  z = a * LOG ( 1.0-eta ) / eta
  z = z + polynomial ( b, eta ) / ( 1.0 - eta ) ** 3 ! Compressibility factor P/(rho*kT)
  p = z * density                                    ! Pressure P / kT
  WRITE ( unit=output_unit, fmt='(a,t20,f15.6)' ) 'P', p
  WRITE ( unit=output_unit, fmt='(a,t20,f15.6)' ) 'Z = P/(rho*kT)', z

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
