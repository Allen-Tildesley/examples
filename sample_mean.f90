PROGRAM sample_mean
  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit

  IMPLICIT NONE

  REAL                          :: v, z_sum
  REAL, DIMENSION(2)            :: r, zeta
  REAL, DIMENSION(2), PARAMETER :: r_0 = [1.0, 2.0]
  REAL,               PARAMETER :: a_0 = PRODUCT(r_0)
  INTEGER                       :: tau, tau_max

  CALL RANDOM_SEED()
  tau_max = 1000000

  z_sum = 0.0
  DO tau = 1, tau_max
     CALL RANDOM_NUMBER ( zeta ) ! uniform in (0,1)
     r = zeta * r_0 ! uniform in xy rectangle
     IF ( r(2) < 2.0-2.0*r(1) ) THEN   ! in xy triangle
        z_sum = z_sum + ( 1.0 + r(2) ) ! value of z
     END IF
  END DO
  v = a_0 * z_sum / REAL ( tau_max ) 
  WRITE ( unit=output_unit, fmt='(a,f10.5)' ) 'Estimate = ', v
END PROGRAM sample_mean
