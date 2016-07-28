PROGRAM sample_mean
  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit
  IMPLICIT NONE
  REAL                          :: volume, sum
  REAL, DIMENSION(2)            :: r_tau, zeta
  REAL, DIMENSION(2), PARAMETER :: r_0 = [1.0, 2.0]
  INTEGER	                :: tau, tau_max

  CALL RANDOM_SEED()
  tau_max = 1000000

  sum = 0.0
  DO tau = 1, tau_max
     CALL RANDOM_NUMBER ( zeta ) ! uniform in (0,1)
     r_tau = zeta * r_0 ! uniform in xy rectangle
     IF ( r_tau(2) < ( 2.0 - 2.0*r_tau(1) ) ) THEN ! within xy triangle
        sum = sum + ( 1.0 + r_tau(2) ) ! value of z
     END IF
  END DO
  volume = 2.0 * sum / REAL ( tau_max ) 
  WRITE ( unit=output_unit, fmt='(a,f10.5)' ) 'Estimate of volume = ', volume
END PROGRAM sample_mean
