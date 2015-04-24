PROGRAM sample_mean
  IMPLICIT NONE
  REAL                          :: volume, sum
  REAL, DIMENSION(2)            :: r_tau, zeta
  REAL, DIMENSION(2), PARAMETER :: r_0 = [1.0, 2.0]
  INTEGER	                :: tau, tau_max

  CALL random_SEED()
  tau_max = 1000000

  sum = 0.0
  DO tau = 1, tau_max
     CALL random_NUMBER ( zeta ) ! uniform in (0,1)
     r_tau = zeta * r_0 ! uniform in xy rectangle
     IF ( r_tau(2) < ( 2.0 - 2.0*r_tau(1) ) ) THEN ! within xy triangle
        sum = sum + ( 1.0 + r_tau(2) ) ! value of z
     END IF
  END DO
  volume = 2.0 * sum / REAL ( tau_max ) 
  WRITE(*,'(''Estimate of volume = '',f10.5)') volume
END PROGRAM sample_mean
