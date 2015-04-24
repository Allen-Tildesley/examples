PROGRAM hit_and_miss
  IMPLICIT NONE
  REAL                          :: volume
  REAL, DIMENSION(3)            :: r_tau, zeta
  REAL, DIMENSION(3), PARAMETER :: r_0 = [1.0, 2.0, 3.0]
  REAL,               PARAMETER :: vol_0 = PRODUCT(r_0)
  INTEGER	                      :: tau, tau_shot, tau_hit

  CALL random_SEED()
  tau_hit = 0
  tau_shot = 1000000

  DO tau = 1, tau_shot
     CALL random_NUMBER ( zeta(:) ) ! uniform in range (0,1)
     r_tau = zeta * r_0             ! uniform in vol_0
     IF ( r_tau(2) < ( 2.0 - 2.0*r_tau(1) ) .AND. &
          &    r_tau(3) < ( 1.0 + r_tau(2) )           ) THEN
        tau_hit = tau_hit + 1
     END IF
  END DO
  volume = vol_0 * REAL ( tau_hit ) / REAL ( tau_shot ) 
  WRITE(*,'(''Estimate of volume = '',f10.5)') volume
END PROGRAM hit_and_miss
