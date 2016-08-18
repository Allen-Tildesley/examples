! hit_and_miss.f90
! Estimates volume of polyhedron by simple MC
PROGRAM hit_and_miss

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit

  IMPLICIT NONE
  REAL                          :: v
  REAL, DIMENSION(3)            :: r, zeta
  REAL, DIMENSION(3), PARAMETER :: r_0 = [1.0, 2.0, 3.0]
  REAL,               PARAMETER :: v_0 = PRODUCT(r_0)
  INTEGER                       :: tau, tau_shot, tau_hit

  CALL RANDOM_SEED()
  tau_hit  = 0
  tau_shot = 1000000

  DO tau = 1, tau_shot
     CALL RANDOM_NUMBER ( zeta(:) ) ! uniform in range (0,1)
     r = zeta * r_0                 ! uniform in v_0
     IF (   r(2) < ( 2.0 - 2.0*r(1) ) .AND. &
          & r(3) < ( 1.0 + r(2) )   ) THEN ! in polyhedron
        tau_hit = tau_hit + 1
     END IF
  END DO
  v = v_0 * REAL ( tau_hit ) / REAL ( tau_shot ) 
  WRITE ( unit=output_unit, fmt='(a,f10.5)' ) 'Estimate = ', v
END PROGRAM hit_and_miss
