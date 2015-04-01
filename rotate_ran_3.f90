! rotate_ran_3.f90
! Random rotation of unit vector of linear molecule
SUBROUTINE rotate_ran_3 ( dot_min, e )
  IMPLICIT NONE
  REAL,               INTENT(in)    :: dot_min ! determines max rotation angle
  REAL, DIMENSION(3), INTENT(inout) :: e       ! unit vector

  ! Ref: Marsaglia, Ann Maths Stat 43, 645 (1972)
  ! Uses a rejection technique to create a trial orientation
  ! subject to the constraint that the cosine of the angle
  ! turned through is greater than ( 1.0 - dot_min )

  REAL, DIMENSION(3) :: e0
  REAL, DIMENSION(2) :: zeta
  REAL               :: zeta_sq, f

  e0 = e ! store old unit vector

  DO

     DO
        CALL random_NUMBER ( zeta ) ! Two uniform random numbers between 0 and 1
        zeta = 2.0 * zeta - 1.0     ! now each between -1 and 1
        zeta_sq = SUM ( zeta**2 )   ! squared magnitude
        IF ( zeta_sq < 1.0 ) EXIT   ! now inside unit disk
     END DO

     f = 2.0 * SQRT ( 1.0 - zeta_sq )
     e = [ zeta(1) * f, zeta(2) * f, 1.0 - 2.0 * zeta_sq ] ! on surface of unit sphere
     IF ( dot_min + dot_PRODUCT ( e, e0 ) > 1.0 ) EXIT     ! close enough to e0

  END DO

END SUBROUTINE rotate_ran_3



