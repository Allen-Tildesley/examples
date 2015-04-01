! rotate_ran_1.f90
! Random change in polar angles of linear molecule
SUBROUTINE rotate_ran_1 ( dphi_max, dcos_max, e  )
  IMPLICIT NONE
  REAL,               INTENT(in)    :: dphi_max ! max change in phi
  REAL,               INTENT(in)    :: dcos_max ! max change in cos(theta)
  REAL, DIMENSION(3), INTENT(inout) :: e        ! unit vector

  REAL, PARAMETER    :: twopi = 8.0 * ATAN(1.0)
  REAL               :: cos_theta, sin_theta, phi ! angles
  REAL, DIMENSION(2) :: zeta ! random numbers

  ! In this example, the angles phi and theta are calculated from the
  ! unit vector e, and vice versa, but it is also possible to simply
  ! provide the angles directly as arguments.
  
  ! There is a slight cheat: if the new cos(theta) is outside (-1,+1)
  ! we bring it inside by wrapping around the boundary. The resulting
  ! distribution is uniform, but the rotation is large (+/-pi) not small!

  ! We use the built-in Fortran random number generator for simplicity
  CALL RANDOM_NUMBER ( zeta ) ! two uniform random numbers between 0 and 1
  zeta = 2.0*zeta - 1.0       ! now each between -1 and +1

  ! convert unit vector to angles
  cos_theta = e(3)
  phi       = ATAN2 ( e(2), e(1) )

  phi = phi + zeta(1) * dphi_max
  phi = phi - ANINT ( phi / twopi ) * twopi ! bring within range
  cos_theta = cos_theta + zeta(2) * dcos_max
  cos_theta = cos_theta - ANINT ( cos_theta / 2.0 ) * 2.0 ! bring within range
  sin_theta = SQRT ( 1.0 - cos_theta ** 2 )

  ! convert euler angles to unit vector
  e = [ COS ( phi ) * sin_theta, SIN ( phi ) * sin_theta, cos_theta ]
END SUBROUTINE rotate_ran_1
