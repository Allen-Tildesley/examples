! pbc_to.f90
SUBROUTINE pbc ( r )
  ! Periodic boundary conditions for truncated octahedron

  ! The box is centred at the origin. 
  ! The containing cube is of unit length.
  ! The axes pass through the centres of 
  ! the six square faces of the truncated octahedron.
  ! Ref: Adams DJ, CCP5 quarterly, 10, 30 (1983)

  IMPLICIT NONE
  REAL, DIMENSION(3), INTENT(inout) :: r ! Argument

  REAL               :: corr
  REAL, PARAMETER    :: r75 = 4.0 / 3.0

  r(:) = r(:) - ANINT ( r(:) )
  corr = 0.5 * AINT ( r75 * SUM ( ABS ( r(:) ) ) )
  r(:) = r(:) - SIGN ( corr, r(:) )

END SUBROUTINE pbc
