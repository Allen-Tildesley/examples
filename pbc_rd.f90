! pbc_rd.f90
SUBROUTINE pbc ( r )
  ! Periodic boundary conditions for rhombic dodecahedron

  ! The box is centred at the origin. 
  ! The x and y axes join the centres of opposite faces of the dodecahedron.
  ! The z axis joins opposite vertices of the dodecahedron.
  ! The diagonal of the rhombic face is of unit length.
  ! The side of the containing cube is sqrt(2.0).
  ! The x and y axes pass through the cube edges.
  ! The z axis passes through the cube faces.
  ! Ref: Smith W, CCP5 quarterly, 10, 37 (1983).

  IMPLICIT NONE
  REAL, DIMENSION(3), INTENT(inout) :: r ! Argument

  REAL            :: corr
  REAL, PARAMETER :: rt2 = SQRT(2.0), rrt2 = 1.0 / rt2

  r(1) = r(1) - ANINT ( r(1) )
  r(2) = r(2) - ANINT ( r(2) )
  r(3) = r(3) - rt2 * ANINT ( rrt2 * r(3) )
  corr = 0.5 * AINT ( ABS ( r(1) ) + ABS ( r(2) ) + rt2 * ABS ( r(3) ) )
  r(1) = r(1) - SIGN ( corr, r(1) )
  r(2) = r(2) - SIGN ( corr, r(2) )
  r(3) = r(3) - SIGN ( corr, r(3) ) * rt2

END SUBROUTINE pbc
