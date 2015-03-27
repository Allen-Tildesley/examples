! pbc_rh.f90
SUBROUTINE pbc ( r )
  ! Periodic boundary conditions for rhombus

  ! Periodic corrections are applied in two dimensions x, y.
  ! In most applications the molecules will be confined in the
  ! z direction by real walls rather than by periodic boundaries.
  ! The box is centred at the origin. 
  ! The x axis lies along the side of the rhombus, which is of unit length
  ! The acute angle of the rhombus is sixty degrees
  ! Ref: Talbot J, private communication (1987)

  IMPLICIT NONE
  REAL, DIMENSION(3), INTENT(inout) :: r ! Argument

  REAL, PARAMETER :: rt3 = SQRT(3.0), rrt3 = 1.0 / rt3
  REAL, PARAMETER :: rt32 = rt3 / 2.0, rrt32 = 1.0 / rt32

  r(1) = r(1) - ANINT ( r(1) - rrt3 * r(2) ) - ANINT ( rrt32 * r(2) ) * 0.5
  r(2) = r(2) - ANINT ( rrt32 * r(2) ) * rt32

END SUBROUTINE pbc



