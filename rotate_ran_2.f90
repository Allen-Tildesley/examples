! rotate_ran_2.f90
! Random rotation of unit vector of linear molecule
SUBROUTINE rotate_ran_2 ( delta_max, e )
  USE utility_module, only : random_integer
  IMPLICIT NONE
  REAL,               INTENT(in)    :: delta_max ! max rotation angle
  REAL, DIMENSION(3), INTENT(inout) :: e         ! unit vector

  INTEGER            :: axis        ! axis of rotation
  REAL               :: delta, c, s ! rotation angle
  REAL, DIMENSION    :: zeta        ! random number
  REAL, DIMENSION(3) :: e0

  ! Ref: Barker and Watts, Chem Phys Lett 3, 144 (1969)
  ! Selects a Cartesian axis at random

  ! We use the built-in Fortran random number generator for simplicity
  CALL RANDOM_NUMBER ( zeta )              ! uniform random number between 0 and 1
  axis  = random_integer (1,3)             ! random axis choice 1 = x, 2 = y, 3 = z
  delta = ( 2.0 * zeta - 1.0 ) * delta_max ! uniform random angle
  c     = COS(delta)
  s     = SIN(delta)

  e0 = e

  SELECT CASE ( axis )

  CASE ( 1 ) ! rotate about x
     e(2) = c*e0(2) + s*e0(3)
     e(3) = c*e0(3) - s*e0(2)

  CASE ( 2 ) ! rotate about y
     e(3) = c*e0(3) + s*e0(1)
     e(1) = c*e0(1) - s*e0(3)

  CASE default ! rotate about z
     e(1) = c*e0(1) + s*e0(2)
     e(2) = c*e0(2) - s*e0(1)

  END SELECT

END SUBROUTINE rotate_ran_2
