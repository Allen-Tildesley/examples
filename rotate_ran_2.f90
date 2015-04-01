! rotate_ran_2.f90
! Random rotation of unit vector of linear molecule
SUBROUTINE rotate_ran_2 ( delta_max, e )
  IMPLICIT NONE
  REAL,               INTENT(in)    :: delta_max ! max rotation angle
  REAL, DIMENSION(3), INTENT(inout) :: e         ! unit vector

  INTEGER            :: axis  ! axis of rotation
  REAL               :: delta ! rotation angle
  REAL, DIMENSION(2) :: zeta  ! random numbers
  REAL, DIMENSION(3) :: e0

  ! Ref: Barker and Watts, Chem Phys Lett 3, 144 (1969)
  ! Selects a Cartesian axis at random

  ! We use the built-in Fortran random number generator for simplicity
  CALL RANDOM_NUMBER ( zeta )                 ! two uniform random numbers between 0 and 1
  axis  = FLOOR ( 3.0 * zeta (1) ) + 1        ! random axis choice 1 = x, 2 = y, 3 = z
  delta = ( 2.0 * zeta(2) - 1.0 ) * delta_max ! uniform random angle

  e0 = e

  SELECT CASE ( axis )

  CASE ( 1 ) ! rotate about x
     e(2) = COS(delta)*e0(2) + SIN(delta)*e0(3)
     e(3) = COS(delta)*e0(3) - SIN(delta)*e0(2)

  CASE ( 2 ) ! rotate about y
     e(3) = COS(delta)*e0(3) + SIN(delta)*e0(1)
     e(1) = COS(delta)*e0(1) - SIN(delta)*e0(3)

  CASE default ! rotate about z
     e(1) = COS(delta)*e0(1) + SIN(delta)*e0(2)
     e(2) = COS(delta)*e0(2) - SIN(delta)*e0(1)

  END SELECT

END SUBROUTINE rotate_ran_2
