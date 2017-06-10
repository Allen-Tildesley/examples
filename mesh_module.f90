! mesh_module.f90
! Provides function to convert charges to charge density on a 3-d mesh
MODULE mesh_module

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             !
  ! to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     !
  ! published by Oxford University Press ("the publishers").                                       !
  !                                                                                                !
  ! LICENCE                                                                                        !
  ! Creative Commons CC0 Public Domain Dedication.                                                 !
  ! To the extent possible under law, the authors have dedicated all copyright and related         !
  ! and neighboring rights to this software to the PUBLIC domain worldwide.                        !
  ! This software is distributed without any warranty.                                             !
  ! You should have received a copy of the CC0 Public Domain Dedication along with this software.  !
  ! If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               !
  !                                                                                                !
  ! DISCLAIMER                                                                                     !
  ! The authors and publishers make no warranties about the software, and disclaim liability       !
  ! for all uses of the software, to the fullest extent permitted by applicable law.               !
  ! The authors and publishers do not recommend use of this software for any purpose.              !
  ! It is made freely available, solely to clarify points made in the text. When using or citing   !
  ! the software, you should not imply endorsement by the authors or publishers.                   !
  !------------------------------------------------------------------------------------------------!

  ! The function mesh_function assigns a set of charges to a cubic mesh using the
  ! triangular-shaped cloud distribution described by Hockney and Eastwood (1988)
  ! It is assumed that the charges are in a box of unit length.
  ! The charge mesh is indexed from 0 to sc-1 in each coordinate direction

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: mesh_function, sharpen

CONTAINS

  FUNCTION mesh_function ( r, q, sc ) RESULT ( rho )
    IMPLICIT NONE
    REAL,    DIMENSION(:,:),      INTENT(in) :: r   ! Charge positions (3,n) 
    REAL,    DIMENSION(:),        INTENT(in) :: q   ! Charges (n)
    INTEGER,                      INTENT(in) :: sc  ! Dimension of mesh
    REAL,    DIMENSION(0:sc-1,0:sc-1,0:sc-1) :: rho ! Returns mesh cell charge density

    REAL,    DIMENSION(3)      :: dr ! Charge positions relative to mesh point
    REAL,    DIMENSION(3,-1:1) :: v  ! Weights in each coordinate direction
    INTEGER, DIMENSION(3)      :: nr ! Mesh point index

    INTEGER :: n, i, i1, i2, i3, n1, n2, n3
    REAL    :: h, q1, q2, q3

    n = SIZE ( q )
    IF ( ANY ( SHAPE(r) /= [3,n] ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,4i15)' ) 'r shape error', SHAPE(r), 3, n
       STOP 'Error in mesh_function'
    END IF

    h = 1.0 / REAL( sc ) ! Mesh spacing

    rho = 0.0 ! zero mesh

    DO i = 1, n ! Loop over charges

       nr(:) = NINT ( r(:,i) * sc )       ! Nearest mesh point indices
       nr(:) = MODULO ( nr, sc)           ! With periodic boundaries
       dr(:) = r(:,i) - REAL( nr(:) ) * h ! Vector to charge from mesh cell centre
       dr(:) = dr(:) - ANINT( dr(:) )     ! Periodic boundaries
       dr(:) = dr(:) / h                  ! Normalise by mesh cell size

       ! Weights for three point assignment scheme
       v(:,-1) = 0.5 * ( 0.5 - dr(:) ) ** 2
       v(:, 0) = 0.75 - dr(:)**2
       v(:,+1) = 0.5 * ( 0.5 + dr(:) ) ** 2

       DO i1 = -1, 1 ! Loop over x-neighbours
          q1 = q(i) * v(1,i1)            ! Charge times x-weight
          n1 = MODULO ( nr(1) + i1, sc ) ! Neighbour index with periodic boundaries

          DO i2 = -1, 1 ! Loop over y-neighbours
             q2 = q1 * v(2,i2)              ! Charge times xy-weight
             n2 = MODULO ( nr(2) + i2, sc ) ! Neighbour index with periodic boundaries

             DO i3 = -1, 1 ! Loop over z-neighbours
                q3 = q2 * v(3,i3)                  ! Charge times xyz-weight
                n3 = MODULO ( nr(3) + i3, sc )     ! Neighbour index with periodic boundaries
                rho(n1,n2,n3) = rho(n1,n2,n3) + q3 ! Mesh cell share of charge
             END DO ! End loop over z-neighbours

          END DO ! End loop over y-neighbours

       END DO ! End loop over x-neighbours

    END DO ! End loop over charges

    rho = rho / (h**3) ! Convert charges to charge densities

  END FUNCTION mesh_function

  FUNCTION sharpen ( x ) RESULT ( f )
    REAL             :: f ! Returns sharpening factor for particle-mesh Ewald influence function
    REAL, INTENT(in) :: x ! Argument

    ! In the particle-mesh Ewald method, the accuracy may be improved by optimizing the influence function G
    ! used to calculate energies and forces from the Fourier-transformed charge density.
    ! To illustrate this, we use a simple sharpening function, eqn (2.22) of
    ! V Ballenegger, JJ Cerda, C Holm, J Chem Theo Comp 8, 936 (2012)
    ! We put this routine here because it assumes that we use the P=3 triangular-shaped cloud distribution
    ! in assigning the charge distribution to a cubic mesh
    ! The argument x is pi*i/sc, where i indexes the k-vector,
    ! or equivalently 0.5*k*h where k is the k-value and h is the mesh spacing
    ! This function is applied to x, y, and z components separately

    REAL            :: u ! The U function (in one coordinate) of the text, and of Balleneger et al (2012)
    REAL, PARAMETER :: tol = 1.0e-3

    IF ( ABS(x) < tol ) THEN
       u = 1 - 0.5*x**2 ! Taylor series
    ELSE
       u = ( SIN(x)/x )**3
    END IF

    f = 1.0 / u**2

  END FUNCTION sharpen

END MODULE mesh_module
