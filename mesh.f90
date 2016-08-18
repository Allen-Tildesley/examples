! mesh.f90
! Assignment of charges to a 3-d mesh
PROGRAM mesh

  ! This program assigns a set of charges to a cubic mesh using the
  ! triangular shape cloud distribution described by Hockney and Eastwood (1988)
  ! The charges are positioned in a box of unit length.
  ! The charge mesh is indexed from 0 to sc-1 in each coordinate direction

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit
  IMPLICIT NONE

  INTEGER                             :: n   ! number of charges
  INTEGER                             :: sc  ! dimension of mesh
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rho ! (0:sc-1,0:sc-1,0:sc-1) mesh cell charge density
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: r   ! (3,n) charge positions in range (0,1)
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: dr  ! (3,n) charge positions relative to mesh point
  REAL, ALLOCATABLE, DIMENSION(:)     :: q   ! (n) charges

  REAL,    DIMENSION(3,-1:1) :: v  ! weights in each coordinate direction
  INTEGER, DIMENSION(3)      :: nr ! mesh point index

  INTEGER :: i, i1, i2, i3, n1, n2, n3
  REAL    :: h, q1, q2, q3

  ! Example values for illustration
  n  = 4                ! number of charges
  sc = 10               ! dimension of mesh
  h  = 1.0 / REAL( sc ) ! mesh spacing

  ALLOCATE ( rho(0:sc-1,0:sc-1,0:sc-1) ) ! C-style indexing is convenient here
  ALLOCATE ( r(3,n), dr(3,n), q(n) )

  WRITE ( unit=output_unit, fmt='(a)' ) '3-D mesh assignment for unit boxlength, coordinates in range (0,1)'

  ! For illustration we choose random charge positions with coordinates in range (0,1)
  ! In a real application, we would convert positions into this range
  CALL RANDOM_SEED() ! same random number sequence every time
  CALL RANDOM_NUMBER ( r )

  ! For illustration we choose +1 and -1 charges, alternately
  q(1::2) = 1.0
  q(2::2) = -1.0

  rho = 0.0 ! zero mesh

  DO i = 1, n ! Loop over charges

     nr(:)   = NINT ( r(:,i) * sc )       ! nearest mesh point indices
     nr(:)   = MODULO ( nr, sc)           ! with periodic boundaries
     dr(:,i) = r(:,i) - REAL( nr(:) ) * h ! vector to charge from mesh cell centre
     dr(:,i) = dr(:,i) - ANINT( dr(:,i) ) ! periodic boundaries
     dr(:,i) = dr(:,i) / h                ! normalise by mesh cell size

     ! weights for three point assignment scheme
     v(:,-1) = 0.5 * ( 0.5 - dr(:,i) ) ** 2
     v(:, 0) = 0.75 - dr(:,i)**2
     v(:,+1) = 0.5 * ( 0.5 + dr(:,i) ) ** 2

     DO i1 = -1, 1 ! Loop over x-neighbours
        q1 = q(i) * v(1,i1)            ! charge times x-weight
        n1 = MODULO ( nr(1) + i1, sc ) ! neighbour index with periodic boundaries

        DO i2 = -1, 1 ! Loop over y-neighbours
           q2 = q1 * v(2,i2)              ! charge times xy-weight
           n2 = MODULO ( nr(2) + i2, sc ) ! neighbour index with periodic boundaries

           DO i3 = -1, 1 ! Loop over z-neighbours
              q3 = q2 * v(3,i3)                  ! charge times xyz-weight
              n3 = MODULO ( nr(3) + i3, sc )     ! neighbour index with periodic boundaries
              rho(n1,n2,n3) = rho(n1,n2,n3) + q3 ! mesh cell share of charge
           END DO ! End loop over z-neighbours

        END DO ! End loop over y-neighbours

     END DO ! End loop over x-neighbours

  END DO ! End loop over charges

  rho = rho / (h**3) ! Convert charges to charge densities

  ! Output mesh charge density
  DO n3 = 0, sc-1
     WRITE( unit=output_unit, fmt='(a,i5)' ) 'z-layer ', n3
     DO n2 = 0, sc-1
        WRITE( unit=output_unit, fmt='(*(f10.4))') rho(:,n2,n3)
     END DO
  END DO

  ! Finally check integrated charge density
  WRITE( unit=output_unit, fmt='(a,2f10.6)') 'Total charge = ', SUM ( q ), SUM ( rho )*(h**3)

  DEALLOCATE ( r, dr, q, rho )
  
END PROGRAM mesh
