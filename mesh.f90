! mesh.f90
! Assignment of charges to a 3-d mesh
PROGRAM meshup3d

  !     **************************************************************************************
  !     ** Assignment of charges to a three-dimensional mesh                                **
  !     **                                                                                  **
  !     ** This program assigns a set of charges to a cubic mesh using the triangular       **
  !     ** shape cloud distribution described by Hockney and Eastwood (1988)                **
  !     ** (section 5-3 and section 8-7-3). The charges are positioned in a box of unit     **
  !     ** length centred at the origin. The charge mesh is indexed from 0 to ncell-1       **
  !     ** in each dimension and there are ncell x ncell x ncell mesh cells. There are      **
  !     ** boundary faces at -1 and ncell that are used to implement the periodic boundary  **
  !     ** conditions during the assignment. These boundary layers need to be extended for  **
  !     ** higher order distributions algorithms involving more neighbouring mesh points.   **
  !     **                                                                                  **
  !     ** principal variables:                                                             **
  !     **                                                                                  **
  !     ** integer n                              number of atoms                           **
  !     ** integer ncell                          number of mesh cells in each dimension    **
  !     ** real    rx(n),ry(n),rz(n)              position of charges (-0.5 to 0.5)         **
  !     ** real    xp, yp, zp                     position of the cente of a mesh cell      **
  !     ** real    q(n)                           charge on each atom                       **
  !     ** real    h                              length of cubic mesh cell                 **
  !     ** real    dc(3,n)                        cartesian displacement of charge from     **
  !     **                                        the nearest mesh point                    **
  !     ** real    w(-1:ncell,-1:ncell,-1:ncell)  weights at the centre of the mesh cells   **
  !     ** real    rhom(-1:ncell,-1:cell,-1:cell) charge density at mesh cells              **
  !     ** integer nx,ny,nx                       nearest mesh point to a charge            **
  !     **                                                                                  **
  !     ** routines referenced:                                                             **
  !     **                                                                                  **
  !     ** subroutine faces                                                                 **
  !     **    adds the buffer faces of the mesh cube into their periodic images             **
  !     **************************************************************************************


  IMPLICIT NONE

  COMMON/block1/ w

  INTEGER, PARAMETER :: n= 4, ncell = 10, ncellm1 = ncell - 1

  REAL        ::  rx(n), ry(n), rz(n), q(n)
  REAL        ::  xi, yi, zi, qi, dc(3,n)
  REAL        ::  w(-1:ncell, -1:ncell, -1:ncell)
  REAL        ::  rhom(-1:ncell, -1:ncell, -1:ncell)
  REAL        ::  v(3, -1:1)
  INTEGER     ::  shift(0:ncell)
  INTEGER     ::  nx, ny, nz, n1, n2, n3
  INTEGER     ::  i, j, k, l, m
  REAL        ::  xp, yp, zp, h, p1, p2, p3, qq, sumw



  WRITE(*,'(1h1,'' a 3-D mesh assignment                        ''/ )')
  WRITE(*,'(1h1,'' for unit boxlength, centred at the origin    ''/ )')

  DO i = 1, n

     WRITE(*,'('' enter coordinates for atom '', i3)')  i
     READ(*,*) rx(i)
     READ(*,*) ry(i)
     READ(*,*) rz(i)
     WRITE(*,'('' enter charge for atom '', i3)')  i
     READ(*,*) q(i)

  ENDDO

  !     ** set up mesh **

  h = 1.0 / REAL( ncell )

  DO i = 0, ncell - 1

     shift(i) = i

  ENDDO

  shift(ncell) = 0


  !     ** zero mesh weights including the buffer zones **

  DO i = -1, ncell

     DO j = -1, ncell

        DO k = -1, ncell

           w(i,j,k) = 0.0

        ENDDO

     ENDDO

  ENDDO

  !     ** loop over all atoms (charges) begins **

  DO i = 1, n

     !        ** shift origin changing coordinates to 0 to 1 **

     xi = rx(i) + 0.5
     yi = ry(i) + 0.5
     zi = rz(i) + 0.5
     qi = q(i) / 8.0


     !        ** compute nearest mesh point with periodic boundaries **
     !        ** compute distance from charge to mesh cell centre with sign **
     !        ** then normalise this distance by the mesh cell size **

     nx       = shift( NINT(xi * ncell ))
     xp       = xi - REAL( nx ) * h
     dc(1,i)  = (xp - ANINT( xp )) / h

     ny       = shift( NINT(yi * ncell ))
     yp       = yi - REAL( ny ) * h
     dc(2,i)  = (yp - ANINT( yp )) / h

     nz       = shift( NINT(zi * ncell ))
     zp       = zi - REAL( nz ) * h
     dc(3,i)  = (zp - ANINT( zp )) / h

     !       ** compute weights for three point assignment scheme **

     DO j = 1, 3

        v(j,-1) = (0.5 - dc(j,i)) ** 2
        v(j, 0) = 1.5 - 2.0 * dc(j,i) * dc(j,i)
        v(j,+1) = (0.5 + dc(j,i)) ** 2

     ENDDO

     !        ** distribute charges over 27 nearest-neighbour mesh points **

     sumw = 0.0
     DO k = -1, 1

        p1 = qi * v(1,k)
        n1 = nx + k

        DO l = -1, 1

           p2 = p1 * v(2,l)
           n2 = ny + l

           DO m = -1, 1

              p3 = p2 * v(3,m)
              n3 = nz + m

              w(n1,n2,n3) = w(n1,n2,n3) + p3

           ENDDO

        ENDDO

     ENDDO


     !        ** end loop over charges **

  ENDDO

  !     ** adjust charge on mesh to include periodic buffer zones **

  CALL boundaries

  !       ** write out the final mesh: to be removed when happy

  DO k = 0, ncellm1

     WRITE(*,'('' layer '', i3)') k

     WRITE(*,'(10f8.5)')(( w(i,j,k), j= 0, ncellm1), i= 0, ncellm1)
     WRITE(*,'(/)')

  ENDDO

  !   ** as a check calculate the charge density on mesh **

  sumw = 0.0

  DO i = 0, ncellm1

     DO j = 0, ncellm1

        DO k = 0, ncellm1

           sumw        = sumw + w(i,j,k)
           rhom(i,j,k) = w(i,j,k)/ (h**3)

        ENDDO

     ENDDO

  ENDDO

  WRITE(*,'('' sum of all weights ='',e15.4)') sumw

END PROGRAM meshup3d


SUBROUTINE boundaries

  !     **************************************************************************************
  !     ** Subroutine boundaries adds the weights of face -1 into face ncell-1 and          **
  !     ** the weights of face ncell into face 0, considering the periodic boundaries       **
  !     ** of the mesh. It completes the six edge swaps and the eight corner swaps          **
  !     **                                                                                  **
  !     **                                                                                  **
  !     ** principal variables:                                                             **
  !     **                                                                                  **
  !     ** integer ncell                          number of mesh cells in each dimension    **
  !     ** real    w(-1:ncell,-1:ncell,-1:ncell)  weights at the centre of the mesh cells   **
  !     **************************************************************************************

  IMPLICIT NONE
  COMMON/block1/ w


  INTEGER, PARAMETER :: n = 4, ncell = 10, ncellm1 = ncell - 1

  REAL           :: w(-1: ncell , -1 : ncell, -1: ncell)
  INTEGER        :: i, j, k

  !     ** add top and bottom buffer faces into faces 0 and ncell-1 **

  DO j = 0, ncellm1

     DO k = 0 , ncellm1

        w(0,j,k)          = w(0,j,k)        + w(ncell,j,k)
        w(ncellm1,j,k)    = w(ncellm1,j,k)  + w(-1,j,k)

     ENDDO

  ENDDO

  !     ** add left and right buffer faces into faces 0 and ncell-1 **

  DO i = 0, ncellm1

     DO k = 0, ncellm1

        w(i,0,k)          = w(i,0,k)        + w(i,ncell,k)
        w(i,ncellm1,k)    = w(i,ncellm1,k)  + w(i,-1,k)

     ENDDO

  ENDDO

  !     ** add front and back buffer faces into faces 0 and ncell-1 **

  DO i = 0, ncellm1

     DO j = 0, ncellm1

        w(i,j,0)          = w(i,j,0)         + w(i,j,ncell)
        w(i,j,ncellm1)    = w(i,j,ncellm1)   + w(i,j,-1)

     ENDDO

  ENDDO

  !     ** perform the twelve edge swap **

  DO i = 0, ncellm1

     w(i,0,0)             = w(i,0,0)              + w(i,ncell,ncell)
     w(i,ncellm1,ncellm1) = w(i,ncellm1,ncellm1)  + w(i,-1,-1)
     w(i,0,ncellm1)       = w(i,0,ncellm1)        + w(i,ncell,-1)
     w(i,ncellm1,0)       = w(i,ncellm1,0)        + w(i,-1,ncell)

  ENDDO

  DO j = 0, ncellm1

     w(0,j,0)             = w(0,j,0)              + w(ncell,j,ncell)
     w(ncellm1,j,ncellm1) = w(ncellm1,j,ncellm1)  + w(-1,j,-1)
     w(0,j,ncellm1)       = w(0,j,ncellm1)        + w(ncell,j,-1)
     w(ncellm1,j,0)       = w(ncellm1,j,0)        + w(-1,j,ncell)

  ENDDO

  DO k = 0, ncellm1

     w(0,0,k)             = w(0,0,k)              + w(ncell,ncell,k)
     w(ncellm1,ncellm1,k) = w(ncellm1,ncellm1,k)  + w(-1,-1,k)
     w(0,ncellm1,k)       = w(0,ncellm1,k)        + w(ncell,-1,k)
     w(ncellm1,0,k)       = w(ncellm1,0,k)        + w(-1,ncell,k)

  ENDDO

  !     ** perform the eight corner swaps **

  w(0,0,0)                   = w(0,0,0)                    + w(ncell,ncell,ncell)
  w(ncellm1,ncellm1,ncellm1) = w(ncellm1,ncellm1,ncellm1)  + w(-1,-1,-1)
  w(0,ncellm1,ncellm1)       = w(0,ncellm1,ncellm1)        + w(ncell,-1,-1)
  w(ncellm1,0,0)             = w(ncellm1,0,0)              + w(-1,ncell,ncell)
  w(ncellm1,ncellm1,0)       = w(ncellm1,ncellm1,0)        + w(-1,-1,ncell)
  w(0,0,ncellm1)             = w(0,0,ncellm1)              + w(ncell,ncell,-1)
  w(0,ncellm1,0)             = w(0,ncellm1,0)              + w(ncell,-1,ncell)
  w(0,ncellm1,ncellm1)       = w(0,ncellm1,ncellm1)        + w(ncell,-1,-1)

END SUBROUTINE boundaries
