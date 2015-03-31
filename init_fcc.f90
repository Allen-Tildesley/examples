! init_fcc.f90
SUBROUTINE init_fcc

        COMMON / block1 / rx, ry, rz, ex, ey, ez

c    *******************************************************************
c    ** sets up the alpha fcc lattice for n linear molecules.         **
c    **                                                               **
c    ** the simulation box is a unit cube centred at the origin.      **
c    ** n should be an INTEGER of the form ( 4 * ( nc ** 3 ) ),       **
c    ** WHERE nc is the number of fcc unit cells in each direction.   **
c    ** see figure 5.10 for a diagram of the lattice and a            **
c    ** definition of the four orientational sublattices.             **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                    number of molecules              **
c    ** REAL    rx(n),ry(n),rz(n)    molecular positions              **
c    ** REAL    ex(n),ey(n),ez(n)    unit vectors giving orientations **
c    ** REAL    rroot3               1.0 / SQRT ( 3.0 )               **
c    *******************************************************************

        INTEGER     n, nc
        REAL        rroot3

        PARAMETER ( nc = 3, n = 4 * nc ** 3 )
        PARAMETER ( rroot3 = 0.5773503 )

        REAL        rx(n), ry(n), rz(n), ex(n), ey(n), ez(n)
        REAL        cell, cell2
        INTEGER     i, ix, iy, iz, iref, m

c    *******************************************************************

c    ** calculate the side of the unit cell **

        cell  = 1.0 / REAL ( nc )
        cell2 = 0.5 * cell

c    ** build the unit cell **

c    ** sublattice a **

        rx(1) =  0.0
        ry(1) =  0.0
        rz(1) =  0.0
        ex(1) =  rroot3
        ey(1) =  rroot3
        ez(1) =  rroot3

c    ** sublattice b **

        rx(2) =  cell2
        ry(2) =  cell2
        rz(2) =  0.0
        ex(2) =  rroot3
        ey(2) = -rroot3
        ez(2) = -rroot3

c    ** sublattice c **

        rx(3) =  0.0
        ry(3) =  cell2
        rz(3) =  cell2
        ex(3) = -rroot3
        ey(3) =  rroot3
        ez(3) = -rroot3

c    ** sublattice d **

        rx(4) =  cell2
        ry(4) =  0.0
        rz(4) =  cell2
        ex(4) = -rroot3
        ey(4) = -rroot3
        ez(4) =  rroot3

c    ** construct the lattice from the unit cell **

        m = 0

        DO 99 iz = 1, nc

           DO 98 iy = 1, nc

              DO 97 ix = 1, nc

                 DO 96 iref = 1, 4

                    rx(iref+m) = rx(iref) + cell * REAL ( ix - 1 )
                    ry(iref+m) = ry(iref) + cell * REAL ( iy - 1 )
                    rz(iref+m) = rz(iref) + cell * REAL ( iz - 1 )

                    ex(iref+m) = ex(iref)
                    ey(iref+m) = ey(iref)
                    ez(iref+m) = ez(iref)

96               CONTINUE

                 m = m + 4

97            CONTINUE

98         CONTINUE

99      CONTINUE

c    ** shift centre of box to the origin **

        DO 100 i = 1, n

           rx(i) = rx(i) - 0.5
           ry(i) = ry(i) - 0.5
           rz(i) = rz(i) - 0.5

100     CONTINUE

        RETURN
        END



