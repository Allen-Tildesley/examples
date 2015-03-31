! order_trans.f90
SUBROUTINE order_trans ( klatx, klaty, klatz, param )

        COMMON / block1 / rx, ry, rz, vx, vy, vz, fx, fy, fz

c    *******************************************************************
c    ** calculation of translational order PARAMETER (melting factor).**
c    **                                                               **
c    ** classically, the order PARAMETER is a normalized sum of       **
c    ** cosine terms which should be unity in the perfect lattice     **
c    ** and fluctuate around zero for a disordered system.            **
c    ** however, this is not origin-INDEPENDENT: WITH an unsuitable   **
c    ** choice of origin it could vanish even in a perfect lattice.   **
c    ** accordingly, we calculate here a quantity that is INDEPENDENT **
c    ** of the origin of coordinates.                                 **
c    ** it should be unity in a lattice for which a reciprocal vector **
c    ** (klatx,klaty,klatz) is supplied.                              **
c    ** it should be positive but small, of order SQRT(n) in a        **
c    ** disordered system.                                            **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                 number of molecules                 **
c    ** REAL    rx(n),ry(n),rz(n) molecular coordinates               **
c    ** REAL    vx(n),vy(n),vz(n) molecular velocities (not used)     **
c    ** REAL    fx(n),fy(n),fz(n) molecular forces (not used)         **
c    ** REAL    klatx,klaty,klatz reciproc. vector of initial lattice **
c    ** REAL    param             RESULT: order PARAMETER             **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )

        REAL        klatx, klaty, klatz, param
        REAL        rx(n), ry(n), rz(n)
        REAL        vx(n), vy(n), vz(n)
        REAL        fx(n), fy(n), fz(n)

        INTEGER     i
        REAL        sinsum, cossum

c    *******************************************************************

        sinsum = 0.0
        cossum = 0.0

        DO 100 i = 1, n

           cossum = cossum + COS (  klatx * rx(i)
     :                            + klaty * ry(i)
     :                            + klatz * rz(i) )
           sinsum = sinsum + SIN (  klatx * rx(i)
     :                            + klaty * ry(i)
     :                            + klatz * rz(i) )

100     CONTINUE

        cossum = cossum / REAL ( n )
        sinsum = sinsum / REAL ( n )
        param  = SQRT ( cossum ** 2 + sinsum ** 2 )

        RETURN
        END



