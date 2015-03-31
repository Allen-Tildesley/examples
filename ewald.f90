! ewald.f90
MODULE ewald

********************************************************************************
** fiche f.22.  routines to perform the ewald sum                             **
** this fortran code is intended to illustrate points made in the text.       **
** to our knowledge it works correctly.  however it is the responsibility of  **
** the user to test it, IF it is to be used in a research application.        **
********************************************************************************

c    *******************************************************************
c    ** REAL-space and reciprocal-space parts of ewald sum for ions.  **
c    **                                                               **
c    ** references:                                                   **
c    **                                                               **
c    ** woodcock and singer, trans. faraday soc. 67, 12, 1971.        **
c    ** de leeuw et al., proc. roy. soc. a 373, 27, 1980.             **
c    ** heyes, j. chem. phys. 74, 1924, 1981.                         **
c    ** see also fincham, mdions, ccp5 PROGRAM library.               **
c    **                                                               **
c    ** routines supplied:                                            **
c    **                                                               **
c    ** SUBROUTINE setup ( kappa )                                    **
c    **    sets up the wavevectors for USE in the ewald sum           **
c    ** SUBROUTINE rwald ( kappa, vr )                                **
c    **    calculates the r-space part of the sum                     **
c    ** SUBROUTINE kwald ( kappa, vk )                                **
c    **    calculates the k-space part of the sum                     **
c    ** REAL FUNCTION ERFC ( x )                                      **
c    **    returns the complementary error FUNCTION                   **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER  totk         the total number of k-vectors stored    **
c    ** INTEGER  maxk         maximum possible number of k-vectors    **
c    ** INTEGER  kmax         max INTEGER component of the k-vector   **
c    ** INTEGER  ksqmax       max square mod of the k-vector required **
c    ** REAL     vr           energy from r-space sum                 **
c    ** REAL     vk           energy from k-space sum                 **
c    ** REAL     kvec(maxk)   array used to store k-vectors           **
c    ** REAL     kappa        width of cancelling distribution        **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** setup is called once at the beginning of the simulation       **
c    ** to calculate all the k-vectors required in the ewald sum.     **
c    ** these vectors are used throughout the simulation in the       **
c    ** SUBROUTINE kwald to calculate the k-space contribution to the **
c    ** potential energy at each configuration. the self term is      **
c    ** subtracted from the k-space contribution in kwald.            **
c    ** the surface term for simulations in vacuum is not included.   **
c    ** routine rwald returns the r-space contribution to the ewald   **
c    ** sum and is called for each configuration in the simulation.   **
c    ** a cubic box and unit box length are assumed throughout.       **
c    *******************************************************************



        SUBROUTINE setup ( kappa )

        COMMON / block2 / kvec

c    *******************************************************************
c    ** routine to set up the wave-vectors for the ewald sum.         **
c    **                                                               **
c    ** the wavevectors must fit into a box of unit length.           **
c    ** in this example we allow a maximum of 1000 wavevectors.       **
c    *******************************************************************

        INTEGER     maxk
        PARAMETER ( maxk = 1000 )

        REAL        kvec(maxk), kappa

        INTEGER     kmax, ksqmax, ksq, kx, ky, kz, totk
        REAL        twopi, b, rkx, rky, rkz, rksq
        PARAMETER ( kmax = 5, ksqmax = 27 , twopi = 6.2831853 )

c    *******************************************************************

        b = 1.0 / 4.0 / kappa / kappa

c    ** loop over k-vectors. note kx is non-negative **

        totk = 0

        DO 100 kx = 0, kmax

           rkx = twopi * REAL ( kx )

           DO 99 ky = -kmax, kmax

              rky = twopi * REAL ( ky )

              DO 98 kz = -kmax, kmax

                 rkz = twopi * REAL ( kz )

                 ksq = kx * kx + ky * ky + kz * kz

                 IF ( ( ksq .LT. ksqmax ) .AND. ( ksq .NE. 0 ) ) THEN

                    totk = totk + 1

                    IF ( totk .GT. maxk ) STOP 'kvec is too small'

                    rksq = rkx * rkx + rky * rky + rkz * rkz
                    kvec(totk) = twopi * EXP ( -b * rksq ) / rksq

                 ENDIF

98            CONTINUE

99         CONTINUE

100     CONTINUE

        WRITE( *, ' ( '' ewald sum setup complete ''     ) ' )
        WRITE( *, ' ( '' number of wavevectors is '', i5 ) ' ) totk

        RETURN
        END



        SUBROUTINE rwald ( kappa, vr )

        COMMON / block1 / rx, ry, rz, z

c    *******************************************************************
c    ** calculates r-space part of potential energy by ewald method.  **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                     number of ions                  **
c    ** REAL    rx(n),ry(n),rz(n)     positions of ions               **
c    ** REAL    z(n)                  ionic charges                   **
c    ** REAL    vr                    r-space potential energy        **
c    **                                                               **
c    ** routine referenced:                                           **
c    **                                                               **
c    ** REAL FUNCTION ERFC ( x )                                      **
c    **    returns the complementary error FUNCTION                   **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 216 )

        REAL        rx(n), ry(n), rz(n), z(n)
        REAL        kappa, vr

        REAL        rxi, ryi, rzi, zi, rxij, ryij, rzij
        REAL        rijsq, rij, krij, erfc, vij

        INTEGER     i, j

c    *******************************************************************

        vr = 0.0

        DO 100 i = 1, n - 1

           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)
           zi  = z(i)

           DO 99 j = i + 1, n

              rxij = rxi - rx(j)
              ryij = ryi - ry(j)
              rzij = rzi - rz(j)

              rxij = rxij - ANINT ( rxij )
              ryij = ryij - ANINT ( ryij )
              rzij = rzij - ANINT ( rzij )

              rijsq = rxij * rxij + ryij * ryij + rzij * rzij
              rij   = SQRT ( rijsq )
              krij  = kappa * rij
              vij   = zi * z(j) * ERFC ( krij ) / rij

              vr    = vr + vij

99         CONTINUE

100     CONTINUE

        RETURN
        END



        SUBROUTINE kwald ( kappa, vk )

        COMMON / block1 / rx, ry, rz, z
        COMMON / block2 / kvec

c    *******************************************************************
c    ** calculates k-space part of potential energy by ewald method.  **
c    **                                                               **
c    ** the self term is subtracted.                                  **
c    ** in one coordinate direction (x), symmetry is used to reduce   **
c    ** the sum to INCLUDE ONLY positive k-vectors.                   **
c    ** the negative vectors in this direction are included by USE    **
c    ** of the multiplicative variable 'factor'.                      **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                   number of ions                    **
c    ** REAL    rx(n),ry(n),rz(n)   positions of ions                 **
c    ** REAL    z(n)                ionic charges                     **
c    ** REAL    vk                  k-space potential energy          **
c    ** REAL    vks                 self part of k-space sum          **
c    *******************************************************************

        INTEGER     maxk, n
        PARAMETER ( maxk = 1000, n = 216 )

        REAL        kvec(maxk), rx(n), ry(n), rz(n), z(n)
        REAL        kappa, vk
        INTEGER     totk

        INTEGER     kmax, kx, ky, kz, i, ksqmax, ksq
        REAL        twopi, factor, vd, vs, rsqpi
        PARAMETER ( kmax = 5, ksqmax = 27 )
        PARAMETER ( twopi = 6.2831853, rsqpi = 0.5641896 )

        COMPLEX     eikx(1:n, 0:kmax)
        COMPLEX     eiky(1:n, -kmax:kmax)
        COMPLEX     eikz(1:n, -kmax:kmax)
        COMPLEX     eikr(n), sum

c    *******************************************************************

c    ** construct EXP(ik.r) for all ions and k-vectors **

c    ** calculate kx, ky, kz = 0 , -1 and 1 explicitly **

        DO 10 i = 1, n

           eikx(i, 0) = (1.0, 0.0)
           eiky(i, 0) = (1.0, 0.0)
           eikz(i, 0) = (1.0, 0.0)

           eikx(i, 1) = CMPLX ( COS ( twopi * rx(i) ) ,
     :                          SIN ( twopi * rx(i) ) )
           eiky(i, 1) = CMPLX ( COS ( twopi * ry(i) ) ,
     :                          SIN ( twopi * ry(i) ) )
           eikz(i, 1) = CMPLX ( COS ( twopi * rz(i) ) ,
     :                          SIN ( twopi * rz(i) ) )

           eiky(i, -1) = CONJG ( eiky(i, 1) )
           eikz(i, -1) = CONJG ( eikz(i, 1) )

10      CONTINUE

c    ** calculate remaining kx, ky and kz by recurrence **

        DO 12 kx = 2, kmax

           DO 11 i = 1, n

              eikx(i, kx) = eikx(i, kx-1) * eikx(i, 1)

11         CONTINUE

12      CONTINUE

        DO 14 ky = 2, kmax

           DO 13 i = 1, n

              eiky(i,  ky) = eiky(i, ky-1) * eiky(i, 1)
              eiky(i, -ky) = CONJG ( eiky(i, ky) )

13         CONTINUE

14      CONTINUE

        DO 16 kz = 2, kmax

           DO 15 i = 1, n

              eikz(i,  kz) = eikz(i, kz-1) * eikz(i, 1)
              eikz(i, -kz) = CONJG ( eikz(i, kz) )

15         CONTINUE

16      CONTINUE

c    ** sum over all vectors **

        vd   = 0.0
        totk = 0

        DO 24 kx = 0, kmax

           IF ( kx .EQ. 0 ) THEN

              factor = 1.0

           ELSE

              factor = 2.0

           ENDIF

           DO 23 ky = -kmax, kmax

              DO 22 kz = -kmax, kmax

                 ksq = kx * kx + ky * ky + kz * kz

                 IF ( ( ksq .LT. ksqmax ) .AND. ( ksq .NE. 0 ) ) THEN

                    totk = totk + 1
                    sum  = (0.0, 0.0)

                    DO 21 i = 1, n

                       eikr(i) = eikx(i, kx) * eiky(i, ky) * eikz(i, kz)
                       sum     = sum + z(i) * eikr(i)

21                  CONTINUE

                    vd = vd + factor * kvec(totk) * CONJG ( sum ) * sum

                 ENDIF

22            CONTINUE

23         CONTINUE

24      CONTINUE

c    ** calculates self part of k-space sum **

        vs = 0.0

        DO 25 i = 1, n

           vs = vs + z(i) * z(i)

25      CONTINUE

        vs = rsqpi * kappa * vs

c    ** calculate the total k-space potential **

        vk = vd - vs

        RETURN
        END



        REAL FUNCTION ERFC ( x )

c    *******************************************************************
c    ** approximation to the complementary error FUNCTION             **
c    **                                                               **
c    ** reference:                                                    **
c    **                                                               **
c    ** abramowitz and stegun, handbook of mathematical functions,    **
c    **    national bureau of standards, formula 7.1.26               **
c    *******************************************************************

        REAL        a1, a2, a3, a4, a5, p

        PARAMETER ( a1 = 0.254829592, a2 = -0.284496736 )
        PARAMETER ( a3 = 1.421413741, a4 = -1.453152027 )
        PARAMETER ( a5 = 1.061405429, p  =  0.3275911   )

        REAL        t, x, xsq, tp

c    *******************************************************************

        t  = 1.0 / ( 1.0 + p * x )
        xsq = x * x

        tp = t * ( a1 + t * ( a2 + t * ( a3 + t * ( a4 + t * a5 ) ) ) )

        erfc = tp * EXP ( -xsq )

        RETURN
        END



