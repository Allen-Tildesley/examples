! force_lj.f90
SUBROUTINE force ( epslon, sigma, rcut, box, v, vc, w )

        COMMON / block1 / rx, ry, rz, vx, vy, vz, fx, fy, fz

c    *******************************************************************
c    ** force calculation for lennard-jones atoms.                    **
c    **                                                               **
c    ** in this we aim to show how the forces, potential energy and   **
c    ** virial FUNCTION are calculated in a fairly efficient way.     **
c    ** undoubtedly further improvement would be possible on specific **
c    ** machines.                                                     **
c    ** the potential is v(r) = 4*epslon*((sigma/r)**12-(sigma/r)**6) **
c    ** we INCLUDE spherical cutoff and minimum imaging in cubic box. **
c    ** the box length is box.  the cutoff is rcut.                   **
c    ** the routine actually returns two different potential energies.**
c    ** v is calculated using the lennard-jones potential to be used  **
c    ** for calculating the thermodynamic internal energy.            **
c    ** long-range corrections should be applied to this outside the  **
c    ** routine, in the form                                          **
c    **         sr3 = ( sigma / rcut ) ** 3                           **
c    **         sr9 = sr3 ** 3                                        **
c    **         dens = REAL(n) * ( sigma / box ) ** 3                 **
c    **         vlrc = ( 8.0 /9.0 ) * pi * epslon * dens * REAL ( n ) **
c    **      :           * ( sr9 - 3.0 * sr3 )                        **
c    **         wlrc = ( 16.0 / 9.0 ) * pi * epslon * dens * REAL( n )**
c    **      :           * ( 2.0 * sr9 - 3.0 * sr3 )                  **
c    **         v = v + vlrc                                          **
c    **         w = w + wlrc                                          **
c    ** vc is calculated using the shifted lennard-jones potential,   **
c    ** WITH no discontinuity at the cutoff, to be used in assessing  **
c    ** energy conservation.                                          **
c    ** no reduced units are used: for this potential we could set    **
c    ** epslon = 1 and either sigma = 1 or box = 1 to improve speed.  **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                 number of molecules                 **
c    ** REAL    rx(n),ry(n),rz(n) molecular positions                 **
c    ** REAL    vx(n),vy(n),vz(n) molecular velocities (not used)     **
c    ** REAL    fx(n),fy(n),fz(n) molecular forces                    **
c    ** REAL    sigma             pair potential length PARAMETER     **
c    ** REAL    epslon            pair potential energy PARAMETER     **
c    ** REAL    rcut              pair potential cutoff               **
c    ** REAL    box               simulation box length               **
c    ** REAL    v                 potential energy                    **
c    ** REAL    vc                shifted potential                   **
c    ** REAL    w                 virial FUNCTION                     **
c    ** REAL    vij               pair potential between i and j      **
c    ** REAL    wij               negative of pair virial FUNCTION w  **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )

        REAL        sigma, epslon, rcut, box, v, vc, w
        REAL        rx(n), ry(n), rz(n)
        REAL        vx(n), vy(n), vz(n)
        REAL        fx(n), fy(n), fz(n)

        INTEGER     i, j, ncut
        REAL        boxinv, rcutsq, sigsq, eps4, eps24
        REAL        rxi, ryi, rzi, fxi, fyi, fzi
        REAL        rxij, ryij, rzij, rijsq, fxij, fyij, fzij
        REAL        sr2, sr6, sr12, vij, wij, fij

c    *******************************************************************

c    ** calculate useful quantities **

        boxinv = 1.0 / box
        rcutsq = rcut ** 2
        sigsq  = sigma ** 2
        eps4   = epslon * 4.0
        eps24  = epslon * 24.0

c    ** zero forces, potential, virial **

        DO 100 i = 1, n

           fx(i) = 0.0
           fy(i) = 0.0
           fz(i) = 0.0

100     CONTINUE

        ncut = 0
        v    = 0.0
        w    = 0.0

c    ** outer loop begins **

        DO 200 i = 1, n - 1

           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)
           fxi = fx(i)
           fyi = fy(i)
           fzi = fz(i)

c       ** inner loop begins **

           DO 199 j = i + 1, n

              rxij = rxi - rx(j)
              ryij = ryi - ry(j)
              rzij = rzi - rz(j)
              rxij = rxij - ANINT ( rxij * boxinv ) * box
              ryij = ryij - ANINT ( ryij * boxinv ) * box
              rzij = rzij - ANINT ( rzij * boxinv ) * box
              rijsq = rxij ** 2 + ryij ** 2 + rzij ** 2

              IF ( rijsq .LT. rcutsq ) THEN

                 sr2   = sigsq / rijsq
                 sr6   = sr2 * sr2 * sr2
                 sr12  = sr6 ** 2
                 vij   = sr12 - sr6
                 v     = v + vij
                 wij   = vij + sr12
                 w     = w + wij
                 fij   = wij / rijsq
                 fxij  = fij * rxij
                 fyij  = fij * ryij
                 fzij  = fij * rzij
                 fxi   = fxi + fxij
                 fyi   = fyi + fyij
                 fzi   = fzi + fzij
                 fx(j) = fx(j) - fxij
                 fy(j) = fy(j) - fyij
                 fz(j) = fz(j) - fzij
                 ncut  = ncut + 1

              ENDIF

199        CONTINUE

c       ** inner loop ends **

           fx(i) = fxi
           fy(i) = fyi
           fz(i) = fzi

200     CONTINUE

c    ** outer loop ends **

c    ** calculate shifted potential **

        sr2 = sigsq / rcutsq
        sr6 = sr2 * sr2 * sr2
        sr12 = sr6 * sr6
        vij = sr12 - sr6
        vc = v - REAL ( ncut ) * vij

c    ** multiply results by energy factors **

        DO 300 i = 1, n

           fx(i) = fx(i) * eps24
           fy(i) = fy(i) * eps24
           fz(i) = fz(i) * eps24

300     CONTINUE

        v  = v  * eps4
        vc = vc * eps4
        w  = w  * eps24 / 3.0

        RETURN
        END



