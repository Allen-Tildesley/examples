! force_linked_list.f90
SUBROUTINE force ( sigma, rcut, v, w )

********************************************************************************
** fiche f.20.  routines to construct and USE cell linked-list method         **
** this fortran code is intended to illustrate points made in the text.       **
** to our knowledge it works correctly.  however it is the responsibility of  **
** the user to test it, IF it is to be used in a research application.        **
********************************************************************************

c    *******************************************************************
c    ** construction of cell linked-lists and USE in force routine.   **
c    **                                                               **
c    ** references:                                                   **
c    **                                                               **
c    ** quentrec and brot, j. comput. phys. 13, 430, 1975.            **
c    ** hockney and eastwood, computer simulation using particles,    **
c    **    mcgraw hill, 1981.                                         **
c    **                                                               **
c    ** routines supplied:                                            **
c    **                                                               **
c    ** SUBROUTINE maps                                               **
c    **    sets up map of cell structure for USE in force             **
c    ** SUBROUTINE links ( rcut )                                     **
c    **    sets up head of chain array and linked list                **
c    ** SUBROUTINE force ( sigma, rcut, v, w )                        **
c    **    calculates forces using a linked list                      **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** SUBROUTINE maps is called once at the start of a simulation   **
c    ** to establish cell neighbour identities.  at each timestep,    **
c    ** SUBROUTINE links is called to set up the linked list and this **
c    ** is immediately used by SUBROUTINE force.                      **
c    *******************************************************************



        SUBROUTINE maps

        COMMON / block2 / list, head, map

c    *******************************************************************
c    ** routine to set up a list of neighbouring cells                **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER m                  number of cells in each direction  **
c    ** INTEGER mapsiz             size of cell-cell map              **
c    ** INTEGER map(mapsiz)        list of neighbouring cells         **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** this SUBROUTINE sets up a list of the thirteen neighbouring   **
c    ** cells of each of the small cells in the central box. the      **
c    ** effects of the periodic boundary conditions are included.     **
c    ** the SUBROUTINE is called once at the beginning of the         **
c    ** simulation and the map is used in the force SUBROUTINE        **
c    *******************************************************************

        INTEGER     n, m, ncell, mapsiz
        PARAMETER ( n = 1372, m = 5, ncell = m * m * m )
        PARAMETER ( mapsiz = 13 * ncell )

        INTEGER     list(n), head(ncell), map(mapsiz)
        INTEGER     ix, iy, iz, imap, icell

c    *******************************************************************

c    ** statement FUNCTION to give cell index **

        icell ( ix, iy, iz) = 1 + MOD ( ix - 1 + m, m )
     :                          + MOD ( iy - 1 + m, m ) * m
     :                          + MOD ( iz - 1 + m, m ) * m * m

c    ** find half the nearest neighbours of each cell **

        DO 50 iz = 1, m

           DO 40 iy = 1, m

              DO 30 ix = 1, m

                 imap = ( icell ( ix, iy, iz ) - 1 ) * 13

                 map( imap + 1  ) = icell( ix + 1, iy    , iz     )
                 map( imap + 2  ) = icell( ix + 1, iy + 1, iz     )
                 map( imap + 3  ) = icell( ix    , iy + 1, iz     )
                 map( imap + 4  ) = icell( ix - 1, iy + 1, iz     )
                 map( imap + 5  ) = icell( ix + 1, iy    , iz - 1 )
                 map( imap + 6  ) = icell( ix + 1, iy + 1, iz - 1 )
                 map( imap + 7  ) = icell( ix    , iy + 1, iz - 1 )
                 map( imap + 8  ) = icell( ix - 1, iy + 1, iz - 1 )
                 map( imap + 9  ) = icell( ix + 1, iy    , iz + 1 )
                 map( imap + 10 ) = icell( ix + 1, iy + 1, iz + 1 )
                 map( imap + 11 ) = icell( ix    , iy + 1, iz + 1 )
                 map( imap + 12 ) = icell( ix - 1, iy + 1, iz + 1 )
                 map( imap + 13 ) = icell( ix    , iy    , iz + 1 )

30            CONTINUE

40         CONTINUE

50      CONTINUE

        RETURN
        END



        SUBROUTINE links ( rcut )

        COMMON / block1 / rx, ry, rz, vx, vy, vz, fx, fy, fz
        COMMON / block2 / list, head, map

c    *******************************************************************
c    ** routine to set up linked list and the head of chain arrays    **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                  number of atoms                    **
c    ** INTEGER m                  number of cells in each direction  **
c    ** INTEGER ncell              total number of cells (m**3)       **
c    ** INTEGER list(n)            linked list of atoms               **
c    ** INTEGER head(ncell)        head of chain for each cell        **
c    ** REAL    rx(n),ry(n),rz(n)  positions                          **
c    ** REAL    rcut               the cutoff distance for the force  **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** each atom is sorted into one of the m**3 small cells.         **
c    ** the first atom in each cell is placed in the head array.      **
c    ** subsequent atoms are placed in the linked list array.         **
c    ** atom coordinates are assumed to be between -0.5 and +0.5.     **
c    ** the routine is called every timestep before the force routine.**
c    *******************************************************************

        INTEGER     n, m, ncell, mapsiz
        PARAMETER ( n = 1372, m = 5, ncell = m * m * m )
        PARAMETER ( mapsiz = 13 * ncell )

        REAL        rx(n), ry(n), rz(n)
        REAL        vx(n), vy(n), vz(n)
        REAL        fx(n), fy(n), fz(n)
        INTEGER     head(ncell), list(n), map(mapsiz)

        REAL        celli, rcut, cell
        INTEGER     icell, i

c    *******************************************************************

c    ** zero head of chain array **

        DO 10 icell = 1, ncell

           head(icell) = 0

10      CONTINUE

        celli = REAL ( m )
        cell  = 1.0 / celli

        IF ( cell .LT. rcut ) THEN

           STOP ' cell size too small for cutoff '

        ENDIF

c    ** sort all atoms **

        DO 20 i = 1, n

           icell = 1 + INT ( ( rx(i) + 0.5 ) * celli )
     :               + INT ( ( ry(i) + 0.5 ) * celli ) * m
     :               + INT ( ( rz(i) + 0.5 ) * celli ) * m * m

           list(i)     = head(icell)
           head(icell) = i

20      CONTINUE

        RETURN
        END



        SUBROUTINE force ( sigma, rcut, v, w )

        COMMON / block1 / rx, ry, rz, vx, vy, vz, fx, fy, fz
        COMMON / block2 / list, head, map

c    *******************************************************************
c    ** routine to compute forces and potential using a link list     **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                  number of atoms                    **
c    ** INTEGER m                  number of cells in each direction  **
c    ** INTEGER ncell              number of small cells (m**3)       **
c    ** INTEGER mapsiz             size of cell-cell map              **
c    ** INTEGER list(n)            the linked list                    **
c    ** INTEGER head(ncell)        the head of chain array            **
c    ** INTEGER map(mapsiz)        list of neighbouring cells         **
c    ** REAL    rx(n),ry(n),rz(n)  positions                          **
c    ** REAL    fx(n),fy(n),fz(n)  forces                             **
c    ** REAL    sigma              the lj length PARAMETER            **
c    ** REAL    rcut               the cut-off distance               **
c    ** REAL    v                  the potential energy               **
c    ** REAL    w                  the virial                         **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** force is called in an md PROGRAM to calculate the force on    **
c    ** each atom. the routine is written for a liquid of lennard     **
c    ** jones atoms. SUBROUTINE force requires a linked list, set up  **
c    ** using SUBROUTINE links, and the map of the small cells set up **
c    ** using SUBROUTINE maps.                                        **
c    *******************************************************************

        INTEGER     n, m, ncell, mapsiz
        PARAMETER ( n = 1372, m = 5, ncell = m * m * m )
        PARAMETER ( mapsiz = 13 * ncell)

        REAL        rx(n), ry(n), rz(n)
        REAL        vx(n), vy(n), vz(n)
        REAL        fx(n), fy(n), fz(n)
        INTEGER     head(ncell), list(n), map(mapsiz)

        REAL        rcut, sigma, v, w

        REAL        rxi, ryi, rzi, fxij, fyij, fzij, fij, rcutsq
        REAL        sigsq, fxi, fyi, fzi, sr2, sr6, vij, wij
        REAL        rijsq, rxij, ryij, rzij
        INTEGER     icell, jcell0, jcell, i, j, nabor

c    *******************************************************************

        sigsq  = sigma ** 2
        rcutsq = rcut ** 2

c    ** zero forces and potential **

        DO 10 i = 1, n

           fx(i) = 0.0
           fy(i) = 0.0
           fz(i) = 0.0

10      CONTINUE

        v = 0.0
        w = 0.0

c    ** loop over all cells **

        DO 5000 icell = 1, ncell

           i = head(icell)

c       ** loop over all molecules in the cell **

1000       IF ( i .GT. 0 ) THEN

              rxi = rx(i)
              ryi = ry(i)
              rzi = rz(i)
              fxi = fx(i)
              fyi = fy(i)
              fzi = fz(i)

c          ** loop over all molecules below i in the current cell **

              j = list(i)

2000          IF ( j .GT. 0 ) THEN

                 rxij  = rxi - rx(j)
                 ryij  = ryi - ry(j)
                 rzij  = rzi - rz(j)

                 rxij  = rxij - ANINT ( rxij )
                 ryij  = ryij - ANINT ( ryij )
                 rzij  = rzij - ANINT ( rzij )
                 rijsq = rxij * rxij + ryij * ryij + rzij * rzij

                 IF ( rijsq .LT. rcutsq ) THEN

                    sr2   = sigsq / rijsq
                    sr6   = sr2 * sr2 * sr2
                    vij   = sr6 * ( sr6 - 1.0 )
                    wij   = sr6 * ( sr6 - 0.5 )
                    v     = v + vij
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

                 ENDIF

                 j = list(j)

                 go to 2000

              ENDIF

c          ** loop over neighbouring cells **

              jcell0 = 13 * (icell - 1)

              DO 4000 nabor = 1, 13

                 jcell = map ( jcell0 + nabor )

c             ** loop over all molecules in neighbouring cells **

                 j = head(jcell)

3000             IF ( j .NE. 0 ) THEN

                    rxij  = rxi - rx(j)
                    ryij  = ryi - ry(j)
                    rzij  = rzi - rz(j)

                    rxij  = rxij - ANINT ( rxij )
                    ryij  = ryij - ANINT ( ryij )
                    rzij  = rzij - ANINT ( rzij )
                    rijsq = rxij * rxij + ryij * ryij + rzij * rzij

                    IF ( rijsq. lt. rcutsq ) THEN

                       sr2   = sigsq / rijsq
                       sr6   = sr2 * sr2 * sr2
                       vij   = sr6 * ( sr6 - 1.0 )
                       wij   = sr6 * ( sr6 - 0.5 )
                       v     = v + vij
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

                    ENDIF

                    j = list(j)

                    go to 3000

                 ENDIF

4000          CONTINUE

              fx(i) = fxi
              fy(i) = fyi
              fz(i) = fzi

              i = list(i)

              go to 1000

           ENDIF

5000    CONTINUE

c    ** incorporate energy factors **

        v = 4.0  * v
        w = 48.0 * w / 3.0

        DO 50 i = 1, n

           fx(i) = 48.0 * fx(i)
           fy(i) = 48.0 * fy(i)
           fz(i) = 48.0 * fz(i)

50      CONTINUE

        RETURN
        END



