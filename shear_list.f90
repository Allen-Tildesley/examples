! shear_list.f90
module shear_list
! Placeholder for f90 module to handle cell linked lists in sheared boundaries

********************************************************************************
** FICHE F.32.  CELL LINKED-LISTS IN SHEARED BOUNDARIES.                      **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** ROUTINES TO IMPLEMENT CELL LINKED-LISTS IN SHEARED BOUNDARIES.**
C    **                                                               **
C    ** ROUTINES PROVIDED:                                            **
C    **                                                               **
C    ** SUBROUTINE MAPS                                               **
C    **    SETS UP MAP OF CELL NEIGHBOURS FOR BULK OF SIMULATION BOX  **
C    ** SUBROUTINE TOPMAP ( STRAIN )                                  **
C    **    SETS UP MAP OF NEIGHBOURS FOR TOP LAYER OF CELLS           **
C    ** SUBROUTINE LINKS ( RCUT )                                     **
C    **    CONSTRUCTS LINK-LIST GIVEN MAP OF CELL NEIGHBOURS          **
C    ** SUBROUTINE FORCE ( SIGMA, RCUT, STRAIN, V, W, WXY )           **
C    **    CALCULATES FORCES, POTENTIAL, VIRIAL ETC. USING LIST       **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF ATOMS                     **
C    ** REAL    RX(N),RY(N),RZ(N) ATOMIC POSITIONS                    **
C    ** REAL    VX(N),VY(N),VZ(N) ATOMIC VELOCITIES                   **
C    ** REAL    FX(N),FY(N),FZ(N) ATOMIC FORCES                       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** SUBROUTINE MAPS IS CALLED ONCE AT THE START OF THE SIMULATION **
C    ** TO DEFINE CELL NEIGHBOURS FOR ALL BUT THE TOP LAYER OF CELLS. **
C    ** AT EACH TIME STEP, SUBROUTINE TOPMAP IS CALLED FOR THE TOP    **
C    ** LAYER, SUBROUTINE LINKS TO ESTABLISH THE ATOM NEIGHBOUR LIST, **
C    ** AND THEN THE FORCE SUBROUTINE.                                **
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** THE PROGRAM ASSUMES A BOX OF UNIT LENGTH AND TAKES THE        **
C    ** LENNARD-JONES POTENTIAL WITH UNIT WELL-DEPTH.                 **
C    ** SUMMARY FOR BOX LENGTH L, ATOMIC MASS M, AND LENNARD-JONES    **
C    ** POTENTIAL PARAMETERS SIGMA AND EPSILON:                       **
C    **                                                               **
C    **                OUR PROGRAM            LENNARD-JONES SYSTEM    **
C    ** LENGTH         L                      SIGMA                   **
C    ** MASS           M                      M                       **
C    ** ENERGY         EPSILON                EPSILON                 **
C    ** TIME           SQRT(M*L**2/EPSILON)   SQRT(M*SIGMA**2/EPSILON)**
C    ** VELOCITY       SQRT(EPSILON/M)        SQRT(EPSILON/M)         **
C    ** PRESSURE       EPSILON/L**3           EPSILON/SIGMA**3        **
C    *******************************************************************



        SUBROUTINE MAPS

        COMMON / BLOCK2 / LIST, HEAD, MAP

C    *******************************************************************
C    ** CONSTRUCTS MAP OF CELL NEIGHBOURS.                            **
C    **                                                               **
C    ** THIS SUBROUTINE SETS UP A LIST OF THE THIRTEEN NEIGHBOURING   **
C    ** CELLS OF EACH OF THE SMALL CELLS IN THE CENTRAL BOX. THE      **
C    ** EFFECTS OF THE PERIODIC BOUNDARY CONDITIONS ARE INCLUDED.     **
C    ** HOWEVER THE TOP LAYER (IY = M) IS TACKLED SEPARATELY BECAUSE  **
C    ** OF THE SHEARED BOUNDARY CONDITIONS, IN SUBROUTINE TOPMAP.     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
C    ** INTEGER MAPSIZ             SIZE OF CELL-CELL MAP              **
C    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE SUBROUTINE IS CALLED ONCE AT THE BEGINNING OF THE         **
C    ** SIMULATION AND THE MAP IS USED IN THE FORCE SUBROUTINE        **
C    *******************************************************************

        INTEGER     N, M, NCELL, MAPSIZ, M3
        PARAMETER ( N = 1372 )
        PARAMETER ( M = 5, NCELL = M * M * M, MAPSIZ = 16 * NCELL )
        PARAMETER ( M3 = M * 3 )

        INTEGER     LIST(N), HEAD(NCELL), MAP(MAPSIZ)
        INTEGER     IX, IY, IZ, IMAP, ICELL

C    *******************************************************************

C    ** STATEMENT FUNCTION TO GIVE CELL INDEX **

        ICELL ( IX, IY, IZ ) = 1 + MOD ( IX - 1 + M3, M )
     :                           + MOD ( IY - 1 + M3, M ) * M
     :                           + MOD ( IZ - 1 + M3, M ) * M * M

C    ** FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **

        DO 50 IZ = 1, M

           DO 40 IY = 1, M - 1

              DO 30 IX = 1, M

                 IMAP = ( ICELL ( IX, IY, IZ ) - 1 ) * 16

                 MAP( IMAP + 1  ) = ICELL ( IX + 1, IY    , IZ     )
                 MAP( IMAP + 2  ) = ICELL ( IX + 1, IY + 1, IZ     )
                 MAP( IMAP + 3  ) = ICELL ( IX    , IY + 1, IZ     )
                 MAP( IMAP + 4  ) = ICELL ( IX - 1, IY + 1, IZ     )
                 MAP( IMAP + 5  ) = ICELL ( IX + 1, IY    , IZ - 1 )
                 MAP( IMAP + 6  ) = ICELL ( IX + 1, IY + 1, IZ - 1 )
                 MAP( IMAP + 7  ) = ICELL ( IX    , IY + 1, IZ - 1 )
                 MAP( IMAP + 8  ) = ICELL ( IX - 1, IY + 1, IZ - 1 )
                 MAP( IMAP + 9  ) = ICELL ( IX + 1, IY    , IZ + 1 )
                 MAP( IMAP + 10 ) = ICELL ( IX + 1, IY + 1, IZ + 1 )
                 MAP( IMAP + 11 ) = ICELL ( IX    , IY + 1, IZ + 1 )
                 MAP( IMAP + 12 ) = ICELL ( IX - 1, IY + 1, IZ + 1 )
                 MAP( IMAP + 13 ) = ICELL ( IX    , IY    , IZ + 1 )
                 MAP( IMAP + 14 ) = 0
                 MAP( IMAP + 15 ) = 0
                 MAP( IMAP + 16 ) = 0

30            CONTINUE

40         CONTINUE

50      CONTINUE

        RETURN
        END



        SUBROUTINE TOPMAP ( STRAIN )

        COMMON / BLOCK2 / LIST, HEAD, MAP

C    *******************************************************************
C    ** CALCULATES CELL NEIGHBOUR MAP FOR TOP LAYER.                  **
C    **                                                               **
C    ** THIS SUBROUTINE SUPPLEMENTS THE LIST OF NEIGHBOURING CELLS    **
C    ** FOR THE TOP LAYER (IY = M) WITH SHEARED BOUNDARY CONDITIONS   **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
C    ** INTEGER MAPSIZ             SIZE OF CELL-CELL MAP              **
C    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
C    ** REAL    STRAIN             THE X-DISPLACEMENT OF NEXT BOX UP  **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE SUBROUTINE IS CALLED AT EVERY TIMESTEP IN THE SIMULATION  **
C    ** JUST BEFORE THE FORCE SUBROUTINE                              **
C    *******************************************************************

        INTEGER     N, M, NCELL, MAPSIZ, M3
        PARAMETER ( N = 1372 )
        PARAMETER ( M = 5, NCELL = M * M * M, MAPSIZ = 16 * NCELL )
        PARAMETER ( M3 = M * 3 )

        REAL        STRAIN

        INTEGER     LIST(N), HEAD(NCELL), MAP(MAPSIZ)
        INTEGER     IX, IY, IZ, IMAP, ICELL, IIX

C    *******************************************************************

C    ** STATEMENT FUNCTION TO GIVE CELL INDEX **

        ICELL ( IX, IY, IZ ) = 1 + MOD ( IX - 1 + M3, M )
     :                           + MOD ( IY - 1 + M3, M ) * M
     :                           + MOD ( IZ - 1 + M3, M ) * M * M

C    ** CALCULATE X OFFSET IN CELL LENGTHS WHERE STRAIN **
C    ** IS BETWEEN -1/2 AND +1/2 AND BOX LENGTH = 1.0   **
C    ** ADDING 1.0 SIMPLY GUARANTEES A POSITIVE RESULT  **

        STRAIN = STRAIN - ANINT ( STRAIN )
        IIX   = INT ( ( STRAIN + 1.0 ) * REAL ( M ) )

C    ** FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **

        IY = M

        DO 50 IZ = 1, M

           DO 30 IX = 1, M

              IMAP = ( ICELL ( IX, IY, IZ ) - 1 ) * 16

              MAP( IMAP + 1  ) = ICELL ( IX + 1      , IY    , IZ     )
              MAP( IMAP + 2  ) = ICELL ( IX + 1 - IIX, IY + 1, IZ     )
              MAP( IMAP + 3  ) = ICELL ( IX     - IIX, IY + 1, IZ     )
              MAP( IMAP + 4  ) = ICELL ( IX - 1 - IIX, IY + 1, IZ     )
              MAP( IMAP + 5  ) = ICELL ( IX + 1      , IY    , IZ - 1 )
              MAP( IMAP + 6  ) = ICELL ( IX + 1 - IIX, IY + 1, IZ - 1 )
              MAP( IMAP + 7  ) = ICELL ( IX     - IIX, IY + 1, IZ - 1 )
              MAP( IMAP + 8  ) = ICELL ( IX - 1 - IIX, IY + 1, IZ - 1 )
              MAP( IMAP + 9  ) = ICELL ( IX + 1      , IY    , IZ + 1 )
              MAP( IMAP + 10 ) = ICELL ( IX + 1 - IIX, IY + 1, IZ + 1 )
              MAP( IMAP + 11 ) = ICELL ( IX     - IIX, IY + 1, IZ + 1 )
              MAP( IMAP + 12 ) = ICELL ( IX - 1 - IIX, IY + 1, IZ + 1 )
              MAP( IMAP + 13 ) = ICELL ( IX          , IY    , IZ + 1 )
              MAP( IMAP + 14 ) = ICELL ( IX - 2 - IIX, IY + 1, IZ     )
              MAP( IMAP + 15 ) = ICELL ( IX - 2 - IIX, IY + 1, IZ - 1 )
              MAP( IMAP + 16 ) = ICELL ( IX - 2 - IIX, IY + 1, IZ + 1 )

30         CONTINUE

50      CONTINUE

        RETURN
        END



        SUBROUTINE LINKS ( RCUT )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / LIST, HEAD, MAP

C    *******************************************************************
C    ** SUBROUTINE TO SET UP LINKED LIST AND THE HEAD OF CHAIN ARRAYS **
C    **                                                               **
C    ** EACH ATOM IS SORTED INTO ONE OF THE M**3 SMALL CELLS.         **
C    ** THE FIRST ATOM IN EACH CELL IS PLACED IN THE HEAD ARRAY.      **
C    ** SUBSEQUENT ATOMS ARE PLACED IN THE LINKED LIST ARRAY.         **
C    ** ATOM COORDINATES ARE ASSUMED TO BE BETWEEN -0.5 AND +0.5.     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                  NUMBER OF ATOMS                    **
C    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
C    ** INTEGER NCELL              TOTAL NUMBER OF CELLS (M**3)       **
C    ** INTEGER LIST(N)            LINKED LIST OF ATOMS               **
C    ** INTEGER HEAD(NCELL)        HEAD OF CHAIN FOR EACH CELL        **
C    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
C    ** REAL    RCUT               THE CUTOFF DISTANCE FOR THE FORCE  **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE ROUTINE IS CALLED EVERY TIMESTEP BEFORE THE FORCE ROUTINE.**
C    *******************************************************************

        INTEGER     N, M, NCELL, MAPSIZ
        PARAMETER ( N = 1372 )
        PARAMETER ( M = 5, NCELL = M * M * M, MAPSIZ = 16 * NCELL )

        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        INTEGER     HEAD(NCELL), LIST(N), MAP(MAPSIZ)

        REAL        CELLI, RCUT, CELL
        INTEGER     ICELL, I

C    *******************************************************************

C    ** ZERO HEAD OF CHAIN ARRAY **

        DO 10 ICELL = 1, NCELL

           HEAD(ICELL) = 0

10      CONTINUE

        CELLI = REAL ( M )
        CELL  = 1.0 / CELLI

        IF ( CELL. LT. RCUT ) THEN

           WRITE(*,'('' CELL SIZE TOO SMALL FOR CUTOFF '')')
           STOP

        ENDIF

C    ** SORT ALL ATOMS **

        DO 20 I = 1, N

           ICELL = 1 + INT ( ( RX(I) + 0.5 ) * CELLI )
     :               + INT ( ( RY(I) + 0.5 ) * CELLI ) * M
     :               + INT ( ( RZ(I) + 0.5 ) * CELLI ) * M * M

           LIST(I)     = HEAD(ICELL)
           HEAD(ICELL) = I

20      CONTINUE

        RETURN
        END



        SUBROUTINE FORCE ( SIGMA, RCUT, STRAIN, V, W, WXY )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / LIST, HEAD, MAP

C    *******************************************************************
C    ** COMPUTES FORCES, ETC. USING A LINK LIST IN SHEARED BOUNDARIES.**
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                  NUMBER OF ATOMS                    **
C    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
C    ** INTEGER NCELL              NUMBER OF SMALL CELLS (M**3)       **
C    ** INTEGER MAPSIZ             SIZE OF CELL-CELL MAP              **
C    ** INTEGER LIST(N)            THE LINKED LIST                    **
C    ** INTEGER HEAD(NCELL)        THE HEAD OF CHAIN ARRAY            **
C    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
C    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
C    ** REAL    FX(N),FY(N),FZ(N)  FORCES                             **
C    ** REAL    SIGMA              THE LJ LENGTH PARAMETER            **
C    ** REAL    RCUT               THE CUT-OFF DISTANCE               **
C    ** REAL    STRAIN             X OFFSET OF SUCCESSIVE BOXES       **
C    ** REAL    V                  THE POTENTIAL ENERGY               **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** FORCE IS CALLED IN AN MD PROGRAM TO CALCULATE THE FORCE ON    **
C    ** EACH ATOM. THE ROUTINE IS WRITTEN FOR A LIQUID OF LENNARD     **
C    ** JONES ATOMS. SUBROUTINE FORCE REQUIRES A LINKED LIST SET UP   **
C    ** USING SUBROUTINE LINKS AND THE MAP OF THE SMALL CELLS SET UP  **
C    ** USING SUBROUTINES MAPS AND TOPMAP.                            **
C    *******************************************************************

        INTEGER     N, M, NCELL, MAPSIZ
        PARAMETER ( N = 1372 )
        PARAMETER ( M = 5, NCELL = M * M * M, MAPSIZ = 16 * NCELL )

        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        INTEGER     HEAD(NCELL), LIST(N), MAP(MAPSIZ)

        REAL        RCUT, SIGMA, STRAIN, V, W, WXY

        REAL        RXI, RYI, RZI, FXIJ, FYIJ, FZIJ, RCUTSQ
        REAL        VIJ, WIJ, FIJ
        REAL        SIGSQ, FXI, FYI, FZI, SR2, SR6, SR12
        REAL        RIJSQ, RXIJ, RYIJ, RZIJ, CORY
        INTEGER     ICELL, JCELL0, JCELL, I, J, NABOR

C    *******************************************************************

        SIGSQ  = SIGMA ** 2
        RCUTSQ = RCUT ** 2

C    ** ZERO FORCES AND POTENTIAL **

        DO 10 I = 1, N

           FX(I) = 0.0
           FY(I) = 0.0
           FZ(I) = 0.0

10      CONTINUE

        V    = 0.0
        W    = 0.0
        WXY  = 0.0

C    ** LOOP OVER ALL CELLS **

        DO 5000 ICELL = 1, NCELL

           I = HEAD(ICELL)

C       ** LOOP OVER ALL MOLECULES IN THE CELL **

1000       IF ( I .GT. 0 ) THEN

              RXI = RX(I)
              RYI = RY(I)
              RZI = RZ(I)
              FXI = FX(I)
              FYI = FY(I)
              FZI = FZ(I)

C          ** LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL **

              J = LIST(I)

2000          IF ( J .GT. 0 ) THEN

                 RXIJ  = RXI - RX(J)
                 RYIJ  = RYI - RY(J)
                 RZIJ  = RZI - RZ(J)
                 CORY  = ANINT ( RYIJ )
                 RXIJ  = RXIJ - CORY * STRAIN
                 RXIJ  = RXIJ - ANINT ( RXIJ )
                 RYIJ  = RYIJ - CORY
                 RZIJ  = RZIJ - ANINT ( RZIJ )
                 RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                 IF ( RIJSQ .LT. RCUTSQ ) THEN

                    SR2   = SIGSQ / RIJSQ
                    SR6   = SR2 * SR2 * SR2
                    SR12  = SR6 ** 2
                    VIJ   = SR12 - SR6
                    V     = V + VIJ
                    WIJ   = VIJ + SR12
                    W     = W + WIJ
                    FIJ   = WIJ / RIJSQ
                    FXIJ  = FIJ * RXIJ
                    FYIJ  = FIJ * RYIJ
                    FZIJ  = FIJ * RZIJ
                    FXI   = FXI + FXIJ
                    FYI   = FYI + FYIJ
                    FZI   = FZI + FZIJ
                    FX(J) = FX(J) - FXIJ
                    FY(J) = FY(J) - FYIJ
                    FZ(J) = FZ(J) - FZIJ
                    WXY   = WXY + RXIJ * FYIJ

                 ENDIF

                 J = LIST(J)
                 GOTO 2000

              ENDIF

C          ** LOOP OVER NEIGHBOURING CELLS **

              JCELL0 = 16 * ( ICELL - 1 )

              DO 4000 NABOR = 1, 16

                 JCELL = MAP( JCELL0 + NABOR )

                 IF ( JCELL .GT. 0 ) THEN

C                ** LOOP OVER ALL MOLECULES IN NEIGHBOURING CELLS **

                    J = HEAD(JCELL)

3000                IF ( J .NE. 0 ) THEN

                       RXIJ  = RXI - RX(J)
                       RYIJ  = RYI - RY(J)
                       RZIJ  = RZI - RZ(J)
                       CORY  = ANINT ( RYIJ )
                       RXIJ  = RXIJ - CORY * STRAIN
                       RXIJ  = RXIJ - ANINT ( RXIJ )
                       RYIJ  = RYIJ - CORY
                       RZIJ  = RZIJ - ANINT ( RZIJ )
                       RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                       IF ( RIJSQ. LT. RCUTSQ ) THEN

                          SR2   = SIGSQ / RIJSQ
                          SR6   = SR2 * SR2 * SR2
                          SR12  = SR6 ** 2
                          VIJ   = SR12 - SR6
                          V     = V + VIJ
                          WIJ   = VIJ + SR12
                          W     = W + WIJ
                          FIJ   = WIJ / RIJSQ
                          FXIJ  = FIJ * RXIJ
                          FYIJ  = FIJ * RYIJ
                          FZIJ  = FIJ * RZIJ
                          FXI   = FXI + FXIJ
                          FYI   = FYI + FYIJ
                          FZI   = FZI + FZIJ
                          FX(J) = FX(J) - FXIJ
                          FY(J) = FY(J) - FYIJ
                          FZ(J) = FZ(J) - FZIJ
                          WXY   = WXY + RXIJ * FYIJ

                       ENDIF

                       J = LIST(J)
                       GOTO 3000

                    ENDIF

                 ENDIF

4000          CONTINUE

              FX(I) = FXI
              FY(I) = FYI
              FZ(I) = FZI

              I = LIST(I)
              GOTO 1000

           ENDIF

5000    CONTINUE

C    ** INCORPORATE ENERGY FACTORS **

        DO 6000 I = 1, N

           FX(I) = FX(I) * 24.0
           FY(I) = FY(I) * 24.0
           FZ(I) = FZ(I) * 24.0

6000    CONTINUE

        V   = V   * 4.0
        W   = W   * 24.0 / 3.0
        WXY = WXY * 24.0

        RETURN
        END

