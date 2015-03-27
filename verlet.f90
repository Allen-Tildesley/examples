! verlet.f90
SUBROUTINE FORCE ( RCUT, RLIST, SIGMA, UPDATE, V, W )
! Placeholder for f90 routine

************************************************************************************************************************
** FICHE F.19.  THE VERLET NEIGHBOUR LIST                                                                             **
** This FORTRAN code is intended to illustrate points made in the text.  To our knowledge it works correctly.         **
** However it is the responsibility of the user to test it, if it is to be used in a research application.            **
************************************************************************************************************************

**   *******************************************************************
**   ** Modified by MPA 22/5/91 correcting error in force subroutine  **
C    *******************************************************************
C    ** FORCE ROUTINE USING A VERLET NEIGHBOUR LIST.                  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER  N                 NUMBER OF ATOMS                    **
C    ** REAL     RCUT              CUTOFF DISTANCE FOR THE FORCE      **
C    ** REAL     RLIST             OUTER RADIUS OF THE VERLET LIST    **
C    ** REAL     SIGMA             LENNARD JONES LENGTH PARAMETER     **
C    ** REAL     V                 POTENTIAL ENERGY                   **
C    ** REAL     W                 VIRIAL                             **
C    ** REAL     RX(N),RY(N),RZ(N) ATOM POSITIONS                     **
C    ** REAL     FX(N),FY(N),FZ(N) FORCE ON AN ATOM                   **
C    ** LOGICAL  UPDATE            IF TRUE THE LIST IS UPDATED        **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE FORCE ( RCUT, RLIST, SIGMA, UPDATE, V, W )         **
C    **    CALCULATE FORCES USING THE VERLET LIST, UPDATES THE LIST   **
C    ** SUBROUTINE CHECK ( RCUT, RLIST, UPDATE )                      **
C    **    SETS UPDATE TO TRUE WHEN THE LIST NEEDS TO BE UPDATED      **
C    ** SUBROUTINE SAVE                                               **
C    **    SAVES CURRENT CONFIGURATION FOR FUTURE CHECKING.           **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** AT THE START OF A RUN, FORCE IS CALLED WITH UPDATE SET TRUE.  **
C    ** THIS SETS UP THE INITIAL VERLET LIST AND SAVES THE POSITIONS. **
C    ** THEREAFTER CHECK IS CALLED JUST BEFORE EACH CALL OF FORCE TO  **
C    ** DECIDE WHETHER OR NOT A LIST UPDATE IS NECESSARY.             **
C    ** THESE ROUTINES CAN BE USED IN A CONVENTIONAL MD PROGRAM AND   **
C    ** CAN BE EASILY ADAPTED FOR USE IN AN MC PROGRAM. FORCES IS     **
C    ** SPECIFIC TO A FLUID OF LENNARD JONES ATOMS.                   **
C    *******************************************************************


        COMMON / BLOCK1 / RX, RY, RZ, FX, FY, FZ
        COMMON / BLOCK2 / POINT, LIST

C    *******************************************************************
C    ** CALCULATES THE FORCE ON AN ATOM IN A LENNARD-JONES LIQUID     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL     RCUT              CUTOFF DISTANCE FOR THE FORCE      **
C    ** REAL     RLIST             OUTER RADIUS OF THE VERLET LIST    **
C    ** REAL     SIGMA             LENNARD JONES LENGTH PARAMETER     **
C    ** REAL     V                 POTENTIAL ENERGY                   **
C    ** REAL     W                 VIRIAL                             **
C    ** REAL     RX(N),RY(N),RZ(N) ATOM POSITIONS                     **
C    ** REAL     FX(N),FY(N),FZ(N) FORCE ON AN ATOM                   **
C    ** INTEGER  POINT(N)          INDEX TO THE NEIGHBOUR LIST        **
C    ** INTEGER  LIST(MAXNAB)      VERLET NEIGHBOUR LIST              **
C    ** LOGICAL  UPDATE            IF TRUE THE LIST IS UPDATED        **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** FORCE IS CALLED IN TWO MODES. IF UPDATE IS TRUE, THE VERLET   **
C    ** LIST IS UPDATED AND THE FORCES CALCULATED AT THE SAME TIME.   **
C    ** THE CURRENT CONFIGURATION IS SAVED FOR FUTURE CHECKING.       **
C    ** IF UPDATE IS FALSE, THE VERLET LIST IS USED TO FIND THE       **
C    ** NEIGHBOURS OF A GIVEN ATOM AND CALCULATE THE FORCES.          **
C    ** THE ATOMS ARE IN A BOX OF UNIT LENGTH, CENTRED AT THE ORIGIN. **
C    *******************************************************************

        INTEGER     N, MAXNAB
        PARAMETER ( N = 108, MAXNAB = N * 25 )

        REAL        RX(N), RY(N), RZ(N), FX(N), FY(N), FZ(N)
        REAL        RLIST, RCUT, SIGMA, V, W
        INTEGER     POINT(N), LIST(MAXNAB)
        LOGICAL     UPDATE

        REAL        RXI, RYI, RZI, FXI, FYI, FZI
        REAL        RIJSQ, SR2, SR6, VIJ, WIJ, FIJ
        REAL        SIGSQ, RCUTSQ, RLSTSQ
        REAL        RXIJ, RYIJ, RZIJ, FXIJ, FYIJ, FZIJ
        INTEGER     I, J, NLIST
        INTEGER     JBEG, JEND, JNAB

C    *******************************************************************

        RLSTSQ = RLIST * RLIST
        SIGSQ  = SIGMA * SIGMA
        RCUTSQ = RCUT * RCUT

C    ** ZERO FORCES **

        DO 10 I = 1, N

           FX(I) = 0.0
           FY(I) = 0.0
           FZ(I) = 0.0

10      CONTINUE

        V = 0.0
        W = 0.0

        IF ( UPDATE ) THEN

C       ** SAVE CURRENT CONFIGURATION, CONSTRUCT **
C       ** NEIGHBOUR LIST AND CALCULATE FORCES   **

           CALL SAVE

           NLIST = 0

           DO 100 I = 1, N - 1

              POINT(I) = NLIST + 1

              RXI      = RX(I)
              RYI      = RY(I)
              RZI      = RZ(I)
              FXI      = FX(I)
              FYI      = FY(I)
              FZI      = FZ(I)

              DO 99 J = I + 1, N

                 RXIJ  = RXI - RX(J)
                 RYIJ  = RYI - RY(J)
                 RZIJ  = RZI - RZ(J)

                 RXIJ  = RXIJ - ANINT ( RXIJ )
                 RYIJ  = RYIJ - ANINT ( RYIJ )
                 RZIJ  = RZIJ - ANINT ( RZIJ )
                 RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                 IF ( RIJSQ .LT. RLSTSQ ) THEN

                    NLIST = NLIST + 1
                    LIST(NLIST) = J

C                ** REMOVE THIS CHECK IF MAXNAB IS APPROPRIATE **

                    IF ( NLIST .EQ. MAXNAB ) STOP 'LIST TOO SMALL'

                    IF ( RIJSQ .LT. RCUTSQ ) THEN

                       SR2   = SIGSQ / RIJSQ
                       SR6   = SR2 * SR2 * SR2
                       VIJ   = SR6 * ( SR6 - 1.0 )
                       WIJ   = SR6 * ( SR6 - 0.5 )
                       V     = V + VIJ
                       W     = W + WIJ
                       FIJ   = WIJ / RIJSQ
                       FXIJ  = RXIJ * FIJ
                       FYIJ  = RYIJ * FIJ
                       FZIJ  = RZIJ * FIJ
                       FXI   = FXI + FXIJ
                       FYI   = FYI + FYIJ
                       FZI   = FZI + FZIJ
                       FX(J) = FX(J) - FXIJ
                       FY(J) = FY(J) - FYIJ
                       FZ(J) = FZ(J) - FZIJ

                    ENDIF

                 ENDIF

99            CONTINUE

              FX(I) = FXI
              FY(I) = FYI
              FZ(I) = FZI

100        CONTINUE

           POINT(N) = NLIST + 1

        ELSE

C       ** USE THE LIST TO FIND THE NEIGHBOURS **

           DO 200 I = 1, N - 1

C          ** ERROR CORRECTED HERE 22/5/91 BY MPA       **
C          ** THANKS TO WERNER LOOSE FOR SPOTTING IT    **
C          ** AMAZING THAT NO-ONE POINTED IT OUT BEFORE **
C          ** ERROR **  JBEG = LIST(I)
C          ** ERROR **  JEND = LIST(I+1) - 1
              JBEG = POINT(I)
              JEND = POINT(I+1) - 1

C          ** CHECK THAT ATOM I HAS NEIGHBOURS **

              IF( JBEG .LE. JEND ) THEN

                 RXI  = RX(I)
                 RYI  = RY(I)
                 RZI  = RZ(I)
                 FXI  = FX(I)
                 FYI  = FY(I)
                 FZI  = FZ(I)

                 DO 199 JNAB = JBEG, JEND

                    J     = LIST(JNAB)

                    RXIJ  = RXI - RX(J)
                    RYIJ  = RYI - RY(J)
                    RZIJ  = RZI - RZ(J)

                    RXIJ  = RXIJ - ANINT( RXIJ )
                    RYIJ  = RYIJ - ANINT( RYIJ )
                    RZIJ  = RZIJ - ANINT( RZIJ )
                    RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                    IF ( RIJSQ .LT. RCUTSQ ) THEN

                       SR2   = SIGSQ / RIJSQ
                       SR6   = SR2 * SR2 * SR2
                       VIJ   = SR6 * ( SR6 - 1.0 )
                       WIJ   = SR6 * ( SR6 - 0.5 )
                       V     = V + VIJ
                       W     = W + WIJ
                       FIJ   = WIJ / RIJSQ
                       FXIJ  = RXIJ * FIJ
                       FYIJ  = RYIJ * FIJ
                       FZIJ  = RZIJ * FIJ
                       FXI   = FXI + FXIJ
                       FYI   = FYI + FYIJ
                       FZI   = FZI + FZIJ
                       FX(J) = FX(J) - FXIJ
                       FY(J) = FY(J) - FYIJ
                       FZ(J) = FZ(J) - FZIJ

                    ENDIF

199              CONTINUE

                 FX(I) = FXI
                 FY(I) = FYI
                 FZ(I) = FZI

              ENDIF

200        CONTINUE

        ENDIF

        V =  4.0 * V
        W = 48.0 * W / 3.0

        DO 300 I = 1, N

           FX(I) = 48.0 * FX(I)
           FY(I) = 48.0 * FY(I)
           FZ(I) = 48.0 * FZ(I)

300     CONTINUE

        RETURN
        END



        SUBROUTINE CHECK ( RCUT, RLIST, UPDATE )

        COMMON / BLOCK1 / RX, RY, RZ, FX, FY, FZ
        COMMON / BLOCK3 / RX0, RY0, RZ0

C    *******************************************************************
C    ** DECIDES WHETHER THE LIST NEEDS TO BE RECONSTRUCTED.           **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL     RX(N),RY(N),RZ(N)     ATOM POSITIONS                 **
C    ** REAL     RX0(N),RY0(N),RZ0(N)  COORDINATES AT LAST UPDATE     **
C    ** REAL     RLIST                 RADIUS OF VERLET LIST          **
C    ** REAL     RCUT                  CUTOFF DISTANCE FOR FORCES     **
C    ** REAL     DISPMX                LARGEST DISPLACEMENT           **
C    ** INTEGER  N                     NUMBER OF ATOMS                **
C    ** LOGICAL  UPDATE                IF TRUE THE LIST IS UPDATED    **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** CHECK IS CALLED TO SET UPDATE BEFORE EVERY CALL TO FORCE.     **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N), FX(N), FY(N), FZ(N)
        REAL        RX0(N), RY0(N), RZ0(N)
        REAL        RCUT, RLIST
        LOGICAL     UPDATE

        REAL        DISPMX
        INTEGER     I

C    *******************************************************************

C    ** CALCULATE MAXIMUM DISPLACEMENT SINCE LAST UPDATE **

        DISPMX = 0.0

        DO 30 I = 1, N

           DISPMX = MAX ( ABS ( RX(I) - RX0(I) ), DISPMX )
           DISPMX = MAX ( ABS ( RY(I) - RY0(I) ), DISPMX )
           DISPMX = MAX ( ABS ( RZ(I) - RZ0(I) ), DISPMX )

30      CONTINUE

C    ** A CONSERVATIVE TEST OF THE LIST SKIN CROSSING **

        DISPMX = 2.0 * SQRT  ( 3.0 * DISPMX ** 2 )

        UPDATE = ( DISPMX .GT. ( RLIST - RCUT ) )

        RETURN
        END



        SUBROUTINE SAVE

        COMMON / BLOCK1 / RX, RY, RZ, FX, FY, FZ
        COMMON / BLOCK3 / RX0, RY0, RZ0

C    *******************************************************************
C    ** SAVES CURRENT CONFIGURATION FOR FUTURE CHECKING.              **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL     RX(N),RY(N),RZ(N)     ATOM POSITIONS                 **
C    ** REAL     RX0(N),RY0(N),RZ0(N)  STORAGE LOCATIONS FOR SAVE     **
C    ** INTEGER  N                     NUMBER OF ATOMS                **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** SAVE IS CALLED WHENEVER THE NEW VERLET LIST IS CONSTRUCTED.   **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N), FX(N), FY(N), FZ(N)
        REAL        RX0(N), RY0(N), RZ0(N)

        INTEGER     I

C    *******************************************************************

        DO 100 I = 1, N

           RX0(I) = RX(I)
           RY0(I) = RY(I)
           RZ0(I) = RZ(I)

100     CONTINUE

        RETURN
        END
