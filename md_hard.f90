! md_hard.f90
PROGRAM md_hard
! Placeholder for f90 program

********************************************************************************
** FICHE F.10.  HARD SPHERE MOLECULAR DYNAMICS PROGRAM                        **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************




        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ
        COMMON / BLOCK2 / COLTIM, PARTNR

C    *******************************************************************
C    ** MOLECULAR DYNAMICS OF HARD SPHERE ATOMS.                      **
C    **                                                               **
C    ** THIS PROGRAM TAKES IN A HARD-SPHERE CONFIGURATION (POSITIONS  **
C    ** AND VELOCITIES), CHECKS FOR OVERLAPS, AND THEN CONDUCTS A     **
C    ** MOLECULAR DYNAMICS SIMULATION RUN FOR A SPECIFIED NUMBER OF   **
C    ** COLLISIONS.  THE PROGRAM IS FAIRLY EFFICIENT, BUT USES NO     **
C    ** SPECIAL NEIGHBOUR LISTS, SO IS RESTRICTED TO A SMALL NUMBER   **
C    ** OF PARTICLES (<500).  IT IS ALWAYS ASSUMED THAT COLLISIONS    **
C    ** CAN BE PREDICTED BY LOOKING AT NEAREST NEIGHBOUR PARTICLES IN **
C    ** THE MINIMUM IMAGE CONVENTION OF PERIODIC BOUNDARIES.          **
C    ** THE BOX IS TAKEN TO BE OF UNIT LENGTH.                        **
C    ** HOWEVER, RESULTS ARE GIVEN IN UNITS WHERE SIGMA=1, KT=1.      **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                    NUMBER OF ATOMS                  **
C    ** REAL    RX(N),RY(N),RZ(N)    ATOM POSITIONS                   **
C    ** REAL    VX(N),VY(N),VZ(N)    ATOM VELOCITIES                  **
C    ** REAL    COLTIM(N)            TIME TO NEXT COLLISION           **
C    ** INTEGER PARTNR(N)            COLLISION PARTNER                **
C    ** REAL    SIGMA                ATOM DIAMETER                    **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE READCN ( CNFILE )                                  **
C    **    READS IN CONFIGURATION                                     **
C    ** SUBROUTINE CHECK ( SIGMA, OVRLAP, E )                         **
C    **    CHECKS CONFIGURATION AND CALCULATES ENERGY                 **
C    ** SUBROUTINE UPLIST ( SIGMA, I )                                **
C    **    SEEKS COLLISIONS WITH J>I                                  **
C    ** SUBROUTINE DNLIST ( SIGMA, I )                                **
C    **    SEEKS COLLISIONS WITH J<I                                  **
C    ** SUBROUTINE BUMP ( SIGMA, I, J, W )                            **
C    **    DOES COLLISION DYNAMICS AND CALCULATES COLLISION VIRIAL    **
C    ** SUBROUTINE WRITCN ( CNFILE )                                  **
C    **    WRITES OUT CONFIGURATION                                   **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        TIMBIG
        PARAMETER ( TIMBIG = 1.0E10 )

        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        COLTIM(N)
        INTEGER     PARTNR(N)
        REAL        SIGMA

        INTEGER     I, J, K, NCOLL, COLL
        REAL        DENSTY, DIJ, TIJ, T, RATE
        REAL        E, EN, ENKT, W, PVNKT1, ACW, TEMP, TBC
        CHARACTER   TITLE*80, CNFILE*30
        LOGICAL     OVRLAP

C    *******************************************************************

        WRITE(*,'(1H1,'' **** PROGRAM SPHERE ****               '')')
        WRITE(*,'(//  '' MOLECULAR DYNAMICS OF HARD SPHERES     '')')
        WRITE(*,'(//  '' RESULTS IN UNITS KT = SIGMA = 1        '')')

C    ** READ IN BASIC SIMULATION PARAMETERS **

        WRITE(*,'('' ENTER RUN TITLE                            '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER REDUCED DENSITY (N/V)*SIGMA**3       '')')
        READ (*,*) DENSTY
        WRITE(*,'('' ENTER NUMBER OF COLLISIONS REQUIRED        '')')
        READ (*,*) NCOLL
        WRITE(*,'('' ENTER CONFIGURATION FILENAME               '')')
        READ (*,'(A)') CNFILE

        WRITE(*,'('' RUN TITLE              '',A)'    ) TITLE
        WRITE(*,'('' REDUCED DENSITY IS     '',F15.5)') DENSTY
        WRITE(*,'('' COLLISIONS REQUIRED    '',I15)'  ) NCOLL
        WRITE(*,'('' CONFIGURATION FILENAME '',A)'    ) CNFILE

        SIGMA = ( DENSTY / REAL(N) ) ** ( 1.0 / 3.0 )

C    ** READ IN CONFIGURATION **

        CALL READCN ( CNFILE )

C    ** CHECK FOR PARTICLE OVERLAPS **
C    ** CALCULATE ENERGY            **

        CALL CHECK ( SIGMA, OVRLAP, E )

        IF ( OVRLAP ) STOP 'PARTICLE OVERLAP IN INITIAL CONFIGURATION'

        EN = E / REAL ( N )
        TEMP = 2.0 * EN / 3.0
        WRITE(*,'('' TEMPERATURE          '',F15.5)') TEMP
        ENKT = EN / TEMP
        WRITE(*,'('' INITIAL E/NKT        '',F15.5)') ENKT

C    ** SET UP INITIAL COLLISION LISTS COLTIM AND PARTNR **

        DO 10 I = 1, N

           COLTIM(I) = TIMBIG
           PARTNR(I) = N

10      CONTINUE

        DO 20 I = 1, N

           CALL UPLIST ( SIGMA, I )

20      CONTINUE

C    ** ZERO VIRIAL ACCUMULATOR **

        ACW = 0.0

        WRITE(*,'(//'' **** START OF DYNAMICS **** '')')

C    *******************************************************************
C    ** MAIN LOOP BEGINS                                              **
C    *******************************************************************

        T = 0.0

        DO 1000 COLL = 1, NCOLL

C       ** LOCATE MINIMUM COLLISION TIME **

           TIJ = TIMBIG

           DO 200 K = 1, N

              IF ( COLTIM(K) .LT. TIJ ) THEN

                 TIJ = COLTIM(K)
                 I   = K

              ENDIF

200        CONTINUE

           J = PARTNR(I)

C       ** MOVE PARTICLES FORWARD BY TIME TIJ **
C       ** AND REDUCE COLLISION TIMES         **
C       ** APPLY PERIODIC BOUNDARIES          **

           T = T + TIJ

           DO 300 K = 1, N

              COLTIM(K) = COLTIM(K) - TIJ
              RX(K) = RX(K) + VX(K) * TIJ
              RY(K) = RY(K) + VY(K) * TIJ
              RZ(K) = RZ(K) + VZ(K) * TIJ
              RX(K) = RX(K) - ANINT ( RX(K) )
              RY(K) = RY(K) - ANINT ( RY(K) )
              RZ(K) = RZ(K) - ANINT ( RZ(K) )

300        CONTINUE

C       ** COMPUTE COLLISION DYNAMICS **

           CALL BUMP ( SIGMA, I, J, W )

           ACW = ACW + W

C       ** RESET COLLISION LISTS FOR     **
C       ** THOSE PARTICLES WHICH NEED IT **

           DO 400 K = 1, N

              IF ( ( K .EQ. I ) .OR. ( PARTNR(K) .EQ. I ) .OR.
     :             ( K .EQ. J ) .OR. ( PARTNR(K) .EQ. J )     ) THEN

                 CALL UPLIST ( SIGMA, K )

              ENDIF

400        CONTINUE

           CALL DNLIST ( SIGMA, I )
           CALL DNLIST ( SIGMA, J )

1000    CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS.                                               **
C    *******************************************************************

        WRITE(*,'(//'' **** END OF DYNAMICS **** '')')

        WRITE(*,'(/'' FINAL COLLIDING PAIR '',2I5)') I, J

C    ** CHECK FOR PARTICLE OVERLAPS **

        CALL CHECK ( SIGMA, OVRLAP, E )

        IF ( OVRLAP ) THEN

           WRITE(*,'('' PARTICLE OVERLAP IN FINAL CONFIGURATION '')')

        ENDIF

C    ** WRITE OUT CONFIGURATION **

        CALL WRITCN ( CNFILE )

C     ** WRITE OUT INTERESTING INFORMATION **

        PVNKT1 = ACW / REAL ( N ) / 3.0 / T / TEMP
        EN = E / REAL ( N )
        ENKT = EN / TEMP
        T = T * SQRT ( TEMP ) / SIGMA
        RATE = REAL ( NCOLL ) / T
        TBC  = REAL ( N ) / RATE / 2.0

        WRITE(*,'('' FINAL TIME IS          '',F15.8)') T
        WRITE(*,'('' COLLISION RATE IS      '',F15.8)') RATE
        WRITE(*,'('' MEAN COLLISION TIME    '',F15.8)') TBC
        WRITE(*,'('' FINAL E/NKT IS         '',F15.8)') ENKT
        WRITE(*,'('' PV/NKT - 1 IS          '',F15.8)') PVNKT1

        STOP
        END



        SUBROUTINE CHECK ( SIGMA, OVRLAP, E )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ
        COMMON / BLOCK2 / COLTIM, PARTNR

C    *******************************************************************
C    ** TESTS FOR PAIR OVERLAPS AND CALCULATES KINETIC ENERGY.        **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        COLTIM(N)
        INTEGER     PARTNR(N)

        REAL        SIGMA, E
        LOGICAL     OVRLAP

        INTEGER     I, J
        REAL        RXI, RYI, RZI, RXIJ, RYIJ, RZIJ, RIJSQ, SIGSQ, RIJ
        REAL        TOL
        PARAMETER ( TOL = 1.0E-4 )

C    *******************************************************************

        SIGSQ  = SIGMA ** 2
        OVRLAP = .FALSE.
        E = 0.0

        DO 100 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)

           DO 99 J = I + 1, N

              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RZIJ = RZI - RZ(J)
              RXIJ = RXIJ - ANINT ( RXIJ )
              RYIJ = RYIJ - ANINT ( RYIJ )
              RZIJ = RZIJ - ANINT ( RZIJ )
              RIJSQ = RXIJ ** 2 + RYIJ ** 2 + RZIJ ** 2

              IF ( RIJSQ .LT. SIGSQ ) THEN

                 RIJ = SQRT ( RIJSQ / SIGSQ )
                 WRITE(*,'('' I,J,RIJ/SIGMA = '',2I5,F15.8)')
     :              I, J, RIJ

                 IF ( ( 1.0 - RIJ ) .GT. TOL ) OVRLAP = .TRUE.

              ENDIF

99         CONTINUE

100     CONTINUE

        DO 200 I = 1, N

           E = E + VX(I) ** 2 + VY(I) ** 2 + VZ(I) ** 2

200     CONTINUE

        E = 0.5 * E

        RETURN
        END



        SUBROUTINE UPLIST ( SIGMA, I )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ
        COMMON / BLOCK2 / COLTIM, PARTNR

C    *******************************************************************
C    ** LOOKS FOR COLLISIONS WITH ATOMS J > I                         **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        TIMBIG
        PARAMETER ( TIMBIG = 1.0E10 )

        INTEGER     I
        REAL        SIGMA
        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        COLTIM(N)
        INTEGER     PARTNR(N)

        INTEGER     J
        REAL        RXI, RYI, RZI, RXIJ, RYIJ, RZIJ
        REAL        VXI, VYI, VZI, VXIJ, VYIJ, VZIJ
        REAL        RIJSQ, VIJSQ, BIJ, TIJ, DISCR, SIGSQ

C    *******************************************************************

        IF ( I .EQ. N ) RETURN

        SIGSQ = SIGMA ** 2
        COLTIM(I) = TIMBIG
        RXI = RX(I)
        RYI = RY(I)
        RZI = RZ(I)
        VXI = VX(I)
        VYI = VY(I)
        VZI = VZ(I)

        DO 100 J = I + 1, N

           RXIJ = RXI - RX(J)
           RYIJ = RYI - RY(J)
           RZIJ = RZI - RZ(J)
           RXIJ = RXIJ - ANINT ( RXIJ )
           RYIJ = RYIJ - ANINT ( RYIJ )
           RZIJ = RZIJ - ANINT ( RZIJ )
           VXIJ = VXI - VX(J)
           VYIJ = VYI - VY(J)
           VZIJ = VZI - VZ(J)
           BIJ  = RXIJ * VXIJ + RYIJ * VYIJ + RZIJ * VZIJ

           IF ( BIJ .LT. 0.0 ) THEN

              RIJSQ = RXIJ ** 2 + RYIJ ** 2 + RZIJ ** 2
              VIJSQ = VXIJ ** 2 + VYIJ ** 2 + VZIJ ** 2
              DISCR = BIJ ** 2 - VIJSQ * ( RIJSQ - SIGSQ )

              IF ( DISCR .GT. 0.0 ) THEN

                 TIJ = ( -BIJ - SQRT ( DISCR ) ) / VIJSQ

                 IF ( TIJ .LT. COLTIM(I) ) THEN

                    COLTIM(I) = TIJ
                    PARTNR(I) = J

                 ENDIF

              ENDIF

           ENDIF

100     CONTINUE

        RETURN
        END



        SUBROUTINE DNLIST ( SIGMA, J )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ
        COMMON / BLOCK2 / COLTIM, PARTNR

C    *******************************************************************
C    ** LOOKS FOR COLLISIONS WITH ATOMS I < J                         **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        TIMBIG
        PARAMETER ( TIMBIG = 1.E10 )

        INTEGER     J
        REAL        SIGMA
        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)
        REAL        COLTIM(N)
        INTEGER     PARTNR(N)

        INTEGER     I
        REAL        RXJ, RYJ, RZJ, RXIJ, RYIJ, RZIJ
        REAL        VXJ, VYJ, VZJ, VXIJ, VYIJ, VZIJ
        REAL        RIJSQ, VIJSQ, BIJ, TIJ, DISCR, SIGSQ

C    *******************************************************************

        IF ( J .EQ. 1 ) RETURN

        SIGSQ = SIGMA ** 2
        RXJ = RX(J)
        RYJ = RY(J)
        RZJ = RZ(J)
        VXJ = VX(J)
        VYJ = VY(J)
        VZJ = VZ(J)

        DO 100 I = 1, J - 1

           RXIJ = RX(I) - RXJ
           RYIJ = RY(I) - RYJ
           RZIJ = RZ(I) - RZJ
           RXIJ = RXIJ - ANINT ( RXIJ )
           RYIJ = RYIJ - ANINT ( RYIJ )
           RZIJ = RZIJ - ANINT ( RZIJ )
           VXIJ = VX(I) - VXJ
           VYIJ = VY(I) - VYJ
           VZIJ = VZ(I) - VZJ
           BIJ  = RXIJ * VXIJ + RYIJ * VYIJ + RZIJ * VZIJ

           IF ( BIJ .LT. 0.0 ) THEN

              RIJSQ = RXIJ ** 2 + RYIJ ** 2 + RZIJ ** 2
              VIJSQ = VXIJ ** 2 + VYIJ ** 2 + VZIJ ** 2
              DISCR = BIJ ** 2 - VIJSQ * ( RIJSQ - SIGSQ )

              IF ( DISCR .GT. 0.0 ) THEN

                 TIJ = ( - BIJ - SQRT ( DISCR ) ) / VIJSQ

                 IF ( TIJ .LT. COLTIM(I) ) THEN

                    COLTIM(I) = TIJ
                    PARTNR(I) = J

                 ENDIF

              ENDIF

           ENDIF

100     CONTINUE

        RETURN
        END




        SUBROUTINE BUMP ( SIGMA, I, J, W )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** COMPUTES COLLISION DYNAMICS FOR PARTICLES I AND J.            **
C    **                                                               **
C    ** IT IS ASSUMED THAT I AND J ARE IN CONTACT.                    **
C    ** THE ROUTINE ALSO COMPUTES COLLISIONAL VIRIAL W.               **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     I, J
        REAL        SIGMA, W
        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)

        REAL        RXIJ, RYIJ, RZIJ, FACTOR
        REAL        DELVX, DELVY, DELVZ, SIGSQ

C    *******************************************************************

        SIGSQ = SIGMA ** 2
        RXIJ = RX(I) - RX(J)
        RYIJ = RY(I) - RY(J)
        RZIJ = RZ(I) - RZ(J)
        RXIJ = RXIJ - ANINT ( RXIJ )
        RYIJ = RYIJ - ANINT ( RYIJ )
        RZIJ = RZIJ - ANINT ( RZIJ )

        FACTOR = ( RXIJ * ( VX(I) - VX(J) ) +
     :             RYIJ * ( VY(I) - VY(J) ) +
     :             RZIJ * ( VZ(I) - VZ(J) ) ) / SIGSQ
        DELVX = - FACTOR * RXIJ
        DELVY = - FACTOR * RYIJ
        DELVZ = - FACTOR * RZIJ

        VX(I) = VX(I) + DELVX
        VX(J) = VX(J) - DELVX
        VY(I) = VY(I) + DELVY
        VY(J) = VY(J) - DELVY
        VZ(I) = VZ(I) + DELVZ
        VZ(J) = VZ(J) - DELVZ

        W = DELVX * RXIJ + DELVY * RYIJ + DELVZ * RZIJ

        RETURN
        END



        SUBROUTINE READCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** READS CONFIGURATION FROM CHANNEL 10                           **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )

        INTEGER     NN
        CHARACTER   CNFILE*(*)
        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)

C    *******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS = 'OLD', FORM = 'UNFORMATTED' )

        READ ( CNUNIT ) NN
        IF ( NN .NE. N ) STOP 'N ERROR IN READCN'
        READ ( CNUNIT ) RX, RY, RZ
        READ ( CNUNIT ) VX, VY, VZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ

C    *******************************************************************
C    ** WRITES CONFIGURATION TO CHANNEL 10                            **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )

        CHARACTER   CNFILE*(*)
        REAL        RX(N), RY(N), RZ(N), VX(N), VY(N), VZ(N)

C    *******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS = 'UNKNOWN', FORM = 'UNFORMATTED' )

        WRITE ( CNUNIT ) N
        WRITE ( CNUNIT ) RX, RY, RZ
        WRITE ( CNUNIT ) VX, VY, VZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



