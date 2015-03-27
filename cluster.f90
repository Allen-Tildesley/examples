! cluster.f90
SUBROUTINE cluster ( RCL, IT, NIT )
! Placeholder for f90 subroutine

********************************************************************************
** FICHE F.34.  AN EFFICIENT CLUSTERING ROUTINE                               **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************




        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** ROUTINE TO IDENTIFY ATOM CLUSTERS IN A CONFIGURATION.         **
C    **                                                               **
C    ** THIS ROUTINE SORTS N ATOMS INTO CLUSTERS DEFINED BY A         **
C    ** CRITICAL CLUSTER RADIUS, AND COUNTS THE NUMBER OF ATOMS IN    **
C    ** THE CLUSTER CONTAINING THE ATOM 'IT'.  THE ATOMS ARE IN A     **
C    ** BOX OF UNIT LENGTH CENTRED AT THE ORIGIN                      **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** STODDARD J COMP PHYS, 27, 291, 1977.                          **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF ATOMS                   **
C    ** INTEGER IT                  AN ATOM IN A PARTICULAR CLUSTER   **
C    ** INTEGER L(N)                A LINKED-LIST                     **
C    ** INTEGER NIT                 NUMBER OF ATOMS IN THAT CLUSTER   **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    RCL                 CRITICAL CLUSTER DISTANCE         **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        RCL
        INTEGER     IT, NIT

        REAL        RCLSQ, RXJK, RYJK, RZJK
        REAL        RJKSQ, RXJ, RYJ, RZJ
        INTEGER     I, J, K, LK, LIT, L(N)

C       ****************************************************************

        RCLSQ = RCL * RCL

C    ** SET UP THE SORTING ARRAY **

        DO 10 I = 1, N

           L(I) = I

10      CONTINUE

C    ** SORT THE CLUSTERS **

        DO 50 I = 1, N - 1

           IF ( I .EQ. L(I) ) THEN

              J   = I
              RXJ = RX(J)
              RYJ = RY(J)
              RZJ = RZ(J)

              DO 20 K = I + 1, N

                 LK = L(K)

                 IF ( LK .EQ. K ) THEN

                    RXJK  = RXJ - RX(K)
                    RYJK  = RYJ - RY(K)
                    RZJK  = RZJ - RZ(K)
                    RXJK  = RXJK - ANINT ( RXJK )
                    RYJK  = RYJK - ANINT ( RYJK )
                    RZJK  = RZJK - ANINT ( RZJK )
                    RJKSQ = RXJK * RXJK + RYJK * RYJK + RZJK * RZJK

                    IF ( RJKSQ .LE. RCLSQ ) THEN

                       L(K) = L(J)
                       L(J) = LK

                    ENDIF

                 ENDIF

20            CONTINUE

              J   = L(J)
              RXJ = RX(J)
              RYJ = RY(J)
              RZJ = RZ(J)

30            IF ( J .NE. I ) THEN

                 DO 40 K = I + 1, N

                    LK = L(K)

                    IF ( LK .EQ. K ) THEN

                       RXJK  = RXJ - RX(K)
                       RYJK  = RYJ - RY(K)
                       RZJK  = RZJ - RZ(K)
                       RXJK  = RXJK - ANINT ( RXJK )
                       RYJK  = RYJK - ANINT ( RYJK )
                       RZJK  = RZJK - ANINT ( RZJK )
                       RJKSQ = RXJK * RXJK + RYJK * RYJK + RZJK * RZJK

                       IF ( RJKSQ .LE. RCLSQ ) THEN

                          L(K) = L(J)
                          L(J) = LK

                       ENDIF

                    ENDIF

40               CONTINUE

                 J   = L(J)
                 RXJ = RX(J)
                 RYJ = RY(J)
                 RZJ = RZ(J)

                 GO TO 30

              ENDIF

           ENDIF

50      CONTINUE

C   **  COUNT THE NUMBER IN A CLUSTER CONTAINING ATOM IT **

        NIT = 1
        LIT = L(IT)

60      IF ( LIT .NE. IT ) THEN

           NIT = NIT + 1
           LIT = L(LIT)

           GO TO 60

        ENDIF

        RETURN
        END



