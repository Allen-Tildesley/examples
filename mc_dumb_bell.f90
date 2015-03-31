! mc_dumb_bell.f90
PROGRAM mc_dumb_bell

        COMMON / block1 / rx, ry, rz, ex, ey, ez

c    *******************************************************************
c    ** constant-nvt monte carlo PROGRAM for hard dumb-bells.         **
c    **                                                               **
c    ** the box is of unit length, -0.5 to +0.5. there are no lookup  **
c    ** tables included.                                              **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                   number of molecules               **
c    ** INTEGER natom               number of atoms per molecule      **
c    ** INTEGER nstep               maximum number of cycles          **
c    ** INTEGER iprint              PRINT interval                    **
c    ** integr  isave               SAVE interval                     **
c    ** INTEGER iratio              max displacement update interval  **
c    ** REAL    rx(n),ry(n),rz(n)   positions                         **
c    ** REAL    ex(n),ey(n),ez(n)   orientations, unit axial vector   **
c    ** REAL    dab(natom)          distance from com to natom        **
c    ** REAL    d                   reduced bond length (d/sigma)     **
c    ** REAL    dens                reduced density                   **
c    ** REAL    sigma               hard sphere diameter              **
c    ** REAL    drmax               reduced maximum displacement      **
c    ** REAL    dotmin              controls angular displacement     **
c    ** LOGICAL ovrlap              true IF dumbbells overlap         **
c    **                                                               **
c    ** routines referenced:                                          **
c    **                                                               **
c    ** SUBROUTINE check ( sigma, dab, ovrlap )                       **
c    **    checks for overlaps in a fluid of hard dumbells            **
c    ** SUBROUTINE orien ( exiold, eyiold, eziold, dotmin, exinew,    **
c    **    :               eyinew, ezinew )                           **
c    **    produces a trial random orientation for a molecule         **
c    ** REAL FUNCTION ranf( dummy )                                   **
c    **    returns a uniform random number between zero and one       **
c    ** SUBROUTINE readcn ( cnfile )                                  **
c    **    reads in a configuration                                   **
c    ** SUBROUTINE test ( rxi, ryi, rzi, i, exi, eyi, ezi, sigma,     **
c    **    :                dab, ovrlap )                             **
c    **    checks for overlaps after the displacement of molecule i   **
c    ** SUBROUTINE writcn ( cnfile )                                  **
c    **    writes out a configuration                                 **
c    *******************************************************************

        INTEGER     n, natom
        PARAMETER ( n = 108, natom = 2 )

        REAL        rx(n), ry(n), rz(n), ex(n), ey(n), ez(n)
        REAL        dab(natom), drmax, dotmin, dens, d, sigma, ratio
        REAL        rxiold, ryiold, rziold, rxinew, ryinew, rzinew
        REAL        exiold, eyiold, eziold, exinew, eyinew, ezinew
        REAL        ranf, dummy, acm, acmmva
        INTEGER     step, i, nstep, iratio, iprint, isave
        LOGICAL     ovrlap
        CHARACTER   title*80, cnfile*80

c    *******************************************************************

c    ** READ input DATA **

        WRITE(*,'(1h1,'' **** program mcbell ****                 '')')
        WRITE(*,'(/   '' constant-nvt monte carlo                 '')')
        WRITE(*,'(    '' for hard dumbells                       ''/)')
        WRITE(*,'('' enter the run title                          '')')
        READ (*,'(a)') title
        WRITE(*,'('' enter number of cycles                       '')')
        READ (*,*) nstep
        WRITE(*,'('' enter number of cycles between output        '')')
        READ (*,*) iprint
        WRITE(*,'('' enter number of cycles between data saves    '')')
        READ (*,*) isave
        WRITE(*,'('' enter interval for update of max. displ.     '')')
        READ (*,*) iratio
        WRITE(*,'('' enter the configuration file name            '')')
        READ (*,'(a)') cnfile
        WRITE(*,'(/'' enter the following in lennard-jones units '',/)')
        WRITE(*,'('' enter the density                            '')')
        READ (*,*) dens
        WRITE(*,'('' enter the maximum displacement               '')')
        READ (*,*) drmax
        WRITE(*,'('' enter the reduced bond length                '')')
        READ (*,*) d

c    ** WRITE input DATA **

        WRITE(*,'(       //1x                    ,a     )') title
        WRITE(*,'('' number of atoms           '',i10   )') n
        WRITE(*,'('' number of cycles          '',i10   )') nstep
        WRITE(*,'('' output frequency          '',i10   )') iprint
        WRITE(*,'('' save frequency            '',i10   )') isave
        WRITE(*,'('' ratio update frequency    '',i10   )') iratio
        WRITE(*,'('' configuration file  name  '',a     )') cnfile
        WRITE(*,'('' density                   '',f10.5 )') dens
        WRITE(*,'('' max. displacement         '',f10.5 )') drmax
        WRITE(*,'('' bond length               '',f10.5 )') d

c    ** set dependent variables **

        sigma  = ( dens / REAL ( n ) ) ** ( 1.0 / 3.0 )
        dab(1) = d * sigma / 2.0
        dab(2) = - dab(1)
        drmax  = drmax * sigma
        dotmin = 0.2

c    ** WRITE out some useful information **

        WRITE( *, '( '' number of molecules   =  '', i10   )' )  n
        WRITE( *, '( '' number of atoms       =  '', i10   )' )  natom
        WRITE( *, '( '' sigma  / box          =  '', f10.5 )' )  sigma
        WRITE( *, '( '' dab(1) / box          =  '', f10.5 )' )  dab(1)
        WRITE( *, '( '' dab(2) / box          =  '', f10.5 )' )  dab(2)
        WRITE( *, '( '' drmax  / box          =  '', f10.5 )' )  drmax
        WRITE( *, '( '' dotmin                =  '', f10.5 )' )  dotmin

c    ** READ in initial configuration **

        CALL readcn ( cnfile )

c    ** check for overlaps in initial configuration **

        CALL check ( sigma, dab, ovrlap )

        IF ( ovrlap ) STOP 'overlap in initial configuration'

c    ** zero accumulators **

        acm    = 0.0
        acmmva = 0.0

        WRITE( *, '(//'' start of markov chain                ''//)')
        WRITE( *, '( ''     acm    ratio     drmax     dotmin  '')')

c    *******************************************************************
c    ** loop over cycles begins                                       **
c    *******************************************************************

        DO 100 step = 1, nstep

c       ** loop over molecules **

           DO 99 i = 1, n

              rxiold = rx(i)
              ryiold = ry(i)
              rziold = rz(i)
              exiold = ex(i)
              eyiold = ey(i)
              eziold = ez(i)

c          ** move i and pickup the central image **

              rxinew = rxiold + ( 2.0 * ranf ( dummy ) - 1.0 ) * drmax
              ryinew = ryiold + ( 2.0 * ranf ( dummy ) - 1.0 ) * drmax
              rzinew = rziold + ( 2.0 * ranf ( dummy ) - 1.0 ) * drmax

              rxinew = rxinew - ANINT ( rxinew )
              ryinew = ryinew - ANINT ( ryinew )
              rzinew = rzinew - ANINT ( rzinew )

c          ** change the orientation of molecule i **

              CALL orien ( exiold, eyiold, eziold, dotmin,
     :                     exinew, eyinew, ezinew         )

c          ** check for acceptance **

              CALL test ( rxinew, ryinew, rzinew, i,
     :                    exinew, eyinew, ezinew, sigma, dab, ovrlap )

              IF ( .NOT. ovrlap ) THEN

c             ** accept move **

                 rx(i) = rxinew
                 ry(i) = ryinew
                 rz(i) = rzinew
                 ex(i) = exinew
                 ey(i) = eyinew
                 ez(i) = ezinew
                 acmmva = acmmva + 1.0

              ENDIF

              acm = acm + 1.0

99         CONTINUE

c       ****************************************************************
c       ** loop over molecules complete                               **
c       ****************************************************************

c       ** perform periodic operations  **

c       ** change maximum displacement **

           IF ( MOD ( step, iratio ) .EQ. 0 ) THEN

              ratio = acmmva / REAL ( n * iratio )

              IF ( ratio .GT. 0.5 ) THEN

                 drmax = drmax * 1.05
                 dotmin = dotmin * 1.025

              ELSE

                 drmax = drmax * 0.95
                 dotmin = dotmin * 0.975

              ENDIF

              acmmva = 0.0

           ENDIF

c       ** WRITE out runtime information **

           IF ( MOD ( step, iprint ) .EQ. 0 ) THEN

              WRITE(*,'(i8,3f10.4)') INT(acm), ratio, drmax, dotmin

           ENDIF

c       ** WRITE out the configuration at intervals **

           IF ( MOD ( step, isave ) .EQ. 0 ) THEN

              CALL writcn ( cnfile )

              CALL check ( sigma, dab, ovrlap )

              IF ( ovrlap ) STOP 'overlap during the run'

           ENDIF

100     CONTINUE

c    *******************************************************************
c    ** ends the loop over cycles                                     **
c    *******************************************************************

c    ** checks for ovrlaps in the FINAL configuration  **

        CALL check ( sigma, dab, ovrlap )

        IF ( ovrlap ) STOP 'overlap in final configuration'

c    ** WRITE out the FINAL configuration from the run **

        CALL writcn ( cnfile )

        STOP
        END



        SUBROUTINE orien ( exiold, eyiold, eziold, dotmin,
     :                     exinew, eyinew, ezinew         )

c    *******************************************************************
c    ** finds a trial random orientation of a linear molecule.        **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** REAL      exiold,eyiold,eziold  old axial vector for i        **
c    ** REAL      eyinew,eyinew,ezinew  NEW axial vector for i        **
c    ** REAL      dot                   dot product of old and NEW    **
c    **                                 axial vectors                 **
c    ** REAL      dotmin                minimum allowed dot product   **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** the method USE a rejection technique to create a trial        **
c    ** orientation of molecule i subject to the constraint that      **
c    ** the cosine of the angle between the old and NEW axial         **
c    ** vectors, dot, is greater than ( 1.0 - dotmin ).               **
c    *******************************************************************

        REAL    exiold, eyiold, eziold, exinew, eyinew, ezinew, dotmin

        REAL    dot, xi1, xi2, xi, xisq
        REAL    ranf, dummy

c    *******************************************************************

c    ** initialise dot **

        dot  = 0.0

c    ** iterative loop **

1000    IF ( ( 1.0 - dot ) .GE. dotmin ) THEN

c       ** initialise xisq **

           xisq = 1.0

c       ** inner iterative loop **

2000       IF ( xisq .GE. 1.0 ) THEN

              xi1  = ranf ( dummy ) * 2.0 - 1.0
              xi2  = ranf ( dummy ) * 2.0 - 1.0
              xisq = xi1 * xi1 + xi2 * xi2

              GOTO 2000

           ENDIF

           xi     = SQRT ( 1.0 - xisq )
           exinew = 2.0 * xi1 * xi
           eyinew = 2.0 * xi2 * xi
           ezinew = 1.0 - 2.0 * xisq
           dot    = exinew * exiold + eyinew * eyiold + ezinew * eziold

           GOTO 1000

        ENDIF

        RETURN
        END



        SUBROUTINE test ( rxi, ryi, rzi, i, exi, eyi, ezi, sigma,
     :                    dab, ovrlap )

        COMMON / block1 / rx, ry,rz, ex, ey, ez

c    *******************************************************************
c    ** checks for overlap of i WITH all other molecules.             **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER i                  the molecule of interest           **
c    ** INTEGER n                  number of molecules                **
c    ** INTEGER natom              number of atoms per molecule       **
c    ** REAL    rxi,ryi,rzi        position of molecule i             **
c    ** REAL    exi,eyi,ezi,       orientation of molecule i          **
c    ** REAL    rx(n),ry(n),rz(n)  molecular positions                **
c    ** REAL    ex(n),ey(n),ez(n)  molecular orientations             **
c    ** REAL    dab(natom)         position of atoms in a molecule    **
c    ** REAL    sigma              reduced atom diameter              **
c    ** LOGICAL ovrlap             true IF molecule i overlaps        **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** called after a trial displacement of molecule i to establish  **
c    ** whether there is an overlap in the trial configuration.       **
c    *******************************************************************

        INTEGER     natom, n
        PARAMETER ( natom = 2, n = 108 )

        REAL        rx(n), ry(n), rz(n)
        REAL        ex(n), ey(n), ez(n)
        REAL        rxi, ryi, rzi, exi, eyi, ezi
        REAL        sigma, dab(natom)
        INTEGER     i
        LOGICAL     ovrlap

        REAL        rxij, ryij, rzij, exj, eyj, ezj
        REAL        rxab, ryab, rzab, dabi, sigsq, rabsq
        INTEGER     j, ia, jb

c    *******************************************************************

        ovrlap = .FALSE.
        sigsq  = sigma * sigma

c    ** loops over molecules except i **

        DO 100 j = 1, n

           IF ( j .NE. i ) THEN

              exj = ex(j)
              eyj = ey(j)
              ezj = ez(j)

              rxij = rxi - rx(j)
              ryij = ryi - ry(j)
              rzij = rzi - rz(j)

              rxij = rxij - ANINT ( rxij )
              ryij = ryij - ANINT ( ryij )
              rzij = rzij - ANINT ( rzij )

c          ** loops over atoms **

              DO 99 ia = 1, natom

                 dabi = dab(ia)

                 DO 98 jb = 1, natom

                    rxab = rxij + exi * dabi + exj * dab(jb)
                    ryab = ryij + eyi * dabi + eyj * dab(jb)
                    rzab = rzij + ezi * dabi + ezj * dab(jb)

                    rabsq = rxab * rxab + ryab * ryab + rzab * rzab

                    IF ( rabsq .LT. sigsq ) THEN

                       ovrlap = .TRUE.
                       RETURN

                    ENDIF

98               CONTINUE

99            CONTINUE

           ENDIF

100     CONTINUE

        RETURN
        END



        SUBROUTINE check ( sigma, dab, ovrlap )

        COMMON / block1 / rx, ry, rz, ex, ey, ez

c    *******************************************************************
c    ** routine to check for overlaps in a fluid of hard dumbbells    **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                 number of molecules                 **
c    ** INTEGER natom             number of atoms per molecule        **
c    ** REAL    rx(n),ry(n),rz(n) molecular positions                 **
c    ** REAL    ex(n),ey(n),ez(n) molecular orientations              **
c    ** REAL    dab(natom)        position of atoms in a molecule     **
c    ** REAL    sigma             reduced atom diameter               **
c    ** LOGICAL ovrlap            true IF two dumbells overlap        **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** called at intervals during the run to check for overlaps. IF  **
c    ** ovrlap is returned WITH a true VALUE THEN there is an error   **
c    ** in the PROGRAM and the execution is stopped.                  **
c    *******************************************************************

        INTEGER     natom, n
        PARAMETER ( natom = 2, n = 108 )

        REAL        rx(n), ry(n), rz(n)
        REAL        ex(n), ey(n), ez(n)
        REAL        sigma, dab(natom)
        LOGICAL     ovrlap

        REAL        rxi, ryi, rzi, rxij, ryij, rzij, exi, eyi, ezi
        REAL        exj, eyj, ezj, rxab, ryab, rzab, dabi, sigsq, rabsq
        INTEGER     i, j, ia, jb

c    *******************************************************************

        ovrlap = .FALSE.
        sigsq  = sigma * sigma

c    ** loops over molecules **

        DO 100 i = 1, n - 1

           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)
           exi = ex(i)
           eyi = ey(i)
           ezi = ez(i)

           DO 99 j = i + 1, n

              rxij = rxi - rx(j)
              ryij = ryi - ry(j)
              rzij = rzi - rz(j)

              rxij = rxij - ANINT ( rxij )
              ryij = ryij - ANINT ( ryij )
              rzij = rzij - ANINT ( rzij )

              exj  = ex(j)
              eyj  = ey(j)
              ezj  = ez(j)

c          ** loops over atoms **

              DO 98 ia = 1, natom

                 dabi = dab(ia)

                 DO 97 jb = 1, natom

                    rxab = rxij + exi * dabi + exj * dab(jb)
                    ryab = ryij + eyi * dabi + eyj * dab(jb)
                    rzab = rzij + ezi * dabi + ezj * dab(jb)

                    rabsq = rxab * rxab + ryab * ryab + rzab * rzab

                    IF ( rabsq .LT. sigsq ) THEN

                       ovrlap = .TRUE.
                       RETURN

                    ENDIF

97               CONTINUE

98            CONTINUE

99         CONTINUE

100     CONTINUE

        RETURN
        END



        SUBROUTINE readcn ( cnfile )

        COMMON / block1 / rx, ry, rz, ex, ey, ez

c    *******************************************************************
c    ** SUBROUTINE to READ in the configuration from unit 10          **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )
        CHARACTER   cnfile*(*)
        REAL        rx(n), ry(n), rz(n), ex(n), ey(n), ez(n)

        INTEGER     cnunit
        PARAMETER ( cnunit = 10 )
        INTEGER     nn

c   ********************************************************************

        OPEN ( unit = cnunit, file = cnfile, status = 'old',
     :         form = 'unformatted'                        )

        READ ( cnunit ) nn
        IF ( nn .NE. n ) STOP 'n error in readcn'
        READ ( cnunit ) rx, ry, rz
        READ ( cnunit ) ex, ey, ez

        CLOSE ( unit = cnunit )

        RETURN
        END



        SUBROUTINE writcn ( cnfile )

        COMMON / block1 / rx, ry, rz, ex, ey, ez

c    *******************************************************************
c    ** SUBROUTINE to WRITE out the configuration to unit 10          **
c    *******************************************************************

        INTEGER      n
        PARAMETER (  n = 108 )
        CHARACTER    cnfile*(*)
        REAL         rx(n), ry(n), rz(n), ex(n), ey(n), ez(n)

        INTEGER      cnunit
        PARAMETER (  cnunit = 10 )

c   ********************************************************************

        OPEN ( unit = cnunit, file = cnfile, status = 'unknown',
     :         form = 'unformatted'                        )

        WRITE ( cnunit ) n
        WRITE ( cnunit ) rx, ry, rz
        WRITE ( cnunit ) ex, ey, ez

        CLOSE ( unit = cnunit )

        RETURN
        END



        REAL FUNCTION ranf ( dummy )

c    *******************************************************************
c    ** returns a uniform random variate in the range 0 to 1.         **
c    **                                                               **
c    **                 ***************                               **
c    **                 **  warning  **                               **
c    **                 ***************                               **
c    **                                                               **
c    ** good random number generators are machine specific.           **
c    ** please USE the one recommended for your machine.              **
c    *******************************************************************

        INTEGER     l, c, m
        PARAMETER ( l = 1029, c = 221591, m = 1048576 )

        INTEGER     seed
        REAL        dummy
        SAVE        seed
        DATA        seed / 0 /

c    *******************************************************************

        seed = MOD ( seed * l + c, m )
        ranf = REAL ( seed ) / m

        RETURN
        END



