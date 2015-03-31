! mc_nvt_lj.f90
PROGRAM mc_nvt_lj

! Placeholder for f90 code

        COMMON / block1 / rx, ry, rz

c    *******************************************************************
c    ** monte carlo simulation PROGRAM in the constant-nvt ensemble.  **
c    **                                                               **
c    ** this PROGRAM takes a configuration of lennard jones atoms     **
c    ** and performs a conventional nvt mc simulation. the box is of  **
c    ** unit length, -0.5 to +0.5 and there are no lookup tables.     **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                   number of molecules               **
c    ** INTEGER nstep               maximum number of cycles          **
c    ** REAL    rx(n),ry(n),rz(n)   positions                         **
c    ** REAL    dens                reduced density                   **
c    ** REAL    temp                reduced temperature               **
c    ** REAL    sigma               reduced lj diameter               **
c    ** REAL    rmin                minimum reduced pair separation   **
c    ** REAL    rcut                reduced cutoff distance           **
c    ** REAL    drmax               reduced maximum displacement      **
c    ** REAL    v                   the potential energy              **
c    ** REAL    w                   the virial                        **
c    ** REAL    pres                the pressure                      **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** the PROGRAM takes in a configuration of atoms                 **
c    ** and runs a monte carlo simulation at the given temperature    **
c    ** for the specified number of cycles.                           **
c    **                                                               **
c    ** units:                                                        **
c    **                                                               **
c    ** the PROGRAM uses lennard-jones units for user input and       **
c    ** output but conducts the simulation in a box of unit length.   **
c    ** for example, for a boxlength l, and lennard-jones parameters  **
c    ** epsilon and sigma, the units are:                             **
c    **                                                               **
c    **     property       lj  units            PROGRAM units         **
c    **                                                               **
c    **     temp           epsilon/k            epsilon/k             **
c    **     pres           epsilon/sigma**3     epsilon/l**3          **
c    **     v              epsilon              epsilon               **
c    **     dens           1/sigma**3           1/l**3                **
c    **                                                               **
c    ** routines referenced:                                          **
c    **                                                               **
c    ** SUBROUTINE sumup ( rcut, rmin, sigma, ovrlap, v, w )          **
c    **    calculates the total potential energy for a configuration  **
c    ** SUBROUTINE energy ( rxi, ryi, rzi, i, rcut, sigma, v, w )     **
c    **    calculates the potential energy of atom i WITH all the     **
c    **    other atoms in the liquid                                  **
c    ** SUBROUTINE readcn (cnfile )                                   **
c    **    reads in a configuration                                   **
c    ** SUBROUTINE writcn ( cnfile )                                  **
c    **    writes out a configuration                                 **
c    ** REAL FUNCTION ranf ( dummy )                                  **
c    **    returns a uniform random number between zero and one       **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )

        REAL        rx(n), ry(n), rz(n)

        REAL        drmax, dens, temp, denslj, sigma, rmin, rcut, beta
        REAL        ranf, dummy, acm, acatma, pi, ratio, sr9, sr3
        REAL        v, vnew, vold, vend, vn, deltv, deltvb, vs
        REAL        w, wend, wnew, wold, pres, deltw, ws, ps
        REAL        vlrc, vlrc6, vlrc12, wlrc, wlrc6, wlrc12
        REAL        rxiold, ryiold, rziold, rxinew, ryinew, rzinew
        REAL        avv, avp, avw, acv, acp, acvsq, acpsq, flv, flp
        INTEGER     step, i, nstep, iprint, isave, iratio
        LOGICAL     ovrlap
        CHARACTER   title*80, cnfile*30

        PARAMETER ( pi = 3.1415927 )

c       ****************************************************************

c    ** READ input DATA **

        WRITE(*,'(1h1,'' **** program mclj ****                  ''/)')
        WRITE(*,'('' constant-nvt monte carlo program            '' )')
        WRITE(*,'('' for lennard jones atoms                      '')')

        WRITE(*,'('' enter the run title                          '')')
        READ (*,'(a)') title
        WRITE(*,'('' enter number of cycles                       '')')
        READ (*,*) nstep
        WRITE(*,'('' enter number of steps between output lines   '')')
        READ (*,*) iprint
        WRITE(*,'('' enter number of steps between data saves     '')')
        READ (*,*) isave
        WRITE(*,'('' enter interval for update of max. displ.     '')')
        READ (*,*) iratio
        WRITE(*,'('' enter the configuration file name            '')')
        READ (*,'(a)') cnfile
        WRITE(*,'('' enter the following in lennard-jones units '',/)')
        WRITE(*,'('' enter the density                            '')')
        READ (*,*) dens
        WRITE(*,'('' enter the temperature                        '')')
        READ (*,*) temp
        WRITE(*,'('' enter the potential cutoff distance          '')')
        READ (*,*) rcut

c    ** WRITE input DATA **

        WRITE(*,'(       //1x                    ,a     )') title
        WRITE(*,'('' number of atoms           '',i10   )') n
        WRITE(*,'('' number of cycles          '',i10   )') nstep
        WRITE(*,'('' output frequency          '',i10   )') iprint
        WRITE(*,'('' save frequency            '',i10   )') isave
        WRITE(*,'('' ratio update frequency    '',i10   )') iratio
        WRITE(*,'('' configuration file  name  '',a     )') cnfile
        WRITE(*,'('' temperature               '',f10.4 )') temp
        WRITE(*,'('' density                   '',f10.4 )') dens
        WRITE(*,'('' potential cutoff          '',f10.4 )') rcut

c    ** READ initial configuration **

        CALL readcn ( cnfile )

c    ** convert input DATA to PROGRAM units **

        beta   = 1.0 / temp
        sigma  = ( dens / REAL ( n ) ) ** ( 1.0 / 3.0 )
        rmin   = 0.70 * sigma
        rcut   = rcut * sigma
        drmax  = 0.15 * sigma
        denslj = dens
        dens   = dens / ( sigma ** 3 )

        IF ( rcut .GT. 0.5 ) STOP ' cut-off too large '

c    ** zero accumulators **

        acv    = 0.0
        acvsq  = 0.0
        acp    = 0.0
        acpsq  = 0.0
        flv    = 0.0
        flp    = 0.0
        acm    = 0.0
        acatma = 0.0

c    ** calculate long range corrections    **
c    ** specific to the lennard jones fluid **

        sr3 = ( sigma / rcut ) ** 3
        sr9 = sr3 ** 3

        vlrc12 =   8.0 * pi * denslj * REAL ( n ) * sr9 / 9.0
        vlrc6  = - 8.0 * pi * denslj * REAL ( n ) * sr3 / 3.0
        vlrc   =   vlrc12 + vlrc6
        wlrc12 =   4.0  * vlrc12
        wlrc6  =   2.0  * vlrc6
        wlrc   =   wlrc12 + wlrc6

c    ** WRITE out some useful information **

        WRITE(*,'('' sigma/box              =  '',f10.4)')  sigma
        WRITE(*,'('' rmin/box               =  '',f10.4)')  rmin
        WRITE(*,'('' rcut/box               =  '',f10.4)')  rcut
        WRITE(*,'('' lrc for <v>            =  '',f10.4)')  vlrc
        WRITE(*,'('' lrc for <w>            =  '',f10.4)')  wlrc

c    ** calculate initial energy and check for overlaps **

        CALL sumup ( rcut, rmin, sigma, ovrlap, v, w )

        IF ( ovrlap ) STOP ' overlap in initial configuration '

        vs = ( v + vlrc ) / REAL ( n )
        ws = ( w + wlrc ) / REAL ( n )
        ps = dens * temp + w + wlrc
        ps = ps * sigma ** 3

        WRITE(*,'('' initial v              =  '', f10.4 )' ) vs
        WRITE(*,'('' initial w              =  '', f10.4 )' ) ws
        WRITE(*,'('' initial p              =  '', f10.4 )' ) ps

        WRITE(*,'(//'' start of markov chain               ''//)')
        WRITE(*,'(''  nmove     ratio       v/n            p''/)')

c    *******************************************************************
c    ** loops over all cycles and all molecules                       **
c    *******************************************************************

        DO 100 step = 1, nstep

           DO 99 i = 1, n

              rxiold = rx(i)
              ryiold = ry(i)
              rziold = rz(i)

c          ** calculate the energy of i in the old configuration **

              CALL energy ( rxiold, ryiold, rziold, i, rcut, sigma,
     :                      vold, wold                            )

c          ** move i and pickup the central image **

              rxinew = rxiold + ( 2.0 * ranf ( dummy ) - 1.0 ) * drmax
              ryinew = ryiold + ( 2.0 * ranf ( dummy ) - 1.0 ) * drmax
              rzinew = rziold + ( 2.0 * ranf ( dummy ) - 1.0 ) * drmax

              rxinew = rxinew - ANINT ( rxinew )
              ryinew = ryinew - ANINT ( ryinew )
              rzinew = rzinew - ANINT ( rzinew )

c          ** calculate the energy of i in the NEW configuration **

              CALL energy ( rxinew, ryinew, rzinew, i, rcut, sigma,
     :                      vnew, wnew                            )

c          ** check for acceptance **

              deltv  = vnew - vold
              deltw  = wnew - wold
              deltvb = beta * deltv

              IF ( deltvb .LT. 75.0 ) THEN

                 IF ( deltv .LE. 0.0 ) THEN

                    v      = v + deltv
                    w      = w + deltw
                    rx(i)  = rxinew
                    ry(i)  = ryinew
                    rz(i)  = rzinew
                    acatma = acatma + 1.0

                 ELSEIF ( EXP ( - deltvb ) .GT. ranf ( dummy ) ) THEN

                    v      = v + deltv
                    w      = w + deltw
                    rx(i)  = rxinew
                    ry(i)  = ryinew
                    rz(i)  = rzinew
                    acatma = acatma + 1.0

                 ENDIF

              ENDIF

              acm = acm + 1.0

c          ** calculate instantaneous values **

              vn     = ( v + vlrc ) / REAL ( n )
              pres   = dens * temp + w + wlrc

c          ** convert pressure to lj units **

              pres   = pres * sigma ** 3

c          ** accumulate averages **

              acv    = acv   + vn
              acp    = acp   + pres
              acvsq  = acvsq + vn * vn
              acpsq  = acpsq + pres * pres

c          *************************************************************
c          ** ends loop over atoms                                    **
c          *************************************************************

99         CONTINUE

c       ** perform periodic operations  **

           IF ( MOD ( step, iratio ) .EQ. 0 ) THEN

c          ** adjust maximum displacement **

              ratio = acatma / REAL ( n * iratio )

              IF ( ratio .GT. 0.5 ) THEN

                 drmax  = drmax  * 1.05

              ELSE

                 drmax  = drmax  * 0.95

              ENDIF

              acatma = 0.0

           ENDIF

           IF ( MOD ( step, iprint ) .EQ. 0 ) THEN

c          ** WRITE out runtime information **

              WRITE(*,'(i8,3f12.6)') INT(acm), ratio, vn, pres

           ENDIF

           IF ( MOD ( step, isave ) .EQ. 0 ) THEN

c          ** WRITE out the configuration at intervals **

              CALL writcn ( cnfile )

           ENDIF

100     CONTINUE

c    *******************************************************************
c    ** ends the loop over cycles                                     **
c    *******************************************************************

        WRITE(*,'(//'' end of markov chain          ''//)')

c    ** checks FINAL VALUE of the potential energy is consistent **

        CALL sumup ( rcut, rmin, sigma, ovrlap, vend, wend )

        IF ( ABS ( vend - v ) .GT. 1.0e-03 ) THEN

           WRITE(*,'('' problem with energy,'')')
           WRITE(*,'('' vend              = '', e20.6)') vend
           WRITE(*,'('' v                 = '', e20.6)') v

        ENDIF

c    ** WRITE out the FINAL configuration from the run **

        CALL writcn ( cnfile )

c    ** calculate and WRITE out running averages **

        avv   = acv / acm
        acvsq = ( acvsq / acm ) - avv ** 2
        avp   = acp / acm
        acpsq = ( acpsq / acm ) - avp ** 2

c    ** calculate fluctuations **

        IF ( acvsq .GT. 0.0 ) flv = SQRT ( acvsq )
        IF ( acpsq .GT. 0.0 ) flp = SQRT ( acpsq )

        WRITE(*,'(/'' averages ''/ )')
        WRITE(*,'('' <v/n>   = '',f10.6)') avv
        WRITE(*,'('' <p>     = '',f10.6)') avp

        WRITE(*,'(/'' fluctuations ''/)')

        WRITE(*,'('' fluctuation in <v/n> = '',f10.6)') flv
        WRITE(*,'('' fluctuation in <p>   = '',f10.6)') flp
        WRITE(*,'(/'' end of simulation '')')

        STOP
        END



        SUBROUTINE sumup ( rcut, rmin, sigma, ovrlap, v, w )

        COMMON / block1 / rx, ry, rz

c    *******************************************************************
c    ** calculates the total potential energy for a configuration.    **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER n                 the number of atoms                 **
c    ** REAL    rx(n(,ry(n),rz(n) the positions of the atoms          **
c    ** REAL    v                 the potential energy                **
c    ** REAL    w                 the virial                          **
c    ** LOGICAL ovrlap            true for substantial atom overlap   **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** the SUBROUTINE returns the total potential energy at the      **
c    ** beginning and END of the run.                                 **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )
        REAL        rx(n), ry(n), rz(n)
        REAL        sigma, rmin, rcut, v, w
        LOGICAL     ovrlap

        REAL        rcutsq, rminsq, sigsq, rxij, ryij, rzij
        REAL        rxi, ryi, rzi, vij, wij, sr2, sr6, rijsq
        INTEGER     i, j

c    *******************************************************************

        ovrlap = .FALSE.
        rcutsq = rcut * rcut
        rminsq = rmin * rmin
        sigsq  = sigma * sigma

        v      = 0.0
        w      = 0.0

c    ** loop over all the pairs in the liquid **

        DO 100 i = 1, n - 1

           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)

           DO 99 j = i + 1, n

              rxij  = rxi - rx(j)
              ryij  = ryi - ry(j)
              rzij  = rzi - rz(j)

c          ** minimum image the pair separations **

              rxij  = rxij - ANINT ( rxij )
              ryij  = ryij - ANINT ( ryij )
              rzij  = rzij - ANINT ( rzij )
              rijsq = rxij * rxij + ryij * ryij + rzij * rzij

              IF ( rijsq .LT. rminsq ) THEN

                 ovrlap = .TRUE.
                 RETURN

              ELSEIF ( rijsq .LT. rcutsq ) THEN

                 sr2 = sigsq / rijsq
                 sr6 = sr2 * sr2 * sr2
                 vij = sr6 * ( sr6 - 1.0 )
                 wij = sr6 * ( sr6 - 0.5 )
                 v   = v + vij
                 w   = w + wij

              ENDIF

99         CONTINUE

100     CONTINUE

        v = 4.0 * v
        w = 48.0 * w / 3.0

        RETURN
        END



        SUBROUTINE energy ( rxi, ryi, rzi, i, rcut, sigma, v, w )

        COMMON / block1 / rx, ry, rz

c    *******************************************************************
c    ** returns the potential energy of atom i WITH all other atoms.  **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER i                 the atom of interest                **
c    ** INTEGER n                 the number of atoms                 **
c    ** REAL    rx(n),ry(n),rz(n) the atom positions                  **
c    ** REAL    rxi,ryi,rzi       the coordinates of atom i           **
c    ** REAL    v                 the potential energy of atom i      **
c    ** REAL    w                 the virial of atom i                **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** this SUBROUTINE is used to calculate the change of energy     **
c    ** during a trial move of atom i. it is called before and        **
c    ** after the random displacement of i.                           **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )
        REAL        rx(n), ry(n), rz(n)
        REAL        rcut, sigma, rxi, ryi, rzi, v, w
        INTEGER     i

        REAL        rcutsq, sigsq, sr2, sr6
        REAL        rxij, ryij, rzij, rijsq, vij, wij
        INTEGER     j

c     ******************************************************************

        rcutsq = rcut * rcut
        sigsq  = sigma * sigma

        v      = 0.0
        w      = 0.0

c    ** loop over all molecules except i  **

        DO 100 j = 1, n

           IF ( i .NE. j ) THEN

              rxij  = rxi - rx(j)
              ryij  = ryi - ry(j)
              rzij  = rzi - rz(j)

              rxij  = rxij - ANINT ( rxij )
              ryij  = ryij - ANINT ( ryij )
              rzij  = rzij - ANINT ( rzij )

              rijsq = rxij * rxij + ryij * ryij + rzij * rzij

              IF ( rijsq .LT. rcutsq ) THEN

                 sr2 = sigsq / rijsq
                 sr6 = sr2 * sr2 * sr2
                 vij = sr6 * ( sr6 - 1.0 )
                 wij = sr6 * ( sr6 - 0.5 )
                 v   = v + vij
                 w   = w + wij

              ENDIF

           ENDIF

100     CONTINUE

        v = 4.0 * v
        w = 48.0 * w / 3.0

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



        SUBROUTINE readcn ( cnfile )

        COMMON / block1 / rx, ry, rz

c    *******************************************************************
c    ** SUBROUTINE to READ in the configuration from unit 10          **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )
        CHARACTER   cnfile*(*)
        REAL        rx(n), ry(n), rz(n)

        INTEGER     cnunit
        PARAMETER ( cnunit = 10 )

        INTEGER     nn

c   ********************************************************************

        OPEN ( unit = cnunit, file = cnfile, status = 'old',
     :         form = 'unformatted'                        )

        READ ( cnunit ) nn
        IF ( nn .NE. n ) STOP 'n error in readcn'
        READ ( cnunit ) rx, ry, rz

        CLOSE ( unit = cnunit )

        RETURN
        END



        SUBROUTINE writcn ( cnfile )

        COMMON / block1 / rx, ry, rz

c    *******************************************************************
c    ** SUBROUTINE to WRITE out the configuration to unit 10          **
c    *******************************************************************

        INTEGER     n
        PARAMETER ( n = 108 )
        CHARACTER   cnfile*(*)
        REAL        rx(n), ry(n), rz(n)

        INTEGER     cnunit
        PARAMETER ( cnunit = 10 )

c   ********************************************************************

        OPEN ( unit = cnunit, file = cnfile, status = 'unknown',
     :         form = 'unformatted'                            )

        WRITE ( cnunit ) n
        WRITE ( cnunit ) rx, ry, rz

        CLOSE ( unit = cnunit )

        RETURN
        END



