! mc_npt_lj
program mcnpt
! Placeholder for f90 code
********************************************************************************
** fiche f.12.  constant-npt monte carlo algorithm.                           **
** this fortran code is intended to illustrate points made in the text.       **
** to our knowledge it works correctly.  however it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************


        common / block1 / rx, ry, rz

c    *******************************************************************
c    ** monte carlo simulation in the constant-npt ensemble.          **
c    **                                                               **
c    ** this program takes a configuration of lennard jones atoms     **
c    ** and performs a monte carlo simulation at constant npt. the    **
c    ** box is in units of sigma the lennard jones diameter.          **
c    ** there are no lookup tables included.                          **
c    **                                                               **
c    ** reference:                                                    **
c    **                                                               **
c    ** mcdonald, chem. phys. lett. 3, 241, 1969.                     **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** integer n                             number of molecules     **
c    ** real    rx(n),ry(n),rz(n)             positions               **
c    ** real    vol                           volume                  **
c    ** real    box                           box length              **
c    ** real    dens                          reduced density         **
c    ** real    temp                          reduced temperature     **
c    ** real    sigma                         reduced lj diameter     **
c    ** real    drmax                         maximum displacement    **
c    ** real    v                             the potential energy    **
c    ** real    w                             the virial              **
c    ** real    presur                        required pressure       **
c    ** real    dboxmx                        max change in box       **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** conducts monte carlo simulation at constant pressure for a    **
c    ** specified number of cycles from a given initial configuration.**
c    **                                                               **
c    ** units:                                                        **
c    **                                                               **
c    ** this program uses the usual reduced lj units. in particular   **
c    ** the box length is in units of sigma.                          **
c    **                                                               **
c    ** routines referenced:                                          **
c    **                                                               **
c    ** subroutine energy ( rxi, ryi, rzi, i, rcut, box, v12, v6,     **
c    **    :                w12, w6 )                                 **
c    **    calculates the energy and virial for atom i in the fluid   **
c    ** real function ranf ( dummy )                                  **
c    **    returns a uniform random number between zero and one       **
c    ** subroutine readcn ( cnfile, box )                             **
c    **    reads in configuration and box variables                   **
c    ** subroutine sumup ( rcut, rmin, ovrlap, box, v12, v6, w12, w6 )**
c    **    calculates potential and virial for a configuration        **
c    ** subroutine writcn ( cnfile, box )                             **
c    **    writes out configuration and box variables                 **
c    *******************************************************************

        integer     n
        parameter ( n = 108 )

        real        rx(n), ry(n), rz(n)

        integer     step, nstep, iprint, iratio, iratb, i
        real        acv, acp, acd, acm, acatma, acboxa, norm
        real        acvsq, acpsq, acdsq
        real        avv, avp, avd
        real        flv, flp, fld
        real        dens, temp, rcut, rmin, presur, box, vol, pres, vn
        real        boxinv, boxnew, ratbox, rat12, rat6, dvol, dpv
        real        delthb, drmax, dboxmx, beta, ranf, dummy, ratio
        real        rrbox, bratio, deltvb, rcutn
        real        rxiold, ryiold, rziold, rxinew, ryinew, rzinew
        real        v12old, v6old, v12new, v6new, vs
        real        w12old, w6old, w12new, w6new, ws, ps
        real        delv12, delv6, delw12, delw6, deltv
        real        v6, w6, v12, w12, vlrc, vlrcn, wlrc, wlrcn
        real        sr3, sr9, vlrc6, wlrc6, vlrc12, wlrc12, pi
        character   title*80, cnfile*30
        logical     ovrlap

        parameter ( pi = 3.1415927 )

c    *******************************************************************

        write(*,'(1h1,'' **** program mcnpt ****               '')')
        write(*,'(//'' monte carlo in a constant-npt ensemble  '')')
        write(*,'(  '' lennard-jones atoms                     '')')

c    ** basic simulation parameters **

        write(*,'('' enter run title                           '')')
        read (*,'(a)') title
        write(*,'('' enter configuration filename              '')')
        read (*,'(a)') cnfile
        write(*,'('' enter number of cycles                    '')')
        read (*,*) nstep
        write(*,'('' enter interval between prints in cycles   '')')
        read (*,*) iprint
        write(*,'('' enter interval for update of maximum      '')')
        write(*,'('' displacement of atoms in cycles           '')')
        read (*,*) iratio
        write(*,'('' enter interval for update of maximum      '')')
        write(*,'('' displacement of the box in cycles         '')')
        read (*,*) iratb
        write(*,'(/'' enter the following in l-j reduced units  ''/)')
        write(*,'('' enter potential cutoff                    '')')
        read (*,*) rcut
        write(*,'('' enter desired pressure                    '')')
        read (*,*) presur
        write(*,'('' enter desired temperature                 '')')
        read (*,*) temp

        write(*,'(//1x,a)') title
        write(*,'('' configuration filename            '',a)') cnfile
        write(*,'('' number of cycles                = '',i6)') nstep
        write(*,'('' print interval                  = '',i6)') iprint
        write(*,'('' ratio update interval for atoms = '',i6)') iratio
        write(*,'('' ratio update interval for box   = '',i6)') iratb
        write(*,'('' potential cutoff                = '',f10.5)')rcut
        write(*,'('' desired pres.                   = '',f10.5)')presur
        write(*,'('' desired temp                    = '',f10.5)')temp

c    ** readcn reads in initial configuration **
c    ** origin should be at centre of box     **

        call readcn ( cnfile, box )

c    ** set dependent parameters **

        vol    = box ** 3
        boxinv = 1.0 / box
        dens   = real ( n ) / vol

        if ( rcut .gt. ( 0.5 * box ) ) stop 'cut-off too large'

        dboxmx = box / 40.0
        drmax  = 0.15
        rmin   = 0.70
        beta   = 1.0 / temp

        write(*,'('' initial density                 = '',f10.5)') dens

c    ** calculate long-range corrections for lj potential.    **
c    ** 6 is for attractive contributions 12 is for repulsive **

        sr3    =    ( 1.0 / rcut ) ** 3
        sr9    =    sr3 ** 3
        vlrc12 =    8.0 * pi * dens * real ( n ) * sr9 / 9.0
        vlrc6  = -  8.0 * pi * dens * real ( n ) * sr3 / 3.0
        wlrc12 =    4.0 * vlrc12
        wlrc6  =    2.0 * vlrc6
        vlrc   =    vlrc12 + vlrc6
        wlrc   =    wlrc12 + wlrc6

c    ** zero accumulators **

        acm    = 0.0
        acatma = 0.0
        acboxa = 0.0

        acv = 0.0
        acp = 0.0
        acd = 0.0

        acvsq = 0.0
        acpsq = 0.0
        acdsq = 0.0

        flv = 0.0
        flp = 0.0
        fld = 0.0

c    ** calculate initial energy and virial **

        call sumup ( rcut, rmin, ovrlap, box, v12, v6, w12, w6 )

        if ( ovrlap ) stop 'overlap in initial configuration'

c    ** calculate the initial energy and virial **

        vs = ( v12 + v6 + vlrc ) / real ( n )
        ws = ( w12 + w6 + wlrc ) / real ( n )
        ps = dens * temp + ( w12 + w6 + wlrc ) / vol

c    ** add long range corrections **
c    ** into the energy and virial **

        v12 = v12 + vlrc12
        v6  = v6  + vlrc6
        w12 = w12 + wlrc12
        w6  = w6  + wlrc6

        write(*,'('' initial v/n                     = '',f10.6)') vs
        write(*,'('' initial w/n                     = '',f10.6)') ws
        write(*,'('' initial p                       = '',f10.6)') ps
        write(*,'(//'' **** start of markov chain ****'')')
        write(*,10001)

c    *******************************************************************
c    ** main loop starts                                              **
c    *******************************************************************

        do 100 step = 1, nstep

c       ** loop over atoms starts **

           do 97 i = 1, n

              rxiold = rx(i)
              ryiold = ry(i)
              rziold = rz(i)

c          ** calculate v for an atom in old state **

              call energy ( rxiold, ryiold, rziold, i, rcut, box,
     :                      v12old, v6old, w12old, w6old )

c          ** move atom i and pickup central image **

              rxinew = rxiold + ( 2.0 * ranf ( dummy ) - 1.0 ) * drmax
              ryinew = ryiold + ( 2.0 * ranf ( dummy ) - 1.0 ) * drmax
              rzinew = rziold + ( 2.0 * ranf ( dummy ) - 1.0 ) * drmax

              rxinew = rxinew - anint ( rxinew * boxinv ) * box
              ryinew = ryinew - anint ( ryinew * boxinv ) * box
              rzinew = rzinew - anint ( rzinew * boxinv ) * box

c          ** calculate v for atom in new state **

              call energy ( rxinew, ryinew, rzinew, i, rcut, box,
     :                      v12new, v6new, w12new, w6new )

c          ** check for acceptance **

              delv12 = v12new - v12old
              delv6  = v6new  - v6old
              delw12 = w12new - w12old
              delw6  = w6new  - w6old
              deltv  = delv12 + delv6
              deltvb = beta * deltv

              if ( deltvb .lt. 75.0 ) then

                 if ( deltv .le. 0.0 ) then

                    v12    = v12 + delv12
                    v6     = v6  + delv6
                    w12    = w12 + delw12
                    w6     = w6  + delw6
                    rx(i)  = rxinew
                    ry(i)  = ryinew
                    rz(i)  = rzinew
                    acatma = acatma + 1.0

                 elseif ( exp ( - deltvb ) .gt. ranf ( dummy ) ) then

                    v12    = v12 + delv12
                    v6     = v6  + delv6
                    w12    = w12 + delw12
                    w6     = w6  + delw6
                    rx(i)  = rxinew
                    ry(i)  = ryinew
                    rz(i)  = rzinew
                    acatma = acatma + 1.0

                 endif

              endif

              vn   = ( v12 + v6 ) / real ( n )
              pres = dens * temp +  ( w12 + w6 )  / vol

c       ** increment accumulators **

              acm = acm + 1.0
              acv = acv + vn
              acp = acp + pres
              acd = acd + dens

              acvsq = acvsq + vn ** 2
              acpsq = acpsq + pres ** 2
              acdsq = acdsq + dens ** 2

97         continue

c       ** ends loop over atoms in one cycle  **

c       ** attempt a box move **

           boxnew = box + ( 2.0 * ranf ( dummy ) - 1.0 ) * dboxmx
           ratbox = box / boxnew
           rrbox  = 1.0 / ratbox
           rcutn  = rcut * rrbox

c       ** calculate scaling parameters **

           rat6   = ratbox ** 6
           rat12  = rat6 * rat6

c       ** scale energy, and virial including lrc **

           v12new = v12  * rat12
           v6new  = v6   * rat6
           w12new = w12  * rat12
           w6new  = w6   * rat6

c       ** calculate change in energy and volume **

           deltv  = v12new + v6new - v12 - v6
           dpv    = presur * ( boxnew ** 3 - vol )
           dvol   = 3.0 * temp * real ( n ) * alog ( ratbox )
           delthb = beta * ( deltv + dpv + dvol )

c       ** check for acceptance **

           if ( delthb .lt. 75.0 ) then

              if ( delthb .le. 0.0 ) then

                 v12    = v12new
                 v6     = v6new
                 w12    = w12new
                 w6     = w6new

                 do 98 i = 1, n

                    rx(i) = rx(i) * rrbox
                    ry(i) = ry(i) * rrbox
                    rz(i) = rz(i) * rrbox

98               continue

                 box    = boxnew
                 rcut   = rcutn
                 acboxa = acboxa + 1.0

              elseif ( exp ( - delthb ) .gt. ranf ( dummy ) ) then

                 v12    = v12new
                 v6     = v6new
                 w12    = w12new
                 w6     = w6new

                 do 99 i = 1, n

                    rx(i) = rx(i) * rrbox
                    ry(i) = ry(i) * rrbox
                    rz(i) = rz(i) * rrbox

99               continue

                 box    = boxnew
                 rcutn  = rcut
                 acboxa = acboxa + 1.0

              endif

           endif

           boxinv = 1.0 / box
           vol    = box ** 3
           dens   = real ( n ) / vol

c       ** calculate energy and pressure **

           vn   = ( v12 + v6 ) / real ( n )
           pres = dens * temp + ( w12 + w6 ) / vol

c       ** increment accumulators **

           acm = acm + 1.0
           acv = acv + vn
           acp = acp + pres
           acd = acd + dens

           acvsq = acvsq + vn ** 2
           acpsq = acpsq + pres ** 2
           acdsq = acdsq + dens ** 2

c       ** ends attempted box move **

c       ** perform periodic operations **

           if ( mod ( step, iratio ) .eq. 0 ) then

c          ** adjust maximum displacement for atoms **

              ratio = acatma / real ( n * iratio )

              if ( ratio .gt. 0.5 ) then

                 drmax = drmax * 1.05

              else

                 drmax = drmax * 0.95

              endif

              acatma = 0.0

           endif

           if ( mod ( step, iratb ) .eq. 0 ) then

c          ** adjust maximum displacement for the box **

              bratio = acboxa / real ( iratb )

              if ( bratio .gt. 0.5 ) then

                 dboxmx = dboxmx * 1.05

              else

                 dboxmx = dboxmx * 0.95

              endif

              acboxa = 0.0

           endif

           if ( mod ( step, iprint ) .eq. 0 ) then

c          ** optionally print information **

              write(*,'(1x,i8,5(2x,f10.5))')
     :                 step, vn, pres, dens, ratio, bratio

           endif

100     continue

c    *******************************************************************
c    ** main loop ends                                                **
c    *******************************************************************

        write(*,'(/1x,''**** end of markov chain **** ''//)')

c    ** write out final configuration and boxlength **

        call writcn ( cnfile, box )

c    ** write out final averages **

        norm = real ( acm )
        avv  = acv / norm
        avp  = acp / norm
        avd  = acd / norm

        acvsq = ( acvsq / norm ) - avv ** 2
        acpsq = ( acpsq / norm ) - avp ** 2
        acdsq = ( acdsq / norm ) - avd ** 2

        if ( acvsq .gt. 0.0 ) flv = sqrt ( acvsq )
        if ( acpsq .gt. 0.0 ) flp = sqrt ( acpsq )
        if ( acdsq .gt. 0.0 ) fld = sqrt ( acdsq )

        write(*, 10002)
        write(*,'('' averages'',3(2x,f10.5))') avv, avp, avd
        write(*,'('' flucts  '',3(2x,f10.5))') flv, flp, fld

        stop

10001   format(//1x,'  cycle   ..potent..  ..pressure.. ..density..',
     :              ' ..ratio..  ..bratio..')
10002   format(//1x,'  cycle   ..potent..  ..pressure.. ..density..')

        end



        subroutine sumup ( rcut, rmin, ovrlap, box, v12, v6, w12, w6 )

        common / block1 / rx, ry, rz

c    *******************************************************************
c    ** calculates the total potential energy for a configuration.    **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** integer n                 the number of atoms                 **
c    ** real    rx(n(,ry(n),rz(n) the positions of the atoms          **
c    ** real    v                 the potential energy                **
c    ** real    w                 the virial                          **
c    ** logical ovrlap            true for substantial atom overlap   **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** the subroutine returns the total potential energy at the      **
c    ** beginning and end of the run.                                 **
c    *******************************************************************

        integer     n
        parameter ( n = 108 )
        real        rx(n), ry(n), rz(n)
        real        rmin, rcut, v12, v6, w12, w6, box
        logical     ovrlap

        real        rcutsq, rminsq, vij12, vij6
        real        rxi, ryi, rzi, rxij, ryij, rzij
        real        sr2, sr6, rijsq, boxinv
        integer     i, j

c    *******************************************************************

        ovrlap = .false.
        rcutsq = rcut * rcut
        rminsq = rmin * rmin
        boxinv = 1.0 / box

        v12    = 0.0
        v6     = 0.0
        w12    = 0.0
        w6     = 0.0

c    ** loop over all the pairs in the liquid **

        do 100 i = 1, n - 1

           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)

           do 99 j = i + 1, n

              rxij  = rxi - rx(j)
              ryij  = ryi - ry(j)
              rzij  = rzi - rz(j)

              rxij  = rxij - anint ( rxij * boxinv ) * box
              ryij  = ryij - anint ( ryij * boxinv ) * box
              rzij  = rzij - anint ( rzij * boxinv ) * box

              rijsq = rxij * rxij + ryij * ryij + rzij * rzij

              if ( rijsq .lt. rminsq ) then

                 ovrlap = .true.
                 return

              elseif ( rijsq .lt. rcutsq ) then

                 sr2   = 1.0 / rijsq
                 sr6   = sr2 * sr2 * sr2
                 vij12 = sr6 *  sr6
                 vij6  = - sr6
                 v12   = v12 + vij12
                 v6    = v6  + vij6
                 w12   = w12 + vij12
                 w6    = w6  + vij6 * 0.5

              endif

99         continue

100     continue

        v12 = 4.0 * v12
        v6  = 4.0 * v6
        w12 = 48.0 * w12 / 3.0
        w6  = 48.0 * w6  / 3.0

        return
        end



        subroutine energy ( rxi, ryi, rzi, i, rcut, box,
     :                      v12, v6, w12, w6 )

        common / block1 / rx, ry, rz

c    *******************************************************************
c    ** calculates the potential energy of i with all other atoms     **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** integer i                 the atom of interest                **
c    ** integer n                 the number of atoms                 **
c    ** real    rx(n),ry(n),rz(n) the atom positions                  **
c    ** real    rxi,ryi,rzi       the coordinates of atom i           **
c    ** real    v                 the potential energy of atom i      **
c    ** real    w                 the virial of atom i                **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** this subroutine is used to calculate the change of energy     **
c    ** during a trial move of atom i. it is called before and        **
c    ** after the random displacement of i.                           **
c    *******************************************************************


        integer     n
        parameter ( n = 108 )
        real        rx(n), ry(n), rz(n)
        real        rcut, box, rxi, ryi, rzi, v12, v6, w12, w6
        integer     i

        real        rcutsq, vij12, vij6
        real        rxij, ryij, rzij, rijsq, sr2, sr6, boxinv
        integer     j

c     ******************************************************************

        rcutsq = rcut * rcut
        boxinv = 1.0 / box

        v12    = 0.0
        v6     = 0.0
        w12    = 0.0
        w6     = 0.0

c    ** loop over all molecules except i  **

        do 100 j = 1, n

        if ( i .ne. j ) then

           rxij  = rxi - rx(j)
           ryij  = ryi - ry(j)
           rzij  = rzi - rz(j)

           rxij  = rxij - anint ( rxij * boxinv ) * box
           ryij  = ryij - anint ( ryij * boxinv ) * box
           rzij  = rzij - anint ( rzij * boxinv ) * box

           rijsq = rxij * rxij + ryij * ryij + rzij * rzij

           if ( rijsq .lt. rcutsq ) then

              sr2   = 1.0 / rijsq
              sr6   = sr2 * sr2 * sr2
              vij12 = sr6 * sr6
              vij6  = - sr6
              v12   = v12 + vij12
              v6    = v6  + vij6
              w12   = w12 + vij12
              w6    = w6  + vij6 * 0.5

           endif

        endif

100     continue

        v12 = 4.0 * v12
        v6  = 4.0 * v6
        w12 = 48.0 * w12 / 3.0
        w6  = 48.0 * w6  / 3.0

        return
        end



        subroutine readcn ( cnfile, box )

        common / block1 / rx, ry, rz

c    *******************************************************************
c    ** subroutine to read in the configuration from unit 10          **
c    *******************************************************************

        integer     n
        parameter ( n = 108 )
        character   cnfile*(*)
        real        rx(n), ry(n), rz(n)
        real        box

        integer     cnunit
        parameter ( cnunit = 10 )
        integer     nn

c   ********************************************************************

        open ( unit = cnunit, file = cnfile, status = 'old',
     :         form = 'unformatted'                        )

        read ( cnunit ) nn, box
        if ( nn .ne. n ) stop 'n error in readcn'
        read ( cnunit ) rx, ry, rz

        close ( unit = cnunit )

        return
        end



        subroutine writcn ( cnfile, box )

        common / block1 / rx, ry, rz

c    *******************************************************************
c    ** subroutine to write out the configuration to unit 10          **
c    *******************************************************************

        integer     n
        parameter ( n = 108 )
        character   cnfile*(*)
        real        rx(n), ry(n), rz(n)
        real        box

        integer     cnunit
        parameter ( cnunit = 10 )

c   ********************************************************************

        open ( unit = cnunit, file = cnfile, status = 'unknown',
     :         form = 'unformatted'                        )

        write ( cnunit ) n, box
        write ( cnunit ) rx, ry, rz

        close ( unit = cnunit )

        return
        end



        real function ranf ( dummy )

c    *******************************************************************
c    ** returns a uniform random variate in the range 0 to 1.         **
c    **                                                               **
c    **                 ***************                               **
c    **                 **  warning  **                               **
c    **                 ***************                               **
c    **                                                               **
c    ** good random number generators are machine specific.           **
c    ** please use the one recommended for your machine.              **
c    *******************************************************************

        integer     l, c, m
        parameter ( l = 1029, c = 221591, m = 1048576 )

        integer     seed
        real        dummy
        save        seed
        data        seed / 0 /

c    *******************************************************************

        seed = mod ( seed * l + c, m )
        ranf = real ( seed ) / m

        return
        end



