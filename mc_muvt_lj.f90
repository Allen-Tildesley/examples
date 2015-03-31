! mc_muvt_lj.f90
program mc_muvt_lj

! placeholder for f90 code

********************************************************************************
** fiche f.13.  the heart of a constant mu vt monte carlo program             **
** this fortran code is intended to illustrate points made in the text.       **
** to our knowledge it works correctly.  however it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

c    *******************************************************************
c    ** attempted creations and destructions in grand canonical mc.   **
c    **                                                               **
c    ** these routines allow for a trial destruction or creation in a **
c    ** grand canonical monte carlo program.                          **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** integer n                   number of atoms before the trial  **
c    ** integer ntrial              number of atoms during the trial  **
c    ** integer nmax                maximum number of atoms allowed   **
c    ** integer locate(nmax)        array of active atom indices      **
c    ** real    rxnew,rynew,rznew   position for addition of atom     **
c    ** real    rx(nmax) etc.       positions of current atoms        **
c    ** real    v                   potential energy + lrc            **
c    ** real    w                   virial + lrc                      **
c    ** real    deltv               change in energy                  **
c    ** real    deltw               change in virial                  **
c    ** real    temp                reduced temperature               **
c    ** real    z                   absolute activity coefficient     **
c    ** real    sigma               lennard jones diameter            **
c    ** real    rcut                reduced cutoff distance           **
c    ** real    rmin                reduced minimum separation        **
c    ** logical ovrlap              true for substantial atom overlap **
c    ** logical create              true for an accepted creation     **
c    ** logical ghost               true for an accepted destruction  **
c    **                                                               **
c    ** routines supplied:                                            **
c    **                                                               **
c    ** subroutine in ( temp, z, sigma, rcut, n, v, w, create )       **
c    **    performs a trial creation                                  **
c    ** subroutine out ( temp, z, sigma, rcut, n, v, w, ghost )       **
c    **    performs a trial destruction                               **
c    ** subroutine potin ( rxnew, rynew, rznew, n, sigma, rcut, rmin, **
c    ** :                  deltv, deltw, ovrlap )                     **
c    **    calculates the potential energy change on creation         **
c    ** subroutine potout ( ipull, n, sigma, rcut, deltv, deltw )     **
c    **    calculates the potential energy change on destruction      **
c    **                                                               **
c    ** routines referenced:                                          **
c    **                                                               **
c    ** real function ranf ( dummy )  (given in f.11)                 **
c    **    returns a uniform random variate on zero to one            **
c    ** subroutine add ( rxnew, rynew, rznew, n )                     **
c    **    updates locate after addition (given in f.14)              **
c    ** subroutine remove ( nloc, ipull, n )                          **
c    **    updates locate after removal (given in f.14)               **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** routines in and out should be called with equal probability   **
c    ** in a grand canonical monte carlo simulation. if a trial       **
c    ** creation is accepted then create is set to true. if a trial   **
c    ** destruction is accepted then ghost is set to true. the        **
c    ** routines are written for lennard-jones atoms. the box is of   **
c    ** unit length, all distances are scaled to the box length.      **
c    ** trial inputs which result in a separation of less than        **
c    ** 0.5*sigma are rejected. the long-range corrections are        **
c    ** included in v and w. all accumulators are updated in the main **
c    ** part of the program which is not given here.                  **
c    *******************************************************************

        subroutine in ( temp, z, sigma, rcut, n, v, w, create )

        common / block1 / rx, ry, rz
        common / block2 / locate

c    *******************************************************************
c    ** routine to attempt a trial creation                           **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** real    temp            temperature                           **
c    ** real    z               absolute activity                     **
c    ** real    sigma           lennard-jones diameter                **
c    ** real    rcut            cut-off distance                      **
c    ** real    v               potential energy                      **
c    ** real    w               virial                                **
c    ** integer n               number of atoms before trial creation **
c    ** logical create          true for a successful creation        **
c    *******************************************************************

        integer     nmax
        parameter ( nmax = 500 )

        real        rx(nmax), ry(nmax), rz(nmax)
        real        temp, z, sigma, rcut, v, w
        integer     locate(nmax)
        integer     n
        logical     create

        real        beta, rxnew, rynew, rznew, deltv, deltw, deltcb
        real        ranf, dummy, rmin
        integer     ntrial
        logical     ovrlap

c    *******************************************************************

        create = .false.
        beta   = 1.0 / temp
        rmin   = 0.5 * sigma
        ntrial = n + 1

        if ( ntrial .ge. nmax ) stop 'maximum number of atoms in box'

c    ** generate the position of the trial atom **

        rxnew  = ranf ( dummy ) - 0.5
        rynew  = ranf ( dummy ) - 0.5
        rznew  = ranf ( dummy ) - 0.5

c    ** calculate energy change on addition **

        call potin ( rxnew, rynew, rznew, n, sigma, rcut, rmin, deltv,
     :               deltw, ovrlap )

c    ** check for acceptance **

        if ( .not. ovrlap ) then

           deltcb = beta * deltv - log ( z / real ( ntrial ) )

           if ( deltcb .lt. 75.0 ) then

               if ( deltcb .le. 0.0 ) then

                 create = .true.

                 call add ( rxnew, rynew, rznew, n )

                 v    = v + deltv
                 w    = w + deltw
                 n    = ntrial

              else if ( exp ( - deltcb ) .gt. ranf ( dummy ) ) then

                 create = .true.

                 call add ( rxnew, rynew, rznew, n )

                 v    = v + deltv
                 w    = w + deltw
                 n    = ntrial

              endif

           endif

        endif

        return
        end



        subroutine out ( temp, z, sigma, rcut, n, v, w, ghost )

        common / block1 / rx, ry, rz
        common / block2 / locate

c    *******************************************************************
c    ** routine to attempt a trial destruction                        **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** real    temp         temperature                              **
c    ** real    z            absolute activity                        **
c    ** real    sigma        lennard-jones diameter                   **
c    ** real    rcut         cut-off distance                         **
c    ** real    v            potential energy                         **
c    ** real    w            virial                                   **
c    ** integer n            number of atoms before trial destruction **
c    ** logical ghost        true for a successful destruction        **
c    *******************************************************************

        integer     nmax
        parameter ( nmax = 500 )

        real        rx(nmax), ry(nmax), rz(nmax)
        real        temp, z, sigma, rcut, v, w
        integer     locate(nmax)
        integer     n
        logical     ghost

        real        beta, deltv, deltw, deltdb, ranf, dummy
        integer     ntrial, nloc, ipull
        logical     ovrlap

c    *******************************************************************

        ghost  = .false.
        beta   = 1.0 / temp
        ntrial = n - 1

        if ( ntrial .eq. 1 ) stop 'only one atom remains'

c    ** pick a random element from the active part of locate **

        nloc  = int ( real ( ntrial ) * ranf ( dummy ) ) + 1
        ipull = locate(nloc)

c    ** calculate energy change on removal of atom ipull **

        call potout ( ipull, n, sigma, rcut, deltv, deltw )

c    ** check for acceptance **

        deltdb = beta * deltv - log ( real ( n ) / z )

        if ( deltdb .lt. 75.0 ) then

           if ( deltdb .lt. 0.0 ) then

              ghost = .true.

              call remove ( nloc, ipull, n )

              v = v + deltv
              w = w + deltw
              n = ntrial

           else if ( exp( -deltdb ) .gt. ranf ( dummy ) ) then

              ghost = .true.

              call remove ( nloc, ipull, n )

              v = v + deltv
              w = w + deltw
              n = ntrial

           endif

        endif

        return
        end



        subroutine potin ( rxi, ryi, rzi, n, sigma, rcut, rmin, deltv,
     :                     deltw, ovrlap )

        common / block1 / rx, ry, rz
        common / block2 / locate

c    *******************************************************************
c    ** returns the potential energy change on adding an atom.        **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** integer n                 the number of atoms before addition **
c    ** real    rxi,ryi,rzi       the coordinates of the added atom   **
c    ** real    deltv             the chance in potential             **
c    ** real    deltw             the change in virial                **
c    ** real    sigma             lj diameter                         **
c    ** real    rcut              cutoff distance for potential       **
c    ** real    rmin              minimum allowed approach of atoms   **
c    ** logical ovrlap            true for substantial atom overlap   **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** this subroutine is used to calculate the change of energy     **
c    ** during a trial addition of an atom to the fluid. the long     **
c    ** range corrections is included.                                **
c    *******************************************************************

        integer     nmax
        parameter ( nmax = 500 )
        real        rx(nmax), ry(nmax), rz(nmax)
        real        rcut, rmin, sigma, rxi, ryi, rzi, deltv, deltw
        integer     n, ntrial, locate(nmax)
        logical     ovrlap

        real        rcutsq, rminsq, sigsq, sr2, sr6, sr3, sr9
        real        rxij, ryij, rzij, rijsq, vij, wij, sigcub, pi
        real        vlrc0, wlrc0
        integer     j, jin

        parameter ( pi = 3.14159265 )

c     ******************************************************************

        ovrlap = .false.
        rcutsq = rcut * rcut
        rminsq = rmin * rmin
        sigsq  = sigma * sigma

c    ** calculate long range corrections **
c    ** note: specific to lennard-jones  **

        sigcub = sigsq * sigma
        sr3    = ( sigma / rcut ) ** 3
        sr9    = sr3 ** 3
        vlrc0  = ( 8.0  / 9.0 ) * pi * sigcub * (     sr9 - 3.0*sr3 )
        wlrc0  = ( 16.0 / 9.0 ) * pi * sigcub * ( 2.0*sr9 - 3.0*sr3 )

c    ** zero accumulators **

        deltv  = 0.0
        deltw  = 0.0

c    ** loop over all atoms  **

        do 100 j = 1, n

c       ** pick active atoms from the array locate **

           jin   = locate(j)

           rxij  = rxi - rx(jin)
           ryij  = ryi - ry(jin)
           rzij  = rzi - rz(jin)

           rxij  = rxij - anint ( rxij )
           ryij  = ryij - anint ( ryij )
           rzij  = rzij - anint ( rzij )

           rijsq = rxij * rxij + ryij * ryij + rzij * rzij

           if ( rijsq .lt. rminsq) then

              ovrlap = .true.
              return

           elseif ( rijsq .lt. rcutsq ) then

              sr2   = sigsq / rijsq
              sr6   = sr2 * sr2 * sr2
              vij   = sr6 * ( sr6 - 1.0 )
              wij   = sr6 * ( sr6 - 0.5 )
              deltv = deltv + vij
              deltw = deltw + wij

           endif

100     continue

        deltv = 4.0  * deltv
        deltw = 48.0 * deltw / 3.0

c    ** add change in long range correction **

        deltv = deltv + ( 2.0 * real ( n ) + 1.0 ) * vlrc0
        deltw = deltw + ( 2.0 * real ( n ) + 1.0 ) * wlrc0

        return
        end



        subroutine potout ( ipull, n, sigma, rcut, deltv, deltw )

        common / block1 / rx, ry, rz
        common / block2 / locate

c    *******************************************************************
c    ** returns the potential energy change when an atom is removed.  **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** integer n                 the number of atoms before removal  **
c    ** integer ipull             the atom to be removed              **
c    ** integer locate(nmax)      array of active atom indices        **
c    ** real    rx(nmax) etc.     the atom positions                  **
c    ** real    deltv             the chance in potential             **
c    ** real    deltw             the change in virial                **
c    ** real    sigma             lj diameter                         **
c    ** real    rcut              cutoff distance for potential       **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** this subroutine is used to calculate the change of energy     **
c    ** during a trial deletion of an atom from the fluid. the long   **
c    ** range corrections is included.                                **
c    *******************************************************************

        integer     nmax
        parameter ( nmax = 500 )
        real        rx(nmax), ry(nmax), rz(nmax)
        real        rcut, sigma, deltv, deltw
        integer     n, ipull, locate(nmax)

        real        rcutsq, sigsq, sr2, sr6, sr3, sr9, rxi, ryi, rzi
        real        rxij, ryij, rzij, rijsq, vij, wij, sigcub, pi
        real        vlrc0, wlrc0
        integer     j, jin

        parameter ( pi = 3.14159265 )

c     ******************************************************************

        rcutsq = rcut * rcut
        sigsq  = sigma * sigma
        sigcub = sigsq * sigma

c    ** calculate long range corrections **
c    ** note: specific to lennard-jones  **

        sr3    = ( sigma / rcut ) ** 3
        sr9    = sr3 ** 3
        vlrc0  = ( 8.0  / 9.0 ) * pi * sigcub * (     sr9 - 3.0*sr3 )
        wlrc0  = ( 16.0 / 9.0 ) * pi * sigcub * ( 2.0*sr9 - 3.0*sr3 )

c    ** zero accumulators **

        deltv  = 0.0
        deltw  = 0.0

        rxi    = rx(ipull)
        ryi    = ry(ipull)
        rzi    = rz(ipull)

c    ** loop over all atoms  except ipull **

        do 100 j = 1, n

c       ** pick active atoms from locate **

           jin = locate(j)

           if ( jin .ne. ipull ) then

              rxij  = rxi - rx(jin)
              ryij  = ryi - ry(jin)
              rzij  = rzi - rz(jin)

              rxij  = rxij - anint ( rxij )
              ryij  = ryij - anint ( ryij )
              rzij  = rzij - anint ( rzij )

              rijsq = rxij * rxij + ryij * ryij + rzij * rzij

              if ( rijsq .lt. rcutsq ) then

                 sr2   = sigsq / rijsq
                 sr6   = sr2 * sr2 * sr2
                 vij   = sr6 * ( sr6 - 1.0 )
                 wij   = sr6 * ( sr6 - 0.5 )
                 deltv = deltv + vij
                 deltw = deltw + wij

              endif

           endif

100     continue

        deltv =  4.0 * deltv
        deltw = 48.0 * deltw / 3.0

c    ** add change in long range correction **

        deltv = deltv + ( 2.0 * real ( n ) - 1.0 ) * vlrc0
        deltw = deltw + ( 2.0 * real ( n ) - 1.0 ) * wlrc0

c    ** change sign of deltv and deltw for a removal **

        deltv = - deltv
        deltw = - deltw

        return
        end



********************************************************************************
** fiche f.14.  algorithm to handle indices in constant mu vt monte carlo     **
** this fortran code is intended to illustrate points made in the text.       **
** to our knowledge it works correctly.  however it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

c    *******************************************************************
c    ** index-handling in grand canonical monte carlo simulation.     **
c    **                                                               **
c    ** routines supplied:                                            **
c    **                                                               **
c    ** subroutine add ( rxnew, rynew, rznew, n )                     **
c    **    adds an atom to the array locate.                          **
c    ** subroutine remove ( nloc, ipull, n )                          **
c    **    removes an atom from the array locate.                     **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** integer n                   number of atoms before trial      **
c    ** integer nmax                maximum number of atoms           **
c    ** integer ipull               index of atom for removal         **
c    ** integer nloc                position of n in locate           **
c    ** integer ntrial              number of atoms during trial      **
c    ** integer locate(nmax)        array of active atom indices      **
c    ** real    rxnew,rynew,rznew   position for addition of an atom  **
c    ** real    rx(nmax), etc.      positions of atoms                **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** routine add is called after a successful trial addition.      **
c    ** routine remove is called after a successful trial removal.    **
c    ** the array locate is updated in each case.                     **
c    *******************************************************************



        subroutine add ( rxnew, rynew, rznew, n )

        common / block1 / rx, ry, rz
        common / block2 / locate

c    *******************************************************************
c    ** subroutine to add an atom to the array locate.                **
c    **                                                               **
c    ** there are n atoms in the simulation before the new addition   **
c    *******************************************************************

        integer     nmax
        parameter ( nmax = 500 )

        real        rx(nmax), ry(nmax), rz(nmax)
        integer     locate(nmax)

        real        rxnew, rynew, rznew
        integer     n, ipull

        integer     inew, ntrial

c    *******************************************************************

        ntrial = n + 1
        inew = locate(ntrial)

        if ( inew .eq. 0 ) then

c       ** atom requires a new number **

           locate(ntrial) = ntrial
           inew           = ntrial

        endif

c    ** fit new atom into the array **

        rx(inew) = rxnew
        ry(inew) = rynew
        rz(inew) = rznew

        return
        end



        subroutine remove ( nloc, ipull, n )

        common / block1 / rx, ry, rz
        common / block2 / locate

c    *******************************************************************
c    ** subroutine to remove an atom from the array locate.           **
c    **                                                               **
c    ** there are n atoms in the simulation before the removal.       **
c    ** element ipull of locate is to be destroyed.                   **
c    *******************************************************************

        integer     nmax
        parameter ( nmax = 500 )

        real        rx(nmax), ry(nmax), rz(nmax)
        integer     locate(nmax)
        integer     n, ipull, nloc

        integer     k

c    *******************************************************************

        if ( nloc .lt. n ) then

c       ** close up the array locate after the removal **

           do 10 k = nloc + 1, n

              locate(k - 1) = locate(k)

10         continue

c       ** place the ghost atom ipull just outside the active **
c       ** range of the array locate for future use           **

           locate(n) = ipull

        endif

        return
        end



