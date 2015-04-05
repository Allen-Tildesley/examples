! mc_muvt_lj_module.f90  (used by mc_muvt_lj.f90)
! Monte Carlo simulation, constant-muVT ensemble, Lennard-Jones atoms
module mc_muvt_lj_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nmax, n, r, locate

  INTEGER                              :: nmax, n ! max and actual number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r       ! positions (3,nmax)
  integer, DIMENSION(:),   ALLOCATABLE :: locate  ! array of active atom indices (nmax)

!*****************************************************************
! attempted creations and destructions in grand canonical mc.   **
!                                                               **
! these routines allow for a trial destruction or creation in a **
! grand canonical monte carlo program.                          **
!                                                               **
! principal variables:                                          **
!                                                               **
! integer n                   number of atoms before the trial  **
! integer ntrial              number of atoms during the trial  **
! real    rxnew,rynew,rznew   position for addition of atom     **
! real    v                   potential energy + lrc            **
! real    w                   virial + lrc                      **
! real    deltv               change in energy                  **
! real    deltw               change in virial                  **
! real    temp                reduced temperature               **
! real    z                   absolute activity coefficient     **
! real    sigma               lennard jones diameter            **
! real    rcut                reduced cutoff distance           **
! real    rmin                reduced minimum separation        **
! logical ovrlap              true for substantial atom overlap **
! logical create              true for an accepted creation     **
! logical ghost               true for an accepted destruction  **
!                                                               **
! routines supplied:                                            **
!                                                               **
! subroutine in ( temp, z, sigma, rcut, n, v, w, create )       **
!    performs a trial creation                                  **
! subroutine out ( temp, z, sigma, rcut, n, v, w, ghost )       **
!    performs a trial destruction                               **
! subroutine potin ( rxnew, rynew, rznew, n, sigma, rcut, rmin, **
! :                  deltv, deltw, ovrlap )                     **
!    calculates the potential energy change on creation         **
! subroutine potout ( ipull, n, sigma, rcut, deltv, deltw )     **
!    calculates the potential energy change on destruction      **
!                                                               **
! routines referenced:                                          **
!                                                               **
! real function ranf ( dummy )  (given in f.11)                 **
!    returns a uniform random variate on zero to one            **
! subroutine add ( rxnew, rynew, rznew, n )                     **
!    updates locate after addition (given in f.14)              **
! subroutine remove ( nloc, ipull, n )                          **
!    updates locate after removal (given in f.14)               **
!                                                               **
! usage:                                                        **
!                                                               **
! routines in and out should be called with equal probability   **
! in a grand canonical monte carlo simulation. if a trial       **
! creation is accepted then create is set to true. if a trial   **
! destruction is accepted then ghost is set to true. the        **
! routines are written for lennard-jones atoms. the box is of   **
! unit length, all distances are scaled to the box length.      **
! trial inputs which result in a separation of less than        **
! 0.5*sigma are rejected. the long-range corrections are        **
! included in v and w. all accumulators are updated in the main **
! part of the program which is not given here.                  **
!*****************************************************************

        subroutine in ( temp, z, sigma, rcut, n, v, w, create )

        common / block1 / rx, ry, rz
        common / block2 / locate

!*****************************************************************
! routine to attempt a trial creation                           **
!                                                               **
! principal variables:                                          **
!                                                               **
! real    temp            temperature                           **
! real    z               absolute activity                     **
! real    sigma           lennard-jones diameter                **
! real    rcut            cut-off distance                      **
! real    v               potential energy                      **
! real    w               virial                                **
! integer n               number of atoms before trial creation **
! logical create          true for a successful creation        **
!*****************************************************************

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

!*****************************************************************

        create = .false.
        beta   = 1.0 / temp
        rmin   = 0.5 * sigma
        ntrial = n + 1

        if ( ntrial .ge. nmax ) stop 'maximum number of atoms in box'

! generate the position of the trial atom **

        rxnew  = ranf ( dummy ) - 0.5
        rynew  = ranf ( dummy ) - 0.5
        rznew  = ranf ( dummy ) - 0.5

! calculate energy change on addition **

        call potin ( rxnew, rynew, rznew, n, sigma, rcut, rmin, deltv,
     :               deltw, ovrlap )

! check for acceptance **

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

!*****************************************************************
! routine to attempt a trial destruction                        **
!                                                               **
! principal variables:                                          **
!                                                               **
! real    temp         temperature                              **
! real    z            absolute activity                        **
! real    sigma        lennard-jones diameter                   **
! real    rcut         cut-off distance                         **
! real    v            potential energy                         **
! real    w            virial                                   **
! integer n            number of atoms before trial destruction **
! logical ghost        true for a successful destruction        **
!*****************************************************************

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

!*****************************************************************

        ghost  = .false.
        beta   = 1.0 / temp
        ntrial = n - 1

        if ( ntrial .eq. 1 ) stop 'only one atom remains'

! pick a random element from the active part of locate **

        nloc  = int ( real ( ntrial ) * ranf ( dummy ) ) + 1
        ipull = locate(nloc)

! calculate energy change on removal of atom ipull **

        call potout ( ipull, n, sigma, rcut, deltv, deltw )

! check for acceptance **

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

!*****************************************************************
! returns the potential energy change on adding an atom.        **
!                                                               **
! principal variables:                                          **
!                                                               **
! integer n                 the number of atoms before addition **
! real    rxi,ryi,rzi       the coordinates of the added atom   **
! real    deltv             the chance in potential             **
! real    deltw             the change in virial                **
! real    sigma             lj diameter                         **
! real    rcut              cutoff distance for potential       **
! real    rmin              minimum allowed approach of atoms   **
! logical ovrlap            true for substantial atom overlap   **
!                                                               **
! usage:                                                        **
!                                                               **
! this subroutine is used to calculate the change of energy     **
! during a trial addition of an atom to the fluid. the long     **
! range corrections is included.                                **
!*****************************************************************

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

! calculate long range corrections **
! note: specific to lennard-jones  **

        sigcub = sigsq * sigma
        sr3    = ( sigma / rcut ) ** 3
        sr9    = sr3 ** 3
        vlrc0  = ( 8.0  / 9.0 ) * pi * sigcub * (     sr9 - 3.0*sr3 )
        wlrc0  = ( 16.0 / 9.0 ) * pi * sigcub * ( 2.0*sr9 - 3.0*sr3 )

! zero accumulators **

        deltv  = 0.0
        deltw  = 0.0

! loop over all atoms  **

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

! add change in long range correction **

        deltv = deltv + ( 2.0 * real ( n ) + 1.0 ) * vlrc0
        deltw = deltw + ( 2.0 * real ( n ) + 1.0 ) * wlrc0

        return
        end



        subroutine potout ( ipull, n, sigma, rcut, deltv, deltw )

        common / block1 / rx, ry, rz
        common / block2 / locate

!*****************************************************************
! returns the potential energy change when an atom is removed.  **
!                                                               **
! principal variables:                                          **
!                                                               **
! integer n                 the number of atoms before removal  **
! integer ipull             the atom to be removed              **
! integer locate(nmax)      array of active atom indices        **
! real    rx(nmax) etc.     the atom positions                  **
! real    deltv             the chance in potential             **
! real    deltw             the change in virial                **
! real    sigma             lj diameter                         **
! real    rcut              cutoff distance for potential       **
!                                                               **
! usage:                                                        **
!                                                               **
! this subroutine is used to calculate the change of energy     **
! during a trial deletion of an atom from the fluid. the long   **
! range corrections is included.                                **
!*****************************************************************

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

! calculate long range corrections **
! note: specific to lennard-jones  **

        sr3    = ( sigma / rcut ) ** 3
        sr9    = sr3 ** 3
        vlrc0  = ( 8.0  / 9.0 ) * pi * sigcub * (     sr9 - 3.0*sr3 )
        wlrc0  = ( 16.0 / 9.0 ) * pi * sigcub * ( 2.0*sr9 - 3.0*sr3 )

! zero accumulators **

        deltv  = 0.0
        deltw  = 0.0

        rxi    = rx(ipull)
        ryi    = ry(ipull)
        rzi    = rz(ipull)

! loop over all atoms  except ipull **

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

! add change in long range correction **

        deltv = deltv + ( 2.0 * real ( n ) - 1.0 ) * vlrc0
        deltw = deltw + ( 2.0 * real ( n ) - 1.0 ) * wlrc0

! change sign of deltv and deltw for a removal **

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

!*****************************************************************
! index-handling in grand canonical monte carlo simulation.     **
!                                                               **
! routines supplied:                                            **
!                                                               **
! subroutine add ( rxnew, rynew, rznew, n )                     **
!    adds an atom to the array locate.                          **
! subroutine remove ( nloc, ipull, n )                          **
!    removes an atom from the array locate.                     **
!                                                               **
! principal variables:                                          **
!                                                               **
! integer n                   number of atoms before trial      **
! integer nmax                maximum number of atoms           **
! integer ipull               index of atom for removal         **
! integer nloc                position of n in locate           **
! integer ntrial              number of atoms during trial      **
! integer locate(nmax)        array of active atom indices      **
! real    rxnew,rynew,rznew   position for addition of an atom  **
! real    rx(nmax), etc.      positions of atoms                **
!                                                               **
! usage:                                                        **
!                                                               **
! routine add is called after a successful trial addition.      **
! routine remove is called after a successful trial removal.    **
! the array locate is updated in each case.                     **
!*****************************************************************



        subroutine add ( rxnew, rynew, rznew, n )

        common / block1 / rx, ry, rz
        common / block2 / locate

!*****************************************************************
! subroutine to add an atom to the array locate.                **
!                                                               **
! there are n atoms in the simulation before the new addition   **
!*****************************************************************

        integer     nmax
        parameter ( nmax = 500 )

        real        rx(nmax), ry(nmax), rz(nmax)
        integer     locate(nmax)

        real        rxnew, rynew, rznew
        integer     n, ipull

        integer     inew, ntrial

!*****************************************************************

        ntrial = n + 1
        inew = locate(ntrial)

        if ( inew .eq. 0 ) then

c       ** atom requires a new number **

           locate(ntrial) = ntrial
           inew           = ntrial

        endif

! fit new atom into the array **

        rx(inew) = rxnew
        ry(inew) = rynew
        rz(inew) = rznew

        return
        end



        subroutine remove ( nloc, ipull, n )

        common / block1 / rx, ry, rz
        common / block2 / locate

!*****************************************************************
! subroutine to remove an atom from the array locate.           **
!                                                               **
! there are n atoms in the simulation before the removal.       **
! element ipull of locate is to be destroyed.                   **
!*****************************************************************

        integer     nmax
        parameter ( nmax = 500 )

        real        rx(nmax), ry(nmax), rz(nmax)
        integer     locate(nmax)
        integer     n, ipull, nloc

        integer     k

!*****************************************************************

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



SUBROUTINE resize_array
INTEGER,DIMENSION(:),ALLOCATABLE :: tmp_arr

ALLOCATE(tmp_arr(2*SIZE(array)))
tmp_arr(1:SIZE(array))=array
DEALLOCATE(array)
ALLOCATE(array(size(tmp_arr)))
array=tmp_arr

ENDSUBROUTINE resize_array
