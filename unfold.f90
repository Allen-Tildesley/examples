! unfold.f90
SUBROUTINE unfold ( n, rx, ry, rz, rx0, ry0, rz0 )

c    *******************************************************************
c    ** SUBROUTINE to fold trajectories in periodic boundaries.       **
c    **                                                               **
c    ** the unfolding routine undoes the effect of folding.           **
c    ** again we take the unit cube as an example.                    **
c    ** this routine requires that coordinates from successive steps  **
c    ** be supplied and assumes that natural movement across half a   **
c    ** box length in one step will never occur.                      **
c    ** the resulting coordinates would be suitable, for example, to  **
c    ** calculate the diffusion coefficient by einstein's relation.   **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** integer n                    number of molecules              **
c    ** real    rx(n),ry(n),rz(n)    molecular positions at time t    **
c    ** real    rx0(n),ry0(n),ry0(n) molecular positions at time t-dt **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** (1)  read initial coordinates into rx0,ry0,rz0                **
c    ** (2)  write rx0,ry0,rz0 to output file (or use immediately)    **
c    ** (3)  read next step coordinates into rx,ry,rz                 **
c    ** (4)  call unfold ( n, rx, ry, rz, rx0, ry0, rz0 )             **
c    ** (5)  write rx,ry,rz to output file (or use immediately)       **
c    ** (6)  set rx0(i)=rx(i), ry0(i)=ry(i), rz0(i)=rz(i), i = 1,n    **
c    ** (7)  unless data exhausted, go to (3)                         **
c    *******************************************************************

        integer n
        real    rx(n), ry(n), rz(n), rx0(n), ry0(n), rz0(n)

        integer i
        real    dx, dy, dz

c    *******************************************************************

        do 100 i = 1, n

           dx = rx(i) - rx0(i)
           dy = ry(i) - ry0(i)
           dz = rz(i) - rz0(i)
           dx = dx - anint ( dx )
           dy = dy - anint ( dy )
           dz = dz - anint ( dz )
           rx(i) = rx0(i) + dx
           ry(i) = ry0(i) + dy
           rz(i) = rz0(i) + dz

100     continue

        return
        end



