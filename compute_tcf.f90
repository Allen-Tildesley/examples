! compute_tcf.f90
PROGRAM compute_tcf

        COMMON / block1 / storx, story, storz
        COMMON / block2 / vx, vy, vz
        COMMON / block3 / vacf, anorm

c    *******************************************************************
c    ** calculation of time correlation functions.                    **
c    **                                                               **
c    ** this PROGRAM analyzes DATA to calculate a time correlation    **
c    ** FUNCTION in one sweep (low memory requirement). in this       **
c    ** example the velocity auto-correlation FUNCTION is calculated. **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** INTEGER  n                  number of atoms                   **
c    ** INTEGER  nstep              number of steps on the tape       **
c    ** INTEGER  ior                interval for time origins         **
c    ** INTEGER  nt                 correlation length, including t=0 **
c    ** INTEGER  ntimor             number of time origins            **
c    ** INTEGER  nlabel             label for step (1,2,3...nstep)    **
c    ** REAL     vx(n),vy(n),vz(n)  velocities                        **
c    ** REAL     vacf(nt)           the correlation FUNCTION          **
c    ** nstep and nt should be multiples of ior.                      **
c    **                                                               **
c    ** routines referenced:                                          **
c    **                                                               **
c    ** SUBROUTINE store ( j1 )                                       **
c    **    routine to store the DATA for correlation                  **
c    ** SUBROUTINE corr ( j1, j2, it )                                **
c    **    routine to correlate the stored time origins               **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** DATA in file dfile on fortran unit dunit.                     **
c    ** results in file rfile on fortran unit runit.                  **
c    *******************************************************************

        INTEGER     n, nstep, ior, nt, ndim, dunit, runit, ntimor
        INTEGER     fullup
        PARAMETER ( n = 256, nstep = 1000, ior = 4, nt = 200  )
        PARAMETER ( dunit = 10, runit = 11                    )
        PARAMETER ( ndim = nt / ior + 1, ntimor = nstep / ior )
        PARAMETER ( fullup = ndim - 1                         )

        REAL        vx(n), vy(n), vz(n)
        REAL        storx(ndim,n), story(ndim,n), storz(ndim,n)
        REAL        vacf(nt), anorm(nt)
        INTEGER     s(ntimor), tm(ntimor)
        INTEGER     ts, tss, l, nincor, k, ja, ib, in, ia, jo, i
        INTEGER     nlabel
        CHARACTER   dfile * 30
        CHARACTER   rfile * 30

c    *******************************************************************

        WRITE(*,'('' **** program tcorr ****                       '')')
        WRITE(*,'('' calculation of time correlation functions     '')')

c    ** READ in file names **

        WRITE(*,'('' enter data file name                          '')')
        READ (*,'(a)') dfile
        WRITE (*,'('' enter results file name                      '')')
        READ (*,'(a)') rfile

c    ** initialize counters **

        nincor = fullup
        ja = 1
        ia = 1
        ib = 1

c    ** zero arrays **

        DO 5 i = 1, nt

           vacf(i)  = 0.0
           anorm(i) = 0.0

5       CONTINUE

c    ** OPEN DATA file and results file **

        OPEN ( unit = dunit, file = dfile, access = 'sequential',
     :         status = 'old', form = 'unformatted' )

        OPEN ( unit = runit, file = rfile, status = 'new' )

c   ** calculation begins **

        DO 40 l = 1, ntimor

           ja   = ja + 1
           s(l) = ja - 1

           READ ( dunit ) nlabel, vx, vy, vz

           tm(l) = nlabel

c       ** store step as a time origin **

           CALL store ( ja )

c       ** correlate the origins in store **

           DO 10 in = ia, l

              tss = tm(l) - tm(in)
              ts  = tss + 1
              jo  = s(in) + 1
              CALL corr ( jo, ja, ts )

10         CONTINUE

c       ** READ in DATA between time origins. this can  **
c       ** be conveniently stored in element 1 of the   **
c       ** arrays storx etc. and can THEN be correlated **
c       ** WITH the time origins.                       **

           DO 30 k = 1, ior - 1

              READ ( dunit ) nlabel, vx, vy, vz

              CALL store ( 1 )

              DO 20 in = ia, l

                 tss = nlabel - tm(in)
                 ts  = tss + 1
                 jo  = s(in) + 1
                 CALL corr ( jo, 1, ts )

20            CONTINUE

30         CONTINUE

           IF ( l .GE. fullup ) THEN

              IF ( l .EQ. nincor ) THEN

                 nincor = nincor + fullup
                 ja     = 1

              ENDIF

              ia = ia + 1

           ENDIF

40      CONTINUE

        CLOSE ( unit = dunit )

c    ** normalise correlation functions **

        vacf(1) = vacf(1) / anorm(1) / REAL ( n )

        DO 50 i = 2, nt

           vacf(i) = vacf(i) / anorm(i) / REAL ( n ) / vacf(1)

50      CONTINUE

        WRITE ( runit, '('' velocity acf '')')
        WRITE ( runit, '(i6,e15.6)') ( i, vacf(i), i = 1, nt )

        CLOSE ( runit )

        STOP
        END



        SUBROUTINE store ( j1 )

        COMMON/ block1 / storx, story, storz
        COMMON/ block2 / vx, vy, vz

c    *******************************************************************
c    ** SUBROUTINE to store time origins                              **
c    *******************************************************************

        INTEGER     j1
        INTEGER     n, nt, ior, ndim
        PARAMETER ( n = 256, nt = 200, ior = 4 )
        PARAMETER ( ndim = nt / ior + 1        )

        REAL        storx(ndim,n), story(ndim,n), storz(ndim,n)
        REAL        vx(n), vy(n), vz(n)
        INTEGER     i


        DO 10 i = 1, n

           storx(j1,i) = vx(i)
           story(j1,i) = vy(i)
           storz(j1,i) = vz(i)

10      CONTINUE

        RETURN
        END



        SUBROUTINE corr ( j1, j2, it )

        COMMON/ block1 / storx, story, storz
        COMMON/ block3 / vacf, anorm

c    *******************************************************************
c    ** SUBROUTINE to correlate time origins                          **
c    *******************************************************************

        INTEGER     j1, j2, it
        INTEGER     n, nt, ior, ndim
        PARAMETER ( n = 256, nt = 200, ior = 4 )
        PARAMETER ( ndim = nt / ior + 1        )

        REAL        storx(ndim,n), story(ndim,n), storz(ndim,n)
        REAL        vacf(nt), anorm(nt)
        INTEGER     i

c    *******************************************************************

        DO 10 i = 1, n

           vacf(it) = vacf(it) + storx(j1,i) * storx(j2,i)
     :                         + story(j1,i) * story(j2,i)
     :                         + storz(j1,i) * storz(j2,i)

10      CONTINUE

        anorm(it) = anorm(it) + 1.0

        RETURN
        END



