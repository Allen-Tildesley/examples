! cluster.f90
program cluster
! placeholder for f90 subroutine

  use utility_module, only : read_cnf_atoms
  implicit none

  ! Program to identify atom clusters in a configuration
  ! Defines a cluster by a critical separation
  ! Produces a linked list of clusters
  ! Works in units where box = 1 
! Reference:
! Stoddard J Comp Phys, 27, 291, 1977.

  integer n
  real, dimension(:,:), allocatable :: r ! positions (3,n)
  integer, dimension(:), allocatable :: l

  
  REAL     ::   rcl
        INTEGER     it, nit

        REAL        rclsq, rxjk, ryjk, rzjk
        REAL        rjksq, rxj, ryj, rzj
        INTEGER     i, j, k, lk, lit

        rclsq = rcl * rcl

! set up the sorting array **

        DO 10 i = 1, n

           l(i) = i

10      CONTINUE

! sort the clusters **

        DO 50 i = 1, n - 1

           IF ( i .EQ. l(i) ) THEN

              j   = i
              rxj = rx(j)
              ryj = ry(j)
              rzj = rz(j)

              DO 20 k = i + 1, n

                 lk = l(k)

                 IF ( lk .EQ. k ) THEN

                    rxjk  = rxj - rx(k)
                    ryjk  = ryj - ry(k)
                    rzjk  = rzj - rz(k)
                    rxjk  = rxjk - ANINT ( rxjk )
                    ryjk  = ryjk - ANINT ( ryjk )
                    rzjk  = rzjk - ANINT ( rzjk )
                    rjksq = rxjk * rxjk + ryjk * ryjk + rzjk * rzjk

                    IF ( rjksq .LE. rclsq ) THEN

                       l(k) = l(j)
                       l(j) = lk

                    ENDIF

                 ENDIF

20            CONTINUE

              j   = l(j)
              rxj = rx(j)
              ryj = ry(j)
              rzj = rz(j)

30            IF ( j .NE. i ) THEN

                 DO 40 k = i + 1, n

                    lk = l(k)

                    IF ( lk .EQ. k ) THEN

                       rxjk  = rxj - rx(k)
                       ryjk  = ryj - ry(k)
                       rzjk  = rzj - rz(k)
                       rxjk  = rxjk - ANINT ( rxjk )
                       ryjk  = ryjk - ANINT ( ryjk )
                       rzjk  = rzjk - ANINT ( rzjk )
                       rjksq = rxjk * rxjk + ryjk * ryjk + rzjk * rzjk

                       IF ( rjksq .LE. rclsq ) THEN

                          l(k) = l(j)
                          l(j) = lk

                       ENDIF

                    ENDIF

40               CONTINUE

                 j   = l(j)
                 rxj = rx(j)
                 ryj = ry(j)
                 rzj = rz(j)

                 go to 30

              ENDIF

           ENDIF

50      CONTINUE

c   **  count the number in a cluster containing atom it **

        nit = 1
        lit = l(it)

60      IF ( lit .NE. it ) THEN

           nit = nit + 1
           lit = l(lit)

           go to 60

        ENDIF

        RETURN
        END



