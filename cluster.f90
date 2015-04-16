! cluster.f90
PROGRAM cluster
  ! placeholder for f90 subroutine

  USE utility_module, ONLY : read_cnf_atoms
  IMPLICIT NONE

  ! Program to identify atom clusters in a configuration
  ! Defines a cluster by a critical separation
  ! Produces a linked list of clusters
  ! Works in units where box = 1 
  ! Reference:
  ! Stoddard J Comp Phys, 27, 291, 1977.

  INTEGER n
  REAL, DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,n)
  INTEGER, DIMENSION(:), ALLOCATABLE :: list


  REAL     ::   rcl
  INTEGER  ::   it, count

  REAL      ::  rcl_sq, rjk_sq, 
  REAL, DIMENSION(3) :: rjk
  INTEGER  ::   i, j, k

  NAMELIST /cluster_parameters/ rcl

  READ(*,cluster_parameters)

  rcl_sq = rcl * rcl

  CALL read_cnf_atoms ( 'cluster.cnf', n, box )
  WRITE(*,'(''Number of particles'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  ALLOCATE ( r(3,n), list(n) )
  CALL read_cnf_atoms ( 'cluster.cnf', n, box, r )
  r(:,:) = r(:,:) / box ! Convert to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  list(:) = [ (i,i=1,n) ] ! set up the sorting array

  DO i = 1, n - 1 ! Begin outer loop

     IF ( i == list(i) ) THEN

        DO j = i + 1, n ! Begin first inner loop

           IF ( list(j) == j ) THEN
              IF ( in_range ( r(:,i), r(:,j), rcl_sq ) ) THEN
                 list(j) = list(i)
                 list(i) = j
              END IF
           END IF

        END DO ! End first inner loop

        j = list(i)

        DO
           IF ( j == i ) EXIT

           DO k = i + 1, n ! Begin second inner loop

              IF ( list(k) == k ) THEN
                 IF ( in_range ( r(:,j), r(:,k), rcl_sq ) ) THEN
                    list(k) = list(j)
                    list(j) = k
                 END IF
              END IF

           END DO ! End second inner loop

           j = list(j)

        END DO

     END IF

  END DO ! End outer loop

  count = 1
  j = list(it)

  DO 
     IF ( j == it ) EXIT
     count = count + 1
     j = list(j)
  END DO

CONTAINS

  FUNCTION in_range ( rj, rk, rcl_sq )
    LOGICAL                        :: in_range
    REAL, DIMENSION(3), INTENT(in) :: rj, rk
    REAL,               INTENT(in) :: rcl_sq

    REAL, DIMENSION(3) :: rjk
    REAL               :: rjk_sq

    rjk(:) = rj(:) - rk(:)
    rjk(:) = rjk(:) - ANINT ( rjk(:) )
    rjk_sq = SUM ( rjk**2 )

    in_range = ( rjk_sq <= rcl_sq )

  END FUNCTION in_range

END PROGRAM cluster



