! cluster.f90
! Identify atom clusters in a configuration
PROGRAM cluster

  USE utility_module, ONLY : read_cnf_atoms
  IMPLICIT NONE

  ! Defines a cluster by a critical separation
  ! Produces a set of circular linked lists of clusters
  ! Works in units where box = 1 
  ! Reference:
  ! Stoddard J Comp Phys, 27, 291, 1977
  ! This simple algorithm does not scale well to large n

  INTEGER n
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: list ! (n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: done ! (n)

  REAL     ::   rcl, rcl_sq, box
  INTEGER  ::   count, cluster_id
  INTEGER  ::   i, j, k

  NAMELIST /cluster_parameters/ rcl

  READ(*,cluster_parameters)
  WRITE(*,'(''Cluster critical separation'',f15.5)') rcl

  rcl_sq = rcl * rcl

  CALL read_cnf_atoms ( 'cluster.cnfinp', n, box )
  WRITE(*,'(''Number of particles'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box

  ALLOCATE ( r(3,n), list(n), done(n) )

  CALL read_cnf_atoms ( 'cluster.cnf', n, box, r )
  r(:,:) = r(:,:) / box ! Convert to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  list(:) = [ (i,i=1,n) ] ! set up the list

  DO i = 1, n - 1 ! Begin outer loop

     IF ( i == list(i) ) THEN

        j = i

        DO ! Begin inner loop

           DO k = i + 1, n ! Begin innermost loop

              IF ( list(k) == k ) THEN
                 IF ( in_range ( j, k ) ) list([k,j]) = list([j,k])
              END IF

           END DO ! End innermost loop

           j = list(j)
           IF ( j == i ) EXIT

        END DO ! End inner loop

     END IF

  END DO ! End outer loop

  ! for diagnostic purposes, print out the cluster membership
  ! no particular sorting (e.g. by size)
  
  done = 0 ! zero array
  cluster_id = 0

  DO ! Begin loop over remaining clusters
     IF ( ALL ( done > 0 ) ) EXIT
     i = MINLOC ( done, dim = 1 ) ! find first zero
     cluster_id = cluster_id + 1
     WRITE(*,'(''Cluster '',i5,'' members are:'')') cluster_id
     j = i
     done(j) = cluster_id
     count = 1
     WRITE(*,'(i5)',advance='no') j
     
     DO ! Begin loop to find other members of cluster
        j = list(j)
        IF ( j == i ) EXIT ! link list has returned to start
        done(j) = cluster_id
        count = count + 1
        WRITE(*,'(i5)',advance='no') j
     END DO ! End loop to find other members of cluster

     WRITE(*,'(/,''Cluster '',i5,'' contains '',i10,'' members'')') cluster_id, count
  END DO ! End loop over remaining clusters

  DEALLOCATE ( r, list, done )

CONTAINS

  FUNCTION in_range ( j, k )
    LOGICAL             :: in_range
    integer, INTENT(in) :: j, k

    REAL, DIMENSION(3) :: rjk
    REAL               :: rjk_sq

    rjk(:) = r(:,j) - r(:,k)
    rjk(:) = rjk(:) - ANINT ( rjk(:) )
    rjk_sq = SUM ( rjk**2 )

    in_range = ( rjk_sq <= rcl_sq )

  END FUNCTION in_range

END PROGRAM cluster



