! cluster.f90
! Identify atom clusters in a configuration
PROGRAM cluster

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms

  IMPLICIT NONE

  ! Reads an atomic configuration with periodic boundary conditions from inp.cnf
  ! Defines a cluster by a critical separation rcl
  ! Value of rcl read from standard input using a namelist nml
  ! Leave namelist empty to accept supplied default

  ! Produces a set of circular linked lists of clusters
  ! Input data in atomic (e.g. LJ sigma) units
  ! Program works in box = 1 units
  ! Reference: SD Stoddard J Comp Phys, 27, 291 (1978)
  ! This simple algorithm does not scale well to large N

  INTEGER                              :: n    ! Number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r    ! Positions (3,n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: list ! Linked list array (n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: done ! Indicates assignment to cluster (n)

  CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp'
  REAL                        :: rcl, rcl_sq, box
  INTEGER                     :: ioerr, count, cluster_id
  INTEGER                     :: i, j, k

  NAMELIST /nml/ rcl

  rcl = 1.5 ! default value
  
  READ ( unit=input_unit, nml=nml, iostat=ioerr ) ! namelist input
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in cluster'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Cluster separation distance', rcl

  CALL read_cnf_atoms ( filename, n, box )
  WRITE( unit=output_unit, fmt='(a,t40,i15)'  ) 'Number of particles', n
  WRITE( unit=output_unit, fmt='(a,t40,f15.5)') 'Box (in sigma units)', box

  ALLOCATE ( r(3,n), list(n), done(n) )

  CALL read_cnf_atoms ( filename, n, box, r )

  ! Convert to box units and apply periodic boundaries
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) )
  rcl    = rcl / box
  rcl_sq = rcl * rcl

  list(:) = [ (i,i=1,n) ] ! Set up the list

  DO i = 1, n - 1 ! Begin outer loop

     IF ( i == list(i) ) THEN

        j = i

        DO ! Begin inner loop

           DO k = i + 1, n ! Begin innermost loop

              IF ( list(k) == k ) THEN
                 IF ( in_range ( j, k ) ) list([k,j]) = list([j,k]) ! swap elements
              END IF

           END DO ! End innermost loop

           j = list(j)
           IF ( j == i ) EXIT

        END DO ! End inner loop

     END IF

  END DO ! End outer loop

  ! For diagnostic purposes, print out the cluster membership
  ! no particular sorting (e.g. by size)

  done(:)    = 0
  cluster_id = 0

  WRITE ( unit=output_unit, fmt='(a)' ) 'Cluster Members .....'
  DO ! Begin loop over remaining clusters

     IF ( ALL ( done > 0 ) ) EXIT

     i = MINLOC ( done, dim = 1 ) ! find first zero (FINDLOC is not implemented in gfortran at the time of writing)
     cluster_id = cluster_id + 1
     WRITE ( unit=output_unit, fmt='(a,i5,a)', advance='no' ) 'Cluster ', cluster_id, ' = '
     j = i
     done(j) = cluster_id
     WRITE ( unit=output_unit, fmt='(i5)', advance='no') j

     DO ! Begin loop to find other members of cluster
        j = list(j)
        IF ( j == i ) EXIT ! link list has returned to start
        done(j) = cluster_id
        WRITE ( unit=output_unit, fmt='(i5)', advance='no') j
     END DO ! End loop to find other members of cluster
     WRITE ( unit=output_unit, fmt='(1x)' )

  END DO ! End loop over remaining clusters

  ! Count cluster members
  WRITE ( unit=output_unit, fmt='(/,a)' ) 'Cluster Count'
  DO i = 1, cluster_id
     WRITE ( unit=output_unit, fmt='(i7,1x,i5)' ) i, COUNT ( done == i )
  END DO

  DEALLOCATE ( r, list, done )

CONTAINS

  FUNCTION in_range ( j, k )
    IMPLICIT NONE
    LOGICAL             :: in_range ! Returns indicator of whether pair is in-range or not
    INTEGER, INTENT(in) :: j, k     ! Supplied pair of atom indices

    REAL, DIMENSION(3) :: rjk
    REAL               :: rjk_sq

    rjk(:) = r(:,j) - r(:,k)
    rjk(:) = rjk(:) - ANINT ( rjk(:) )
    rjk_sq = SUM ( rjk**2 )

    in_range = ( rjk_sq <= rcl_sq )

  END FUNCTION in_range

END PROGRAM cluster



