! cluster.f90
! Identify atom clusters in a configuration
PROGRAM cluster

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             !
  ! to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     !
  ! published by Oxford University Press ("the publishers").                                       !
  !                                                                                                !
  ! LICENCE                                                                                        !
  ! Creative Commons CC0 Public Domain Dedication.                                                 !
  ! To the extent possible under law, the authors have dedicated all copyright and related         !
  ! and neighboring rights to this software to the PUBLIC domain worldwide.                        !
  ! This software is distributed without any warranty.                                             !
  ! You should have received a copy of the CC0 Public Domain Dedication along with this software.  !
  ! If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               !
  !                                                                                                !
  ! DISCLAIMER                                                                                     !
  ! The authors and publishers make no warranties about the software, and disclaim liability       !
  ! for all uses of the software, to the fullest extent permitted by applicable law.               !
  ! The authors and publishers do not recommend use of this software for any purpose.              !
  ! It is made freely available, solely to clarify points made in the text. When using or citing   !
  ! the software, you should not imply endorsement by the authors or publishers.                   !
  !------------------------------------------------------------------------------------------------!

  ! Reads an atomic configuration with periodic boundary conditions from cluster.inp
  ! Defines a cluster by a critical separation r_cl
  ! Value of r_cl read from standard input using a namelist nml
  ! Leave namelist empty to accept supplied default

  ! Produces a set of circular linked lists of clusters
  ! Input data in atomic (e.g. LJ sigma) units
  ! Program works in the same units
  ! Reference: SD Stoddard J Comp Phys, 27, 291 (1978)
  ! This simple algorithm does not scale well to large N

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms

  IMPLICIT NONE

  INTEGER                              :: n    ! Number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r    ! Positions (3,n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: list ! Linked list array (n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: done ! Indicates assignment to cluster (n)

  CHARACTER(len=11), PARAMETER :: filename = 'cluster.inp'
  REAL                         :: r_cl, r_cl_sq, box
  INTEGER                      :: ioerr, count, cluster_id
  INTEGER                      :: i, j, k

  NAMELIST /nml/ r_cl

  r_cl = 1.5 ! default value
  
  READ ( unit=input_unit, nml=nml, iostat=ioerr ) ! namelist input
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in cluster'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Cluster separation distance', r_cl

  CALL read_cnf_atoms ( filename, n, box ) ! First call to obtain n
  WRITE( unit=output_unit, fmt='(a,t40,i15)'  ) 'Number of particles',  n
  WRITE( unit=output_unit, fmt='(a,t40,f15.6)') 'Box (in sigma units)', box

  ALLOCATE ( r(3,n), list(n), done(n) )

  CALL read_cnf_atoms ( filename, n, box, r ) ! Second call to read in configuration

  r(:,:)  = r(:,:) - ANINT ( r(:,:) / box ) * box ! Apply periodic boundaries

  r_cl_sq = r_cl**2 ! used in in_range function

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

     IF ( ALL ( done > 0 ) ) EXIT ! Loop until all done

     i = MINLOC ( done, dim = 1 ) ! Find first zero (FINDLOC is not implemented in gfortran at the time of writing)
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
    LOGICAL             :: in_range ! Returns indicator of whether pair is in range or not
    INTEGER, INTENT(in) :: j, k     ! Supplied pair of atom indices

    REAL, DIMENSION(3) :: rjk
    REAL               :: rjk_sq

    rjk(:) = r(:,j) - r(:,k)                       ! Separation vector
    rjk(:) = rjk(:) - ANINT ( rjk(:) / box ) * box ! Periodic boundary conditions
    rjk_sq = SUM ( rjk**2 )                        ! Squared separation

    in_range = ( rjk_sq <= r_cl_sq ) ! Determines whether pair is in range

  END FUNCTION in_range

END PROGRAM cluster



