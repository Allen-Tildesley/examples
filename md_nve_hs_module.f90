! md_nve_hs_module.f90
! Collisions and overlap for MD of hard spheres
MODULE md_module

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

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  public :: update, overlap, collide

  ! Public data
  INTEGER,                              PUBLIC :: n       ! Number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r       ! Positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v       ! Velocities (3,n)
  REAL,    DIMENSION(:),   ALLOCATABLE, PUBLIC :: coltime ! Time to next collision (n)
  INTEGER, DIMENSION(:),   ALLOCATABLE, PUBLIC :: partner ! Collision partner (n)

  INTEGER, PARAMETER, PUBLIC :: lt = -1, gt = 1 ! Options for j_range

CONTAINS

  SUBROUTINE introduction
    IMPLICIT NONE
    
    WRITE ( unit=output_unit, fmt='(a)' ) 'Hard sphere potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'   
    WRITE ( unit=output_unit, fmt='(a)' ) 'Energy, kT = 1'   

  END SUBROUTINE introduction
  
  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays
    IMPLICIT NONE
    ALLOCATE ( r(3,n), v(3,n), coltime(n), partner(n) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE
    DEALLOCATE ( r, v, coltime, partner )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE update ( i, box, j_range )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: i       ! Index of atom of interest
    REAL,    INTENT(in) :: box     ! Simulation box length
    INTEGER, intent(in) :: j_range ! Range of j to be considered

    ! If j_range == gt, this routine loops over j > i seeking collisions at times shorter than coltime(i)
    ! Note that coltime(i) is set to a large value at the start, in this case
    ! If j_range == lt, the loop is over j < i, and the comparison is with coltime(j) in each case
    ! We use this approach so as to store information about each collision once only
    ! using the lower of the two indices

    INTEGER            :: j, j1, j2, k
    REAL, DIMENSION(3) :: rij, vij
    REAL               :: rijsq, vijsq, bij, tij, discr

    SELECT CASE ( j_range )
    CASE ( lt ) ! j < i
       j1 = 1
       j2 = i-1
    CASE ( gt ) ! j > i
       j1 = i+1
       j2 = n
       coltime(i) = HUGE(1.0)
    CASE default ! should never happen
       WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error ', j_range
       STOP 'Impossible error in update'
    END SELECT

    DO j = j1, j2 ! Loop over specified range of partners

       rij(:) = r(:,i) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) )
       rij(:) = rij(:) * box ! Now in sigma=1 units
       vij(:) = v(:,i) - v(:,j)
       bij    = DOT_PRODUCT ( rij, vij )

       IF ( bij < 0.0 ) THEN ! Test if collision is possible

          rijsq = SUM ( rij**2 )
          vijsq = SUM ( vij**2 )
          discr = bij ** 2 - vijsq * ( rijsq - 1.0 ) ! sigma**2 = 1.0

          IF ( discr > 0.0 ) THEN ! Test if collision is happening

             tij = ( -bij - SQRT ( discr ) ) / vijsq

             k = MIN(i,j)
             IF ( tij < coltime(k) ) THEN ! Test if collision needs storing

                coltime(k) = tij
                partner(k) = MAX(i,j)

             END IF ! End test if collision needs storing

          END IF ! End test if collision is happening

       END IF ! End test if collision is possible

    END DO ! End loop over specified range of partners

  END SUBROUTINE update

  FUNCTION overlap ( box )
    IMPLICIT NONE
    LOGICAL          :: overlap ! Returns flag indicating any pair overlaps
    REAL, INTENT(in) :: box     ! Simulation box length

    INTEGER            :: i, j
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq, box_sq, rij_mag

    overlap = .FALSE.
    box_sq  = box**2

    DO i = 1, n - 1
       DO j = i + 1, n

          rij(:) = r(:,i) - r(:,j)
          rij(:) = rij(:) - ANINT ( rij(:) )
          rij_sq = SUM ( rij**2 )  ! squared distance
          rij_sq = rij_sq * box_sq ! now in sigma=1 units

          IF ( rij_sq < 1.0 ) THEN
             rij_mag = SQRT(rij_sq)
             WRITE ( unit=error_unit, fmt='(a,2i5,f15.8)' ) 'Warning: i,j,rij = ', i, j, rij_mag
             overlap = .TRUE.
          END IF

       END DO
    END DO

  END FUNCTION overlap

  SUBROUTINE collide ( i, j, box, virial )
    IMPLICIT NONE
    INTEGER, INTENT(in)  :: i, j   ! Colliding atom indices
    REAL,    INTENT(in)  :: box    ! Simulation box length
    REAL,    INTENT(out) :: virial ! Collision contribution to pressure

    ! This routine implements collision dynamics, updating the velocities
    ! The colliding pair (i,j) is assumed to be in contact already
    
    REAL, DIMENSION(3) :: rij, vij
    REAL               :: factor
    
    rij(:) = r(:,i) - r(:,j)
    rij(:) = rij(:) - ANINT ( rij(:) ) ! Separation vector
    rij(:) = rij(:) * box              ! Now in sigma=1 units
    vij(:) = v(:,i) - v(:,j)           ! Relative velocity

    factor = DOT_PRODUCT ( rij, vij )
    vij    = -factor * rij

    v(:,i) = v(:,i) + vij
    v(:,j) = v(:,j) - vij
    virial = DOT_PRODUCT ( vij, rij ) / 3.0

  END SUBROUTINE collide

END MODULE md_module
