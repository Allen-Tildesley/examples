! initialize_chain_module.f90
! Routines to initialize configurations and velocities for chain molecule
MODULE initialize_chain_module

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <dominic.tildesley@epfl.ch> ("the authors"),                          !
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
  PUBLIC :: allocate_arrays, deallocate_arrays
  PUBLIC :: initialize_random, initialize_velocities

  ! Public data
  INTEGER,                           PUBLIC :: n      ! Number of atoms
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r      ! Positions (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v      ! Velocities (3,n)

CONTAINS

  SUBROUTINE allocate_arrays
    IMPLICIT NONE

    ALLOCATE ( r(3,n), v(3,n) )

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, v )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE initialize_random ( bond, soft )
    USE maths_module, ONLY : random_vector
    IMPLICIT NONE
    REAL,    INTENT(in) :: bond ! Chain bond length
    LOGICAL, INTENT(in) :: soft ! Flag for soft interactions (no overlap check)

    ! Chooses chain positions randomly, at desired bond length, avoiding overlap

    REAL               :: diff_sq
    REAL, DIMENSION(3) :: r_cm
    INTEGER            :: i, iter

    REAL,    PARAMETER :: tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 1000

    WRITE ( unit=output_unit, fmt='(a)' ) 'Chain, bonds randomly oriented, avoiding overlaps'

    r(:,1) = [0.0,0.0,0.0]        ! First atom at origin
    r(:,2) = bond*random_vector() ! Second atom at random position (bond length away)

    DO i = 3, n ! Begin loop over atom indices

       iter = 0
       DO ! Loop until non-overlapping position found (N must not be too large!)
          r(:,i) = r(:,i-1) + bond*random_vector() ! Subsequent atoms randomly placed (bond length away)

          IF ( soft ) EXIT ! No overlap test
          IF ( .NOT. overlap ( i, 1, i-2 ) ) EXIT ! Check all so far except bonded neighbour

          iter = iter + 1
          IF ( iter > iter_max ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations ', iter, iter_max
             STOP 'Error in initialize_chain_random'
          END IF

       END DO ! End loop until non-overlapping position found

    END DO ! End loop over atom indices

    r_cm(:) = SUM ( r(:,:), dim=2 ) / REAL(n)               ! Compute centre of mass positions
    r(:,:)  = r(:,:) - SPREAD ( r_cm(:), dim=2, ncopies=n ) ! Shift centre of mass to the origin

    DO i = 1, n-1 ! Loop to confirm bond lengths
       diff_sq = SUM ( (r(:,i)-r(:,i+1))**2 ) - bond**2
       IF ( ABS(diff_sq) > tol ) WRITE ( unit=error_unit, fmt='(a,2i15,f15.8)' ) 'Bond length warning ', i, i+1, diff_sq
    END DO ! End loop to confirm bond lengths

  END SUBROUTINE initialize_random

  SUBROUTINE initialize_velocities ( temperature )
    USE maths_module, ONLY : random_normals, cross_product, solve, outer_product
    IMPLICIT NONE
    REAL, INTENT(in) :: temperature ! reduced temperature

    ! Chooses velocities from Maxwell-Boltzmann (Gaussian) distribution
    ! For simplicity, we just pick each atom velocity randomly and
    ! apply bond constraints afterwards
    ! In between, we take steps to remove linear and angular momentum
    ! since the configuration may be used in MD simulations without periodic boundaries
    ! in which case both these quantities are conserved
    ! NB there is at present no check for a singular inertia tensor in the angular momentum fix!
    ! We assume centre of mass is already at the origin
    ! We assume unit molecular mass and employ Lennard-Jones units
    ! property                  units
    ! energy                    epsilon ( = 1 )
    ! molecular mass            m ( = 1 )
    ! velocity v                sqrt(epsilon/m)

    REAL                 :: temp
    REAL, DIMENSION(3)   :: v_cm, r_cm, ang_mom, ang_vel
    REAL, DIMENSION(3,3) :: inertia
    INTEGER              :: i, xyz
    REAL, PARAMETER      :: tol = 1.e-6

    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Chain velocities at temperature', temperature

    ! Confirm centre-of-mass is at origin
    r_cm(:) = SUM ( r(:,:), dim=2 ) / REAL ( n ) ! Compute centre of mass
    IF ( ANY ( ABS ( r_cm ) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.8)' ) 'Centre of mass error', r_cm
       STOP 'Error in initialize_chain_velocities'
    END IF

    CALL random_normals ( 0.0, SQRT(temperature), v ) ! Choose 3N random velocities

    ! Compute and remove total momentum
    v_cm(:) = SUM ( v(:,:), dim=2 ) / REAL ( n )            ! Compute centre of mass velocity
    v(:,:)  = v(:,:) - SPREAD ( v_cm(:), dim=2, ncopies=n ) ! Set net momentum to zero

    ! Compute total angular momentum and moment of inertia tensor
    ang_mom = 0.0
    inertia = 0.0
    DO i = 1, n
       ang_mom = ang_mom + cross_product ( r(:,i), v(:,i) )
       inertia = inertia - outer_product ( r(:,i), r(:,i) )
       FORALL ( xyz=1:3 ) inertia(xyz,xyz) = inertia(xyz,xyz) + DOT_PRODUCT ( r(:,i), r(:,i) )
    END DO

    ! Solve linear system to get angular velocity
    ang_vel = solve ( inertia, ang_mom )

    ! Remove angular momentum
    DO i = 1, n
       v(:,i) = v(:,i) - cross_product ( ang_vel, r(:,i) )
    END DO

    ! Apply bond constraints (which should not introduce linear or angular momentum)
    CALL rattle_b

    ! Scale velocities to get correct temperature
    ! Number of degrees of freedom is 3*n - (n-1) bonds - 6 for angular and linear momentum
    temp   = SUM(v(:,:)**2) / REAL ( 3*n - (n-1) - 6 )
    v(:,:) = v(:,:) * SQRT ( temperature / temp )

    ! Final check on angular and linear momenta
    v_cm    = 0.0
    ang_mom = 0.0
    DO i = 1, n
       v_cm    = v_cm + v(:,i)
       ang_mom = ang_mom + cross_product ( r(:,i), v(:,i) )
    END DO

    IF ( ANY ( ABS ( v_cm ) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.8)' ) 'Linear momentum error', v_cm
       STOP 'Error in initialize_chain_velocities'
    END IF

    IF ( ANY ( ABS ( ang_mom ) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.8)' ) 'Angular momentum error', ang_mom
       STOP 'Error in initialize_chain_velocities'
    END IF

  END SUBROUTINE initialize_velocities

  SUBROUTINE rattle_b ! A version of velocity Verlet constraint algorithm
    IMPLICIT NONE

    ! This subroutine iteratively adjusts the velocities stored in the array v
    ! to satisfy the time derivatives of the bond constraints

    REAL, DIMENSION(3) :: rij, vij, dv
    LOGICAL            :: done
    REAL               :: g
    INTEGER            :: i, j, iter
    REAL,    PARAMETER :: tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 500
    LOGICAL, DIMENSION(:), ALLOCATABLE :: move, moved

    ALLOCATE ( move(n), moved(n) )

    iter     = 0
    done     = .FALSE.
    moved(:) = .TRUE.

    DO ! Iterative loop until done

       IF ( done ) EXIT

       done    = .TRUE.
       move(:) = .FALSE.

       DO i = 1, n-1 ! Loop over constraints
          j = i + 1 ! Partner atom for this constraint

          IF ( moved(i) .OR. moved(j) ) THEN ! Test whether need to re-examine ij
             vij = v(:,i) - v(:,j)
             rij = r(:,i) - r(:,j)

             ! In the following formulae, inverse masses are all unity
             g = -0.5*DOT_PRODUCT ( rij, vij ) / DOT_PRODUCT ( rij, rij )

             IF ( ABS ( g ) > tol ) THEN ! Test whether constraint already satisfied

                dv      = rij * g     ! Velocity adjustment
                v(:,i)  = v(:,i) + dv ! Adjust velocity i
                v(:,j)  = v(:,j) - dv ! Adjust velocity j
                move(i) = .TRUE.      ! Flag that we moved i
                move(j) = .TRUE.      ! Flag that we moved j
                done    = .FALSE.     ! Flag that we moved something

             END IF ! End test whether constraint already satisfied

          END IF ! End test whether need to re-examine ij

       END DO ! End loop over constraints

       ! Prepare for next iteration
       moved(:) = move(:)
       iter     = iter + 1
       IF ( iter > iter_max ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations', iter, iter_max
          STOP 'Error in rattle_b'
       END IF

    END DO ! End iterative loop until done

    DEALLOCATE ( move, moved )

  END SUBROUTINE rattle_b

  FUNCTION overlap ( i, j1, j2 )
    IMPLICIT NONE
    LOGICAL             :: overlap ! Returns a flag indicating overlap with any j<i
    INTEGER, INTENT(in) :: i       ! Index of atom/molecule to be checked
    INTEGER, INTENT(in) :: j1      ! First j index
    INTEGER, INTENT(in) :: j2      ! Last j index

    ! This routine checks for overlaps of atoms

    INTEGER            :: j
    REAL, DIMENSION(3) :: rij
    REAL               :: rij_sq

    overlap = .FALSE.
    IF ( j2 < j1 ) RETURN

    DO j = j1, j2
       rij    = r(:,i) - r(:,j)
       rij_sq = SUM ( rij**2 )

       IF ( rij_sq < 1.0 ) THEN
          overlap = .TRUE.
          RETURN
       END IF

    END DO

  END FUNCTION overlap

END MODULE initialize_chain_module
