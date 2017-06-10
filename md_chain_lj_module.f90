! md_chain_lj_module.f90
! Force & constraint routines for MD, LJ chain
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
  PUBLIC :: zero_cm, force, spring
  PUBLIC :: rattle_a, rattle_b, worst_bond, milcshake_a, milcshake_b 

  ! Public data
  INTEGER,                              PUBLIC :: n     ! Number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r     ! Positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r_old ! Old positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v     ! Velocities (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: f     ! Non-bonded forces (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: g     ! Spring forces (3,n)

  ! Private data
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r_new                      ! New positions (3,n)
  LOGICAL, DIMENSION(:),   ALLOCATABLE :: move, moved                ! Iterative flags (n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: rij, rij_old, rij_new, vij ! Bond vectors etc (3,n-1)
  REAL,    DIMENSION(:),   ALLOCATABLE :: rijsq, lambda, rhs, rhsold ! Squared bond lengths, constraint equation terms (n-1)
  REAL,    DIMENSION(:),   ALLOCATABLE :: dd, dd_tmp                 ! Diagonal part of constraint matrix (n-1)
  REAL,    DIMENSION(:),   ALLOCATABLE :: dl, dl_tmp                 ! Sub-diagonal part of constraint matrix (n-2)
  REAL,    DIMENSION(:),   ALLOCATABLE :: du, du_tmp                 ! Super-diagonal part of constraint matrix (n-2)

  ! Public derived type
  TYPE, PUBLIC :: potential_type ! A composite variable for interactions comprising
     REAL    :: pot ! the potential energy
     LOGICAL :: ovr ! a flag indicating overlap (i.e. pot too high to use)
   CONTAINS
     PROCEDURE :: add_potential_type
     GENERIC   :: OPERATOR(+) => add_potential_type
  END TYPE potential_type

CONTAINS

  FUNCTION add_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the sum of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot  +   b%pot
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones chain, no cutoff, no shift, no periodic box'   
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'    
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'   
    WRITE ( unit=output_unit, fmt='(a)' ) 'All atomic masses the same m = 1'   
    WRITE ( unit=output_unit, fmt='(a)' ) 'All bond lengths the same'   

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays
    IMPLICIT NONE
    ALLOCATE ( r(3,n), v(3,n), f(3,n), g(3,n), r_old(3,n), r_new(3,n) )
    ALLOCATE ( move(n), moved(n) )
    ALLOCATE ( rij(3,n-1), rij_old(3,n-1), rij_new(3,n-1), vij(3,n-1) )
    ALLOCATE ( rijsq(n-1), lambda(n-1), rhs(n-1), rhsold(n-1) )
    ALLOCATE ( dd(n-1), dd_tmp(n-1) )
    ALLOCATE ( dl(n-2), dl_tmp(n-2), du(n-2), du_tmp(n-2) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE
    DEALLOCATE ( r, v, f, g, r_old, r_new )
    DEALLOCATE ( move, moved )
    DEALLOCATE ( rij, rij_old, rij_new, vij )
    DEALLOCATE ( rijsq, lambda, rhs, rhsold )
    DEALLOCATE ( dd, dd_tmp )
    DEALLOCATE ( dl, dl_tmp, du, du_tmp )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE zero_cm
    IMPLICIT NONE

    ! Routine to set centre-of-mass at the origin and zero the total momentum

    REAL, DIMENSION(3) :: c

    c = SUM ( r, dim=2 ) / REAL(n)
    r = r - SPREAD ( c, dim = 2, ncopies = n )

    c = SUM ( v, dim=2 ) / REAL(n)
    v = v - SPREAD ( c, dim = 2, ncopies = n )

  END SUBROUTINE zero_cm

  FUNCTION worst_bond ( bond ) RESULT ( diff_max )
    IMPLICIT NONE
    REAL             :: diff_max ! Returns max amount by which constraint is violated
    REAL, INTENT(in) :: bond     ! Bond length (the same throughout the chain)

    INTEGER            :: i
    REAL               :: diff, rij_mag
    REAL, DIMENSION(3) :: rij

    diff_max = 0.0
    DO i = 1, n-1
       rij      = r(:,i) - r(:,i+1)      ! Current bond vector
       rij_mag  = SQRT( SUM(rij**2) )    ! Current bond length
       diff     = ABS ( rij_mag - bond ) ! Absolute amount by which constraint is violated
       diff_max = MAX ( diff_max, diff ) ! Find maximum
    END DO
  END FUNCTION worst_bond

  SUBROUTINE force ( total )
    IMPLICIT NONE
    TYPE(potential_type), INTENT(out) :: total ! Composite of pot, ovr

    ! total%pot is the LJ potential energy for the atomic chain (omitting bonded neighbours)
    ! total%ovr is a warning flag that there is an overlap
    ! The Lennard-Jones energy and sigma parameters are taken to be epsilon = 1, sigma = 1
    ! Positions are assumed to be in these units
    ! Forces are calculated in the same units and stored in the array f
    ! NO box, NO periodic boundaries

    INTEGER              :: i, j
    REAL                 :: rij_sq, sr2, sr6, sr12, pair_vir
    REAL, DIMENSION(3)   :: rij, fij
    REAL, PARAMETER      :: sr2_ovr = 1.77 ! Overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    ! Initialize
    f     = 0.0
    total = potential_type (  pot=0.0, ovr=.FALSE. )

    DO i = 1, n - 2 ! Begin outer loop over atoms, stopping 2 short of the end

       DO j = i + 2, n ! Begin inner loop over atoms omitting nearest neighbour

          rij(:) = r(:,i) - r(:,j) ! Separation vector
          rij_sq = SUM ( rij**2 )  ! Squared separation

          sr2      = 1.0 / rij_sq  ! (sigma/rij)**2
          pair%ovr = sr2 > sr2_ovr ! Overlap if too close

          sr6      = sr2 ** 3
          sr12     = sr6 ** 2
          pair%pot = sr12 - sr6           ! LJ pair potential
          pair_vir = pair%pot + sr12      ! LJ pair virial
          fij      = rij * pair_vir * sr2 ! LJ Pair forces

          total  = total  + pair
          f(:,i) = f(:,i) + fij
          f(:,j) = f(:,j) - fij

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    ! Multiply results by numerical factors
    f         = f         * 24.0 ! 24*epsilon
    total%pot = total%pot * 4.0  ! 4*epsilon

  END SUBROUTINE force

  SUBROUTINE spring ( k_spring, bond, total_spr )
    IMPLICIT NONE
    REAL, INTENT(in)  :: k_spring  ! Spring force constant
    REAL, INTENT(in)  :: bond      ! Spring bond length
    REAL, INTENT(out) :: total_spr ! Total harmonic spring potential energy

    ! Calculates bond spring potential for atomic chain
    ! Forces are also calculated and stored in the array g
    ! NO box, NO periodic boundaries

    INTEGER            :: i, j
    REAL               :: rij_sq, rij_mag, pair_pot
    REAL, DIMENSION(3) :: rij, gij

    ! Initialize
    g         = 0.0
    total_spr = 0.0

    DO i = 1, n - 1 ! Begin loop over bonds
       j = i+1 ! Nearest neighbour

       rij(:)  = r(:,i) - r(:,j) ! Separation vector
       rij_sq  = SUM ( rij**2 )  ! Squared separation
       rij_mag = SQRT ( rij_sq ) ! Separation

       pair_pot = (rij_mag-bond) ** 2                ! Spring pair potential without numerical factor
       gij      = rij * ( bond - rij_mag ) / rij_mag ! Spring pair force without numerical factor

       total_spr = total_spr + pair_pot

       g(:,i) = g(:,i) + gij
       g(:,j) = g(:,j) - gij

    END DO ! End loop over bonds

    ! Multiply results by numerical factors
    total_spr = total_spr * 0.5 * k_spring
    g         = g * k_spring

  END SUBROUTINE spring

  SUBROUTINE rattle_a ( dt, bond )
    IMPLICIT NONE
    REAL, INTENT(in) :: dt   ! Time step
    REAL, INTENT(in) :: bond ! Bond length (the same throughout the chain)

    ! This subroutine iteratively adjusts the positions stored in the array r
    ! and the velocities stored in the array v, to satisfy the bond constraints

    ! On entry to this routine we assume:
    ! r_old stores the positions at the start of the step
    ! r stores the positions following the unconstrained drift and
    ! v stores the velocities following the first unconstrained half-kick
    ! On return from this routine, r and v will hold the constrained values

    REAL, DIMENSION(3) :: rij, rij_old, dr
    LOGICAL            :: done
    REAL               :: diffsq, dot, g
    INTEGER            :: i, j, iter

    REAL,    PARAMETER :: tol = 1.0e-9,  tol2 = 2.0 * tol, dot_tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 500

    iter     = 0
    done     = .FALSE.
    moved(:) = .TRUE. ! Ensures that we look at each bond at least once

    DO ! Iterative loop until done

       IF ( done ) EXIT ! done is equivalent to .not.any ( moved )

       done    = .TRUE.
       move(:) = .FALSE.

       DO i = 1, n-1 ! Loop over each constraint in turn
          j = i + 1 ! Partner atom in this constraint

          IF ( moved(i) .OR. moved(j) ) THEN ! Test whether need to re-examine ij

             rij = r(:,i) - r(:,j)             ! Current bond vector
             diffsq = bond**2 - SUM ( rij**2 ) ! Amount by which constraint is violated

             IF ( ABS ( diffsq ) > tol2*bond**2 ) THEN ! Test whether constraint not already satisfied

                rij_old = r_old(:,i) - r_old(:,j) ! Old vector determines direction of constraint force

                dot = DOT_PRODUCT ( rij_old, rij ) ! This should be of the order of bond**2

                IF ( dot < dot_tol * bond**2 ) THEN
                   WRITE ( unit=error_unit, fmt='(a,3f15.6)' ) 'Constraint failure', dot, dot_tol, bond**2
                   STOP 'Error in rattle_a'
                END IF

                ! In the following formulae, inverse masses are all unity
                g       = diffsq / ( 4.0 * dot ) 
                dr      = rij_old * g      ! Position adjustment
                r(:,i)  = r(:,i) + dr      ! Adjust i position
                r(:,j)  = r(:,j) - dr      ! Adjust j position
                v(:,i)  = v(:,i) + dr / dt ! Adjust i velocity
                v(:,j)  = v(:,j) - dr / dt ! Adjust j velocity
                move(i) = .TRUE.           ! Flag that we moved i
                move(j) = .TRUE.           ! Flag that we moved j
                done    = .FALSE.          ! Flag that we moved something

             END IF ! End test whether constraint not already satisfied

          END IF ! End test whether need to re-examine ij

       END DO ! End loop over each constraint in turn

       ! Prepare for next iteration
       moved(:) = move(:)
       iter     = iter + 1
       IF ( iter > iter_max ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations', iter, iter_max
          STOP 'Error in rattle_a'
       END IF

    END DO ! End iterative loop until done

  END SUBROUTINE rattle_a

  SUBROUTINE rattle_b ( dt, bond, wc ) ! Second part of velocity Verlet with constraints
    IMPLICIT NONE
    REAL, INTENT(in)  :: dt   ! Time step
    REAL, INTENT(in)  :: bond ! Bond length (the same throughout the chain)
    REAL, INTENT(out) :: wc   ! Constraint contribution to virial

    ! This subroutine iteratively adjusts the velocities stored in the array v
    ! to satisfy the time derivatives of the bond constraints
    ! Also returns constraint contribution to virial

    ! On entry to this routine we assume:
    ! r stores the positions at the end of the step with constraints applied
    ! v stores the velocities following the second unconstrained half-kick
    ! On return from this routine, v will hold the constrained values

    REAL, DIMENSION(3) :: rij, vij, dv
    LOGICAL            :: done
    REAL               :: dot, g
    INTEGER            :: i, j, iter
    REAL,    PARAMETER :: tol = 1.0e-9,  tol2 = 2.0 * tol, dot_tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 500

    iter     = 0
    done     = .FALSE.
    moved(:) = .TRUE.
    wc       = 0.0

    DO ! Iterative loop until done

       IF ( done ) EXIT

       done    = .TRUE.
       move(:) = .FALSE.

       DO i = 1, n-1 ! Loop over constraints
          j = i + 1 ! Partner atom for this constraint

          IF ( moved(i) .OR. moved(j) ) THEN ! Test whether need to re-examine ij
             vij = v(:,i) - v(:,j)
             rij = r(:,i) - r(:,j)
             dot = DOT_PRODUCT ( rij, vij )

             ! In the following formulae, inverse masses are all unity
             g  = -dot / ( 2.0 * bond**2 )

             IF ( ABS ( g ) > tol ) THEN ! Test whether constraint already satisfied

                wc      = wc + g * bond**2 ! Contribution to virial
                dv      = rij * g          ! Velocity adjustment
                v(:,i)  = v(:,i) + dv      ! Adjust velocity i
                v(:,j)  = v(:,j) - dv      ! Adjust velocity j
                move(i) = .TRUE.           ! Flag that we moved i
                move(j) = .TRUE.           ! Flag that we moved j
                done    = .FALSE.          ! Flag that we moved something

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

    wc = wc / (0.5*dt) / 3.0 ! Scale factors for virial

  END SUBROUTINE rattle_b

  SUBROUTINE milcshake_a ( dt, bond ) ! First part of velocity Verlet algorithm with constraints
    IMPLICIT NONE
    REAL, INTENT(in) :: dt   ! Time step
    REAL, INTENT(in) :: bond ! Bond length (the same throughout the chain)

    ! This subroutine iteratively adjusts the positions stored in the array r
    ! and the velocities stored in the array v, to satisfy the bond constraints
    ! using a tri-diagonal solver: here we use dgtsv from LAPACK.
    ! See AG Bailey, CP Lowe, and AP Sutton, J Comput Phys, 227, 8949 (2008)
    ! and AG Bailey, CP Lowe, and AP Sutton, Comput Phys Commun, 180, 594 (2009)

    ! On entry to this routine we assume:
    ! r_old stores the positions at the start of the step
    ! r stores the positions following the unconstrained drift and
    ! v stores the velocities following the first unconstrained half-kick
    ! On return from this routine, r and v will hold the constrained values

    INTEGER            :: iter, k
    REAL               :: max_error
    LOGICAL            :: info
    REAL,    PARAMETER :: tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 500

    k = n - 1 ! Number of constraints

    r_new = r ! Saves unconstrained positions

    ! Old and new (unconstrained) bond vectors
    rij_old(:,1:k) = r_old(:,1:k) - r_old(:,2:n)
    rij_new(:,1:k) = r_new(:,1:k) - r_new(:,2:n)

    ! Elements of tridiagonal matrix (dot products of old and new vectors)
    ! In this example, all masses are equal to unity
    dl(1:k-1) = -2.0*SUM ( rij_old(:,1:k-1)*rij_new(:,2:k),   dim=1 )       ! k-1 elements of lower-diagonal
    dd(1:k)   =  2.0*SUM ( rij_old(:,1:k)  *rij_new(:,1:k),   dim=1 ) / 0.5 ! k elements of diagonal
    du(1:k-1) = -2.0*SUM ( rij_old(:,2:k)  *rij_new(:,1:k-1), dim=1 )       ! k-1 elements of upper-diagonal

    ! Set up rhs of constraint equation
    rijsq(:)  = SUM ( rij_new(:,:)**2, dim=1 )
    rhs(:)    = bond**2 - rijsq(:)
    rhsold(:) = rhs(:)

    iter = 0

    DO ! Iterative loop until done

       ! Test for done
       max_error = MAXVAL ( ABS ( rijsq - bond**2 ) ) / ( 2.0*bond**2 )
       IF ( max_error <= tol ) EXIT 

       ! Copy tridiagonal elements (may be over-written by solver)
       dl_tmp(:) = dl(:)
       dd_tmp(:) = dd(:)
       du_tmp(:) = du(:)    
       CALL dgtsv( n-1, 1, dl_tmp, dd_tmp, du_tmp, rhs, n-1, info )              
       lambda(:) = rhs(:) ! Store result

       ! Constraint effects on position from lambda multipliers
       r = r_new
       r(:,1:k) = r(:,1:k) + SPREAD(lambda(1:k),dim=1,ncopies=3)*rij_old(:,1:k)
       r(:,2:n) = r(:,2:n) - SPREAD(lambda(1:k),dim=1,ncopies=3)*rij_old(:,1:k)

       ! New bond vectors
       rij(:,1:k) = r(:,1:k) - r(:,2:n)

       ! Prepare for next iteration
       rijsq(:)  = SUM ( rij(:,:)**2, dim=1 ) 
       rhs(:)    = bond**2 - rijsq(:) + rhsold(:)  
       rhsold(:) = rhs(:)               

       iter = iter + 1
       IF ( iter > iter_max ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations', iter, iter_max
          STOP 'Error in milcshake_a'
       END IF

    END DO ! End iterative loop until done

    ! Effect of constraints on velocities
    v(:,1:k) = v(:,1:k) + SPREAD(lambda(1:k),dim=1,ncopies=3)*rij_old(:,1:k)/dt
    v(:,2:n) = v(:,2:n) - SPREAD(lambda(1:k),dim=1,ncopies=3)*rij_old(:,1:k)/dt

  END SUBROUTINE milcshake_a

  SUBROUTINE milcshake_b ( dt, bond, wc ) ! Second part of velocity Verlet algorithm with constraints
    IMPLICIT NONE
    REAL, INTENT(in) :: dt   ! Time step
    REAL, INTENT(in) :: bond ! Bond length (the same throughout the chain)
    REAL, INTENT(out) :: wc  ! Constraint contribution to virial

    ! This subroutine adjusts the velocities stored in the array v
    ! to satisfy the time derivatives of the bond constraints
    ! using a tri-diagonal solver: here we use dgtsv from LAPACK.
    ! See AG Bailey, CP Lowe, and AP Sutton, J Comput Phys, 227, 8949 (2008)
    ! and AG Bailey, CP Lowe, and AP Sutton, Comput Phys Commun, 180, 594 (2009)
    ! Also returns constraint contribution to virial

    ! On entry to this routine we assume:
    ! r stores the positions at the end of the step with constraints applied
    ! v stores the velocities following the second unconstrained half-kick
    ! On return from this routine, v will hold the constrained values

    INTEGER :: k
    LOGICAL :: info

    k = n - 1 ! Number of constraints

    ! Relative velocities and bond vectors
    vij(:,1:k) = v(:,1:k) - v(:,2:n)
    rij(:,1:k) = r(:,1:k) - r(:,2:n)
    rhs(:) = - SUM ( vij*rij, dim=1 )

    ! Elements of tridiagonal matrix (dot products of old and new vectors)
    ! In this example, all masses are equal to unity
    dl(1:k-1) = -SUM ( rij(:,1:k-1)*rij(:,2:k),   dim=1 )       ! k-1 elements of lower-diagonal
    dd(1:k)   =  SUM ( rij(:,1:k)  *rij(:,1:k),   dim=1 ) / 0.5 ! k elements of diagonal
    du(1:k-1) = -SUM ( rij(:,2:k)  *rij(:,1:k-1), dim=1 )       ! k-1 elements of upper-diagonal

    CALL dgtsv ( n-1, 1, dl, dd, du, rhs, n-1, info ) 
    lambda(:) = rhs(:) ! Store result

    ! Effect of constraints on velocities
    v(:,1:k) = v(:,1:k) + SPREAD(lambda(1:k),dim=1,ncopies=3)*rij(:,1:k)
    v(:,2:n) = v(:,2:n) - SPREAD(lambda(1:k),dim=1,ncopies=3)*rij(:,1:k)

    wc = SUM(lambda) * bond**2
    wc = wc / (0.5*dt) / 3.0 ! scale factors for virial

  END SUBROUTINE milcshake_b

END MODULE md_module
