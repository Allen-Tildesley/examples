! md_chain_lj_module.f90
! Force & constraint routines for MD, LJ chain
MODULE md_chain_lj_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, r, v, f, f_spring
  PUBLIC :: allocate_arrays, deallocate_arrays, force, spring
  PUBLIC :: rattle_a, rattle_b, worst_bond, milcshake_a, milcshake_b 

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: v ! velocities (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f ! forces (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f_spring                   ! spring forces (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r_old, r_new               ! old, new positions (3,n)
  LOGICAL, DIMENSION(:),   ALLOCATABLE :: moving, moved              ! iterative flags (n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: rij, rij_old, rij_new, vij ! bond vectors (3,n-1)
  REAL,    DIMENSION(:),   ALLOCATABLE :: rijsq, lambda, rhs, rhsold ! squared bond lengths, constraint equation terms (n-1)
  REAL,    DIMENSION(:),   ALLOCATABLE :: dd, dd_tmp                 ! diagonal part of constraint matrix (n-1)
  REAL,    DIMENSION(:),   ALLOCATABLE :: dl, dl_tmp                 ! sub-diagonal part of constraint matrix (n-2)
  REAL,    DIMENSION(:),   ALLOCATABLE :: du, du_tmp                 ! super-diagonal part of constraint matrix (n-2)

  ! all masses are unity
  ! all bond lengths are the same
  ! no periodic boundary conditions

CONTAINS

  SUBROUTINE allocate_arrays
    ALLOCATE ( r(3,n), v(3,n), f(3,n), f_spring(3,n), r_old(3,n), r_new(3,n) )
    ALLOCATE ( moving(n), moved(n) )
    ALLOCATE ( rij(3,n-1), rij_old(3,n-1), rij_new(3,n-1), vij(3,n-1) )
    ALLOCATE ( rijsq(n-1), lambda(n-1), rhs(n-1), rhsold(n-1) )
    ALLOCATE ( dd(n-1), dd_tmp(n-1) )
    ALLOCATE ( dl(n-2), dl_tmp(n-2), du(n-2), du_tmp(n-2) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, v, f, f_spring, r_old, r_new )
    DEALLOCATE ( moving, moved )
    DEALLOCATE ( rij, rij_old, rij_new, vij )
    DEALLOCATE ( rijsq, lambda, rhs, rhsold )
    DEALLOCATE ( dd, dd_tmp )
    DEALLOCATE ( dl, dl_tmp, du, du_tmp )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( pot )
    REAL, INTENT(out) :: pot          ! total potential energy

    ! Calculates shifted WCA LJ potential for atomic chain (omitting bonded neighbours)
    ! The Lennard-Jones energy and sigma parameters are taken to be epsilon = 1, sigma = 1
    ! Forces are calculated in the same units
    ! NO periodic boundaries

    INTEGER            :: i, j
    REAL               :: rij_sq, sr2, sr6, sr12, potij
    REAL, PARAMETER    :: r_cut_sq = 2.0**(1.0/3.0) ! minimum of LJ potential
    REAL, DIMENSION(3) :: rij, fij

    f     = 0.0
    pot   = 0.0

    DO i = 1, n - 2 ! Begin outer loop over atoms, stopping 2 short of the end

       DO j = i + 2, n ! Begin inner loop over atoms omitting nearest neighbour

          rij(:) = r(:,i) - r(:,j)           ! separation vector
          rij_sq = SUM ( rij**2 )            ! squared separation

          IF ( rij_sq < r_cut_sq ) THEN

             sr2    = 1.0 / rij_sq
             sr6    = sr2 ** 3
             sr12   = sr6 ** 2
             potij  = sr12 - sr6 + 0.25 ! shifted potential
             pot    = pot + potij
             fij    = rij * (2.0*sr12 - sr6) / rij_sq
             f(:,i) = f(:,i) + fij
             f(:,j) = f(:,j) - fij

          END IF

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

    ! Multiply results by numerical factors
    f      = f   * 24.0
    pot    = pot * 4.0

  END SUBROUTINE force

  SUBROUTINE spring ( k_spring, bond, pot_spring )
    REAL, INTENT(in)  :: k_spring   ! spring force constant
    REAL, INTENT(in)  :: bond       ! spring bond length
    REAL, INTENT(out) :: pot_spring ! total spring potential energy

    ! Calculates bond spring potential for atomic chain
    ! NO periodic boundaries

    INTEGER            :: i, j
    REAL               :: rij_sq, rij_mag, potij
    REAL, DIMENSION(3) :: rij, fij

    f_spring   = 0.0
    pot_spring = 0.0

    DO i = 1, n - 1 ! Begin outer loop over bonds
       j = i+1

       rij(:) = r(:,i) - r(:,j)  ! separation vector
       rij_sq = SUM ( rij**2 )   ! squared separation
       rij_mag = SQRT ( rij_sq ) ! separation

       potij         = 0.5*k_spring*(rij_mag-bond)**2 ! spring potential
       pot_spring    = pot_spring + potij
       fij           = rij * k_spring*(bond-rij_mag)/rij_mag
       f_spring(:,i) = f_spring(:,i) + fij
       f_spring(:,j) = f_spring(:,j) - fij

    END DO ! End outer loop over bonds

  END SUBROUTINE spring

  SUBROUTINE rattle_a ( dt, bond ) ! first part of velocity Verlet algorithm with constraints
    IMPLICIT NONE
    REAL, INTENT(in) :: dt, bond

    REAL, DIMENSION(3) :: rij, rij_old, dr
    LOGICAL            :: done
    REAL               :: diffsq, dot, g
    INTEGER            :: i, j, iter

    REAL,    PARAMETER :: tol = 1.0e-9,  tol2 = 2.0 * tol, dot_tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 500

    r_old = r            ! save, since constraints act along the old vectors
    v = v + 0.5 * dt * f ! Kick half-step
    r = r + dt * v       ! Drift step

    iter      = 0
    done      = .FALSE.
    moving(:) = .FALSE.
    moved(:)  = .TRUE.

    DO ! iterative loop

       IF ( done ) EXIT ! done is equivalent to .not.any ( moved )

       done = .TRUE.

       DO i = 1, n-1 ! loop over each constraint in turn
          j = i + 1

          IF ( moved(i) .OR. moved(j) ) THEN ! test whether need to re-examine ij

             rij = r(:,i) - r(:,j)             ! current bond vector
             diffsq = bond**2 - SUM ( rij**2 ) ! amount by which constraint is violated

             IF ( ABS ( diffsq ) > tol2*bond**2 ) THEN ! test whether constraint already satisfied

                rij_old = r_old(:,i) - r_old(:,j)
                dot = DOT_PRODUCT ( rij_old, rij )
                IF ( dot < dot_tol * bond**2 ) THEN
                   WRITE ( unit=error_unit, fmt='(a,3f15.5)' ) 'Constraint failure', dot, dot_tol, bond**2
                   STOP 'Error in rattle_a'
                END IF

                ! in the following formulae, inverse masses are all unity
                g      = diffsq / ( 4.0 * dot ) 
                dr     = rij_old * g      ! position adjustment
                r(:,i) = r(:,i) + dr      ! adjust i position
                r(:,j) = r(:,j) - dr      ! adjust j position
                v(:,i) = v(:,i) + dr / dt ! adjust i velocity
                v(:,j) = v(:,j) - dr / dt ! adjust j velocity
                moving(i)  = .TRUE.       ! flag that we moved i
                moving(j)  = .TRUE.       ! flag that we moved j
                done = .FALSE.            ! flag that we moved something

             END IF ! end test whether constraint already satisfied

          END IF ! end test whether need to re-examine ij

       END DO ! end loop over each constraint in turn

       ! prepare for next iteration
       moved(:)  = moving(:)
       moving(:) = .FALSE.
       iter = iter + 1
       IF ( iter > iter_max ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations', iter, iter_max
          STOP 'Error in rattle_a'
       END IF

    END DO ! end of iterative loop

  END SUBROUTINE rattle_a

  SUBROUTINE rattle_b ( dt, bond, wc ) ! second part of velocity Verlet with constraints
    IMPLICIT NONE
    REAL, INTENT(in)  :: dt, bond
    REAL, INTENT(out) :: wc ! constraint contribution to virial

    REAL, DIMENSION(3) :: rij, vij, dv
    LOGICAL            :: done
    REAL               :: dot, g
    INTEGER            :: i, j, iter
    REAL,    PARAMETER :: tol = 1.0e-9,  tol2 = 2.0 * tol, dot_tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 500

    v = v + 0.5 * dt * f ! Kick half-step

    iter      = 0
    done      = .FALSE.
    moving(:) = .FALSE.
    moved(:)  = .TRUE.
    wc        = 0.0

    DO ! iterative loop

       IF ( done ) EXIT

       done = .TRUE.

       DO i = 1, n-1 ! loop over constraints
          j = i + 1

          IF ( moved(i) .OR. moved(j) ) THEN ! test whether need to re-examine ij
             vij = v(:,i) - v(:,j)
             rij = r(:,i) - r(:,j)
             dot = DOT_PRODUCT ( rij, vij )

             ! in the following formulae, inverse masses are all unity
             g  = -dot / ( 2.0 * bond**2 )

             IF ( ABS ( g ) > tol ) THEN ! test whether constraint already satisfied

                wc = wc + g * bond**2 ! contribution to virial
                dv = rij * g          ! velocity adjustment
                v(:,i) = v(:,i) + dv  ! adjust velocity i
                v(:,j) = v(:,j) - dv  ! adjust velocity j
                moving(i) = .TRUE.    ! flag that we moved i
                moving(j) = .TRUE.    ! flag that we moved j
                done = .FALSE.        ! flag that we moved something

             END IF ! end test whether constraint already satisfied

          END IF ! end test whether need to re-examine ij

       END DO ! end loop over constraints

       moved(:)  = moving(:)
       moving(:) = .FALSE.

       iter = iter + 1
       IF ( iter > iter_max ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations', iter, iter_max
          STOP 'Error in rattle_b'
       END IF

    END DO ! end iterative loop

    wc = wc / (0.5*dt) / 3.0 ! scale factors for virial

  END SUBROUTINE rattle_b

  FUNCTION worst_bond ( bond ) RESULT ( diff_max )
    REAL, INTENT(in) :: bond
    REAL             :: diff_max

    INTEGER            :: i
    REAL               :: diff
    REAL, DIMENSION(3) :: rij

    diff_max = 0.0
    DO i = 1, n-1
       rij = r(:,i) - r(:,i+1)          ! current bond vector
       diff = SQRT( SUM(rij**2) ) -bond ! amount by which constraint is violated
       diff_max = MAX ( diff_max, diff )
    END DO
  END FUNCTION worst_bond

  SUBROUTINE milcshake_a ( dt, bond )
    REAL, INTENT(in) :: dt, bond

    INTEGER :: i, iter
    REAL    :: max_error
    LOGICAL :: done, info
    REAL,    PARAMETER :: tol = 1.0e-9
    INTEGER, PARAMETER :: iter_max = 500

    v = v + 0.5 * dt * f ! Kick half-step
    r_new = r + dt * v   ! Drift step

    ! Old and new (non-constrained) bond vectors
    rij_old(:,1:n-1) = r(:,1:n-1)-r(:,2:n)
    rij_new(:,1:n-1) = r_new(:,1:n-1)-r_new(:,2:n)

    ! Elements of tridiagonal matrix
    ! In this example, all masses are equal to unity
    dd(1) = DOT_PRODUCT ( rij_old(:,1), rij_new(:,1) ) / 0.5
    du(1) = DOT_PRODUCT ( rij_old(:,2), rij_new(:,1) )
    DO i = 2, n-2
       dl(i-1) = DOT_PRODUCT ( rij_old(:,i-1), rij_new(:,i) )
       dd(i)   = DOT_PRODUCT ( rij_old(:,i),   rij_new(:,i) ) / 0.5
       du(i)   = DOT_PRODUCT ( rij_old(:,i+1), rij_new(:,i) )
    ENDDO
    dl(n-2) = DOT_PRODUCT ( rij_old(:,n-3), rij_new(:,n-2) ) 
    dd(n-1) = DOT_PRODUCT ( rij_old(:,n-1), rij_new(:,n-1) ) / 0.5

    dl(:) = -2.0*dl(:)
    dd(:) =  2.0*dd(:)  
    du(:) = -2.0*du(:)

    ! Set up rhs of constraint equation
    rijsq(:)  = SUM ( rij_new(:,:)**2, dim=1 )
    max_error = MAXVAL ( ABS(rijsq-bond**2) )/( 2.0*bond**2 )
    rhs(:)    = bond**2 - rijsq(:)
    rhsold(:) = rhs(:)

    iter = 0
    done = .FALSE.

    DO ! iterative loop
       IF ( done ) EXIT 

       ! Reset tridiagonal elements (may have been over-written by solver)
       dl_tmp(:) = dl(:)
       dd_tmp(:) = dd(:)
       du_tmp(:) = du(:)    
       CALL dgtsv( n-1, 1, dl_tmp, dd_tmp, du_tmp, rhs, n-1, info )              
       lambda(:) = rhs(:)

       ! Constraint effects on position from lambda multipliers
       r(:,1) = r_new(:,1) + lambda(1)*rij_old(:,1)
       DO i = 2, n-1
          r(:,i) = r_new(:,i) + lambda(i)*rij_old(:,i) - lambda(i-1)*rij_old(:,i-1)
       ENDDO
       r(:,n) = r_new(:,n) - lambda(n-1)*rij_old(:,n-1)

       ! New bond vectors
       rij(:,1:n-1) = r(:,1:n-1)-r(:,2:n)

       ! Prepare for next iteration
       rijsq(:)  = SUM ( rij(:,:)**2, dim=1 ) 
       rhs(:)    = bond**2 - rijsq(:) + rhsold(:)  
       rhsold(:) = rhs(:)               

       max_error = MAXVAL ( ABS(rijsq-bond**2) )/( 2.0*bond**2 )
       done      = ( max_error <= tol )                                

       iter = iter + 1
       IF ( iter > iter_max ) then
          WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations', iter, iter_max
          STOP 'Error in milcshake_a'
       END IF

    END DO

    ! Effect of constraints on velocities
    v(:,1) = v(:,1) + lambda(1)*rij_old(:,1)/dt
    DO i = 2, n-1
       v(:,i) = v(:,i) + lambda(i)*rij_old(:,i)/dt - lambda(i-1)*rij_old(:,i-1)/dt
    ENDDO
    v(:,n) = v(:,n) - lambda(n-1)*rij_old(:,n-1)/dt

  END SUBROUTINE milcshake_a

  SUBROUTINE milcshake_b ( dt, bond, wc )
    REAL, INTENT(in)  :: dt, bond
    REAL, INTENT(out) :: wc

    INTEGER :: i
    LOGICAL :: info

    v = v + 0.5 * dt * f ! Kick half-step

    ! Relative velocities and bond vectors
    vij(:,1:n-1) = v(:,1:n-1) - v(:,2:n)
    rij(:,1:n-1) = r(:,1:n-1) - r(:,2:n)
    rhs(:) = - SUM ( vij*rij, dim=1 )

    dd(1) =  DOT_PRODUCT ( rij(:,1), rij(:,1) ) / 0.5
    du(1) = -DOT_PRODUCT ( rij(:,2), rij(:,1) )
    DO i = 2, n-2
       dl(i-1) = -DOT_PRODUCT ( rij(:,i-1), rij(:,i) )
       dd(i)   =  DOT_PRODUCT ( rij(:,i),   rij(:,i) ) / 0.5
       du(i)   = -DOT_PRODUCT ( rij(:,i+1), rij(:,i) )
    ENDDO
    dl(n-2) = -DOT_PRODUCT ( rij(:,n-2), rij(:,n-1) )
    dd(n-1) =  DOT_PRODUCT ( rij(:,n-1), rij(:,n-1) ) / 0.5

    CALL dgtsv ( n-1, 1, dl, dd, du, rhs, n-1, info ) 
    lambda(:) = rhs(:)

    v(:,1) = v(:,1) + lambda(1)*rij(:,1)
    DO i = 2, n-1
       v(:,i) = v(:,i) + ( lambda(i)*rij(:,i) - lambda(i-1)*rij(:,i-1) )
    ENDDO
    v(:,n) = v(:,n) - lambda(n-1)*rij(:,n-1)

    wc = SUM(lambda) * bond**2
    wc = wc / (0.5*dt) / 3.0 ! scale factors for virial

  END SUBROUTINE milcshake_b

END MODULE md_chain_lj_module
