! test_pot_atom.f90
! Test potential, forces for atoms
PROGRAM test_pot_atom

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit

  USE test_pot_module, ONLY : n, force
  USE maths_module,  ONLY : init_random_seed, cross_product

  IMPLICIT NONE

  LOGICAL              :: ok
  INTEGER              :: i, xyz, try
  REAL, DIMENSION(3,n) :: r, f ! positions, forces
  REAL, DIMENSION(3)   :: axis, rsave, ftot, ttot
  REAL                 :: pot, potp, potm, fnum

  CHARACTER(len=2), DIMENSION(3), PARAMETER :: cf = ['Fx','Fy','Fz']

  ! Any of the following parameters could be empirically adjusted
  REAL,    PARAMETER :: delta = 0.00001, d_min = 0.3, d_max = 1.5, pot_max = 10.0
  INTEGER, PARAMETER :: ntry = 1000

  ! Initialize random number generator                          
  CALL init_random_seed

  ! Make a number of attempts to place the atoms
  ok = .FALSE.
  DO try = 1, ntry 
     CALL random_positions ( d_min, d_max, r )
     CALL force ( r, pot )
     IF ( ABS(pot) < pot_max ) THEN
        ok = .TRUE.
        EXIT
     END IF
  END DO
  IF (.NOT.ok) WRITE ( output_unit, fmt='(a)' ) 'Warning: atom placement failed'

  CALL force ( r, pot, f ) ! Calculation of potential, and analytical forces
  WRITE ( output_unit, fmt='(a,t25,f15.5)' ) 'Potential energy = ', pot

  ! Check momentum and angular momentum conservation
  ftot = SUM ( f, dim=2 )
  ttot = 0.0
  DO i = 1, n
     ttot = ttot + cross_product ( r(:,i), f(:,i) )
  END DO
  WRITE ( output_unit, fmt='(a,t25,3es15.3)') 'Total force      = ', ftot
  WRITE ( output_unit, fmt='(a,t25,3es15.3)') 'Total torque     = ', ttot

  WRITE ( output_unit, fmt='(2a10,3a15)') 'Atom', 'Component', 'Exact', 'Numerical', 'Difference'

  DO i = 1, n

     DO xyz = 1, 3 ! Loop to calculate numerical forces
        axis      = 0.0
        axis(xyz) = 1.0 ! Pick axis
        rsave     = r(:,i) ! Save position
        r(:,i)    = rsave + delta*axis ! translate
        CALL force ( r, potp )
        r(:,i)    = rsave - delta*axis ! translate
        CALL force ( r, potm )
        r(:,i)    = rsave ! Restore position
        fnum      = -(potp-potm)/(2.0*delta)
        WRITE ( output_unit, fmt='(i10,a10,2f15.5,es15.3)' ) i, cf(xyz), f(xyz,i), fnum, f(xyz,i)-fnum
     END DO ! End loop to calculate numerical forces

  END DO ! End loop over molecules

CONTAINS

  SUBROUTINE random_positions ( d_min, d_max, r )
    USE maths_module, ONLY : random_orientation_vector
    REAL,                  INTENT(in)  :: d_min, d_max
    REAL, DIMENSION (3,n), INTENT(out) :: r

    INTEGER :: i, j
    REAL    :: d, zeta
    LOGICAL :: ok

    ! This routine places n atoms at positions which are all within
    ! the desired range d_min ... d_max of each other.
    ! Obviously this approach can fail, depending on the value of n
    ! and the supplied distances.
    r(:,1) = [0.0,0.0,0.0] ! first one at origin

    DO i = 2, n ! Loop over remaining atoms

       DO ! loop until successful
          r(:,i) = random_orientation_vector ( ) ! direction of r
          CALL RANDOM_NUMBER ( zeta )
          d      = d_min + (d_max-d_min)*zeta ! magnitude of r
          r(:,i) = r(:,i) * d                 ! within desired range of origin
          ok = .TRUE.
          DO j = 2, i-1 ! check intermediate atoms if any
             d  = SQRT ( SUM((r(:,i)-r(:,j))**2) )
             ok = ok .AND. ( d >= d_min ) .AND. ( d <= d_max )
          END DO ! end check intermediate atoms if any
          IF ( ok ) EXIT
       END DO ! End loop until successful

    END DO ! End loop over remaining atoms

  END SUBROUTINE random_positions

END PROGRAM test_pot_atom
