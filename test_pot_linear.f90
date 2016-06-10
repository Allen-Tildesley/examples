! test_pot_linear.f90
! Test potential, forces, torques for linear molecule
PROGRAM test_pot
  USE test_pot_module, ONLY : n, force
  USE utility_module,  ONLY : init_random_seed, rotate_vector, cross_product
  IMPLICIT NONE

  LOGICAL              :: ok
  INTEGER              :: i, xyz, try
  REAL, DIMENSION(3,n) :: r, e, f, t ! positions, orientations, forces, torques
  REAL, DIMENSION(3)   :: axis, rsave, esave, ftot, ttot
  REAL                 :: pot, potp, potm, fnum, tnum

  CHARACTER(len=2), DIMENSION(3), PARAMETER :: cf = ['Fx','Fy','Fz'], ct = ['Tx','Ty','Tz']

  ! Any of the following parameters could be empirically adjusted
  REAL, PARAMETER      :: delta = 0.00001, d_min = 0.3, d_max = 1.5, pot_max = 10.0
  INTEGER, PARAMETER   :: ntry = 1000

  ! Initialize random number generator                          
  CALL init_random_seed()

  ! Make a number of attempts to place the atoms
  ok = .FALSE.
  DO try = 1, ntry 
     CALL random_orientations ( e )
     CALL random_positions ( d_min, d_max, r )
     CALL force ( r, e, pot )
     IF ( ABS(pot) < pot_max ) THEN
        ok = .TRUE.
        EXIT
     END IF
  END DO
  IF (.NOT.ok) WRITE(*,'(a)') 'Warning: molecule placement failed'

  ! Calculation of potential, and analytical forces and torques
  CALL force ( r, e, pot, f, t )
  WRITE(*,'(a,t25,f15.5)') 'Potential energy = ', pot

  ! Check momentum and angular momentum conservation
  ftot = SUM ( f, dim=2 )
  ttot = SUM ( t, dim=2 )
  DO i = 1, n
     ttot = ttot + cross_product ( r(:,i), f(:,i) )
  END DO
  WRITE(*,'(a,t25,3es15.3)') 'Total force      = ', ftot
  WRITE(*,'(a,t25,3es15.3)') 'Total torque     = ', ttot

  WRITE(*,'(2a10,3a15)') 'Molecule', 'Component', 'Exact', 'Numerical', 'Difference'

  DO i = 1, n

     DO xyz = 1, 3 ! Loop to calculate numerical forces
        axis      = 0.0
        axis(xyz) = 1.0 ! Pick axis
        rsave     = r(:,i) ! Save position
        r(:,i)    = rsave + delta*axis ! translate
        CALL force ( r, e, potp )
        r(:,i)    = rsave - delta*axis ! translate
        CALL force ( r, e, potm )
        r(:,i)    = rsave ! Restore position
        fnum      = -(potp-potm)/(2.0*delta)
        WRITE(*,'(i10,a10,2f15.5,es15.3)') i, cf(xyz), f(xyz,i), fnum, f(xyz,i)-fnum
     END DO ! End loop to calculate numerical forces

     DO xyz = 1, 3 ! Loop to calculate numerical torques
        axis      = 0.0
        axis(xyz) = 1.0 ! Pick axis
        esave     = e(:,i) ! Save orientation
        e(:,i)    = rotate_vector ( delta, axis, esave ) ! rotate
        CALL force ( r, e, potp )
        e(:,i)    = rotate_vector ( -delta, axis, esave ) ! rotate
        CALL force ( r, e, potm )
        e(:,i)    = esave ! Restore orientation
        tnum      = -(potp-potm)/(2.0*delta)
        WRITE(*,'(i10,a10,2f15.5,es15.3)') i, ct(xyz), t(xyz,i), tnum, t(xyz,i)-tnum
     END DO ! End loop to calculate numerical torques on A

  END DO ! End loop over molecules

CONTAINS

  SUBROUTINE random_orientations ( e )
    USE utility_module, ONLY : random_orientation_vector
    REAL, DIMENSION (3,n), INTENT(out) :: e

    INTEGER :: i

    DO i = 1, n
       CALL random_orientation_vector ( e(:,i) )
    END DO

  END SUBROUTINE random_orientations

  SUBROUTINE random_positions ( d_min, d_max, r )
    USE utility_module, ONLY : random_orientation_vector
    REAL,                  INTENT(in)  :: d_min, d_max
    REAL, DIMENSION (3,n), INTENT(out) :: r

    INTEGER :: i, j
    REAL    :: d, ran
    LOGICAL :: ok

    ! This routine places n molecules at positions which are all within
    ! the desired range d_min ... d_max of each other.
    ! Obviously this approach can fail, depending on the value of n
    ! and the supplied distances.
    r(:,1) = [0.0,0.0,0.0] ! first one at origin

    ! Choose r vector randomly in appropriate range
    DO i = 2, n
       DO ! loop until successfully 
          CALL random_orientation_vector ( r(:,i) ) ! direction of r
          CALL RANDOM_NUMBER ( ran )
          d      = d_min + (d_max-d_min)*ran ! magnitude of r
          r(:,i) = r(:,i) * d ! within desired range of origin
          ok = .TRUE.
          DO j = 2, i-1 ! check intermediate atoms if any
             d = SQRT ( SUM((r(:,i)-r(:,j))**2) )
             ok = ok .AND. ( d >= d_min ) .AND. ( d <= d_max )
          END DO ! end check intermediate atoms if any
          IF ( ok ) EXIT
       END DO
    END DO

  END SUBROUTINE random_positions

END PROGRAM test_pot


