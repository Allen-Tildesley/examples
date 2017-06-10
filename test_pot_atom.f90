! test_pot_atom.f90
! Test potential, forces for atoms
PROGRAM test_pot_atom

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

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               test_pot_module, ONLY : n, force
  USE               maths_module,    ONLY : init_random_seed, cross_product

  IMPLICIT NONE

  REAL, DIMENSION(3,n) :: r     ! Positions
  REAL, DIMENSION(3,n) :: f     ! Forces
  REAL, DIMENSION(3)   :: axis  ! Axis of translation
  REAL, DIMENSION(3)   :: rsave ! Saves reference position
  REAL, DIMENSION(3)   :: tot   ! Total force or torque to check conservation

  LOGICAL :: ok
  INTEGER :: i, xyz, try, ntry, npos, ioerr
  REAL    :: pot, potp, potm, fnum
  REAL    :: delta, d_min, d_max, pot_max

  CHARACTER(len=2), DIMENSION(3), PARAMETER :: cf = ['Fx','Fy','Fz']

  NAMELIST /nml/ delta, d_min, d_max, pot_max, ntry, npos

  ! Initialize random number generator (hopefully different every time!)
  CALL init_random_seed

  ! Default values: any of the following could be empirically adjusted
  delta   = 1.e-5 ! Small displacement
  d_min   = 0.3   ! Minimum separation between atoms
  d_max   = 1.5   ! Maximum separation between atoms
  pot_max = 10.0  ! Maximum potential to allow in atom placement
  ntry    = 1000  ! Number of attempts to make in order to place atoms
  npos    = 1000  ! Number of attempts to position each atom

  ! Read parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in test_pot_atom'
  END IF

  ! Write out parameters
  WRITE ( unit=output_unit, fmt='(a,t40,es15.4)' ) 'Displacement delta',      delta
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'  ) 'Min separation d_min',    d_min
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'  ) 'Max separation d_max',    d_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'  ) 'Max potential pot_max',   pot_max
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'    ) 'Max placement tries',     ntry
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'    ) 'Max atom position tries', npos

  ! Make a number of attempts to place the atoms
  ok = .FALSE.
  DO try = 1, ntry 
     CALL random_positions ( d_min, d_max, npos, r )
     CALL force ( r, pot )
     IF ( ABS(pot) < pot_max ) THEN
        ok = .TRUE.
        EXIT
     END IF
  END DO
  IF ( .NOT. ok ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Exceeded allowed number of tries'
     STOP 'Error in test_pot_atom'
  END IF

  CALL force ( r, pot, f ) ! Calculation of potential, and analytical forces
  WRITE ( output_unit, fmt='(/,a,t40,f15.6)' ) 'Potential energy', pot

  ! Check momentum and angular momentum conservation
  tot = SUM ( f, dim=2 )
  WRITE ( output_unit, fmt='(a,t40,3es15.4)') 'Total force', tot
  tot = 0.0
  DO i = 1, n
     tot = tot + cross_product ( r(:,i), f(:,i) )
  END DO
  WRITE ( output_unit, fmt='(a,t40,3es15.4)') 'Total torque', tot

  WRITE ( output_unit, fmt='(/,2a15,t40,3a15)') 'Atom', 'Component', 'Exact', 'Numerical', 'Difference'

  DO i = 1, n ! Loop over atoms

     DO xyz = 1, 3 ! Loop over Cartesian components

        ! Pick axis
        axis      = 0.0
        axis(xyz) = 1.0

        rsave  = r(:,i) ! Save position
        r(:,i) = rsave + delta*axis ! Translate
        CALL force ( r, potp )
        r(:,i) = rsave - delta*axis ! Translate
        CALL force ( r, potm )
        r(:,i) = rsave ! Restore position
        fnum   = -(potp-potm)/(2.0*delta)

        WRITE ( output_unit, fmt='(i15,a15,t40,2f15.6,es15.4)' ) i, cf(xyz), f(xyz,i), fnum, f(xyz,i)-fnum

     END DO ! End loop over Cartesian components

  END DO ! End loop over atoms

CONTAINS

  SUBROUTINE random_positions ( d_min, d_max, npos, r )
    USE maths_module, ONLY : random_vector
    IMPLICIT NONE
    REAL,                  INTENT(in)  :: d_min, d_max
    INTEGER,               INTENT(in)  :: npos
    REAL, DIMENSION (3,n), INTENT(out) :: r

    INTEGER :: i, j, pos_try
    REAL    :: d, zeta
    LOGICAL :: ok

    ! This routine places n atoms at positions which are all within
    ! the desired range d_min ... d_max of each other.
    ! Obviously this approach can fail, depending on the value of n
    ! and the supplied distances.

    r(:,1) = [0.0,0.0,0.0] ! First one at origin

    DO i = 2, n ! Loop over remaining atoms

       pos_try = 0
       DO ! Loop until successful
          r(:,i) = random_vector ( ) ! Direction of r
          CALL RANDOM_NUMBER ( zeta )
          d      = d_min + (d_max-d_min)*zeta ! Magnitude of r
          r(:,i) = r(:,i) * d                 ! Within desired range of origin
          ok = .TRUE.
          DO j = 2, i-1 ! Check intermediate atoms if any
             d  = SQRT ( SUM((r(:,i)-r(:,j))**2) )
             ok = ok .AND. ( d >= d_min ) .AND. ( d <= d_max )
          END DO ! end check intermediate atoms if any
          IF ( ok ) EXIT
          pos_try = pos_try + 1
          IF ( pos_try > npos ) THEN
             WRITE ( unit=error_unit, fmt='(a)') 'Exceeded allowed number of tries'
             STOP 'Error in random_positions'
          END IF
       END DO ! End loop until successful

    END DO ! End loop over remaining atoms

  END SUBROUTINE random_positions

END PROGRAM test_pot_atom
