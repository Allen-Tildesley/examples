! md_lj_mts_module.f90
! Force routine for MD, LJ atoms, multiple timesteps
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
  PUBLIC :: force, hessian

  ! Public data
  INTEGER,                                PUBLIC :: n ! Number of atoms
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: r ! Positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: v ! Velocities (3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: f ! Forces for each shell (3,n,k_max) 

  ! Public derived type
  TYPE, PUBLIC :: potential_type ! A composite variable for interactions comprising
     REAL    :: cut ! the potential energy cut (but not shifted) at r_cut and
     REAL    :: pot ! the potential energy cut-and-shifted at r_cut and
     REAL    :: vir ! the virial and
     REAL    :: lap ! the Laplacian and
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
    c%cut = a%cut  +   b%cut
    c%pot = a%pot  +   b%pot
    c%vir = a%vir  +   b%vir
    c%lap = a%lap  +   b%lap
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut-and-shifted version for dynamics'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut (but not shifted) version also calculated'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Calculated in shells with switching functions'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'   
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'  

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( r_cut )
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(in)  :: r_cut  ! shell cutoff distances

    INTEGER :: k_max
    k_max = SIZE(r_cut)
    ALLOCATE ( r(3,n), v(3,n), f(3,n,k_max) )

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, v, f )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( box, r_cut, lambda, k, total )
    IMPLICIT NONE
    REAL,                 INTENT(in)  :: box    ! Box length
    REAL, DIMENSION(:),   INTENT(in)  :: r_cut  ! Array of shell cutoff distances
    REAL,                 INTENT(in)  :: lambda ! Switch function healing length
    INTEGER,              INTENT(in)  :: k      ! Shell for this force evaluation
    TYPE(potential_type), INTENT(out) :: total  ! Composite of pot, vir, lap etc

    ! total%pot is the nonbonded cut-and-shifted potential energy for whole system
    ! total%cut is the nonbonded cut (but not shifted) potential energy for whole system
    ! total%vir is the corresponding virial
    ! total%lap is the corresponding Laplacian
    ! total%ovr is a warning flag that there is an overlap
    ! This routine also calculates forces and stores them in the array f
    ! Forces are derived from pot, not cut (which has a discontinuity)
    ! If total%ovr is set to .true., the forces etc should not be used

    ! Only separations lying in specified shell are considered
    ! The contributions include a multiplicative switching function
    ! NB total%vir includes derivative of switching function, because it is used in the 
    ! calculation of forces: when we sum over k to get pressure, these extra terms cancel
    ! total%lap does not include terms of this kind because we only need to sum over k to
    ! get the total at the end of a long step, and the extra terms would all vanish

    ! All quantities are calculated in units where sigma = 1 and epsilon = 1
    ! Positions are assumed to be in these units as well

    INTEGER              :: i, j, k_max
    REAL                 :: rij_sq, rij_mag, sr2, sr6, sr12, s, ds, x
    REAL                 :: pot_lj, cut_lj, vir_lj, lap_lj, pot_cut
    REAL, DIMENSION(3)   :: rij, fij
    REAL                 :: rk, rkm, rk_sq, rkm_sq         ! r_cut(k), r_cut(k-1) and squared values
    REAL                 :: rk_l, rkm_l, rk_l_sq, rkm_l_sq ! Same, but for r_cut(k)-lambda etc
    REAL, PARAMETER      :: sr2_ovr = 1.77                 ! Overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    ! Distances for switching function
    k_max = SIZE(r_cut)
    IF ( k < 1 .OR. k > k_max ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'k, k_max error', k, k_max
       STOP 'Error in force'
    END IF

    rk      = r_cut(k)
    rk_l    = rk-lambda
    rk_sq   = rk ** 2
    rk_l_sq = rk_l ** 2
    IF ( k == 1 ) THEN
       rkm   = 0.0
       rkm_l = 0.0
    ELSE
       rkm   = r_cut(k-1)
       rkm_l = rkm-lambda
    END IF
    rkm_sq   = rkm ** 2
    rkm_l_sq = rkm_l ** 2

    ! Calculate shift in potential at outermost cutoff
    sr2     = 1.0 / r_cut(k_max)**2
    sr6     = sr2 ** 3
    sr12    = sr6 ** 2
    pot_cut = 4.0 * ( sr12 - sr6 ) ! NB we include numerical factor

    ! Initialize
    f(:,:,k) = 0.0
    total = potential_type ( pot=0.0, cut=0.0, vir=0.0, lap=0.0, ovr=.FALSE. )

    ! Double loop over atoms
    DO i = 1, n-1
       DO j = i+1, n

          rij(:) = r(:,i) - r(:,j)                     ! Separation vector
          rij(:) = rij(:) - ANINT ( rij(:)/box ) * box ! Periodic boundary conditions
          rij_sq = SUM ( rij**2 )                      ! Squared separation

          IF ( rij_sq <= rk_sq .AND. rij_sq >= rkm_l_sq ) THEN ! Test whether in shell

             sr2 = 1.0 / rij_sq       ! (sigma/rij)**2
             pair%ovr = sr2 > sr2_ovr ! Overlap if too close

             sr6    = sr2 ** 3
             sr12   = sr6 ** 2
             cut_lj = 4.0 *  ( sr12 - sr6 )                ! LJ cut (but not shifted) potential function
             pot_lj = cut_lj - pot_cut                     ! LJ cut-and-shifted potential function
             vir_lj = 24.0 * ( 2.0*sr12 - sr6)             ! -rij_mag*derivative of pot_lj
             lap_lj = 24.0 * ( 22.0*sr12 - 5.0*sr6 ) * sr2 ! LJ pair Laplacian

             ! The following statements implement the switching function S(k,rij_mag) - S(k-1,rij_mag)
             ! It is assumed that r_cut(k-1) and r_cut(k) are at least lambda apart

             IF ( rij_sq < rkm_sq ) THEN ! S(k,rij_mag)=1, S(k-1,rij_mag) varying
                rij_mag  = SQRT ( rij_sq )              ! rij_mag lies between rkm-lambda and rkm
                x        = (rij_mag-rkm)/lambda         ! x lies between -1 and 0
                s        = 1.0 - (2.0*x+3.0)*x**2       ! 1 - S(k-1,rij_mag) lies between 0 and 1
                ds       = 6.0*(x+1.0)*x*rij_mag/lambda ! -rij_mag*derivative of (1-S(k-1,rij_mag))
                pair%pot = s * pot_lj                   ! Potential includes switching function
                pair%vir = s * vir_lj + ds * pot_lj     ! Virial also includes derivative
                pair%cut = s * cut_lj                   ! Cut potential includes switching function
                pair%lap = s * lap_lj                   ! Laplacian includes switching function

             ELSE IF ( rij_sq > rk_l_sq ) THEN ! S(k,rij_mag) varying, S(k-1,rij_mag)=0

                IF ( k == k_max ) THEN ! No switch at outermost cutoff
                   pair%pot = pot_lj   ! Potential unchanged
                   pair%vir = vir_lj   ! Virial unchanged
                   pair%cut = cut_lj   ! Cut potential unchanged
                   pair%lap = lap_lj   ! Laplacian unchanged

                ELSE ! the switch function applies
                   rij_mag  = SQRT ( rij_sq )               ! rij_mag lies between rk-lambda and rk
                   x        = (rij_mag-rk)/lambda           ! x lies between -1 and 0
                   s        = (2.0*x+3.0)*x**2              ! S(k,rij_mag) lies between 1 and 0
                   ds       = -6.0*(x+1.0)*x*rij_mag/lambda ! -rij_mag*derivative of S(k,rij_mag)
                   pair%pot = s * pot_lj                    ! Potential includes switching function
                   pair%vir = s * vir_lj + ds * pot_lj      ! Virial also includes derivative
                   pair%cut = s * cut_lj                    ! Cut potential includes switching function
                   pair%lap = s * lap_lj                    ! Laplacian includes switching function
                END IF

             ELSE ! S(k,rij_mag)=1, S(k-1,rij_mag)=0, rij_mag lies between rkm and rk-lambda
                pair%pot = pot_lj  ! Potential unchanged
                pair%vir = vir_lj  ! Virial unchanged
                pair%cut = cut_lj  ! Cut potential unchanged
                pair%lap = lap_lj  ! Laplacian unchanged
             END IF

             fij = rij * pair%vir / rij_sq ! Pair force

             total    = total    + pair
             f(:,i,k) = f(:,i,k) + fij
             f(:,j,k) = f(:,j,k) - fij

          END IF ! End test whether in shell

       END DO
    END DO
    ! End double loop over atoms

    total%vir = total%vir / 3.0 ! Divide virial by 3
    total%lap = total%lap * 2.0 ! Factor 2 for ij and ji

  END SUBROUTINE force

  FUNCTION hessian ( box, r_cut ) RESULT ( hes )
    IMPLICIT NONE
    REAL             :: hes   ! Returns the total Hessian
    REAL, INTENT(in) :: box   ! Simulation box length
    REAL, INTENT(in) :: r_cut ! Potential cutoff distance

    ! Calculates Hessian function (for 1/N correction to config temp)
    ! This routine is only needed in a constant-energy ensemble
    ! The result is given in units where sigma = 1 and epsilon = 1
    ! It is assumed that positions are in those units too
    ! It is assumed that forces for all shells have already been calculated in array f(3,n,:)
    ! These need to be summed over the last index to get the total fij between two atoms

    INTEGER            :: i, j
    REAL               :: r_cut_sq, rij_sq
    REAL               :: sr2, sr6, sr8, sr10, rf, ff, v1, v2
    REAL, DIMENSION(3) :: rij, fij

    r_cut_sq = r_cut ** 2

    hes = 0.0

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO j = i + 1, n ! Begin inner loop over atoms

          rij(:) = r(:,i) - r(:,j)                       ! Separation vector
          rij(:) = rij(:) - ANINT ( rij(:) / box ) * box ! Periodic boundary conditions
          rij_sq = SUM ( rij**2 )                        ! Squared separation

          IF ( rij_sq < r_cut_sq ) THEN ! Check within cutoff

             fij(:) = SUM ( f(:,i,:) - f(:,j,:), dim=2 ) ! Difference in forces

             ff   = DOT_PRODUCT(fij,fij)
             rf   = DOT_PRODUCT(rij,fij)
             sr2  = 1.0 / rij_sq
             sr6  = sr2 ** 3
             sr8  = sr6 * sr2
             sr10 = sr8 * sr2
             v1   = 24.0 * ( 1.0 - 2.0 * sr6 ) * sr8
             v2   = 96.0 * ( 7.0 * sr6 - 2.0 ) * sr10
             hes  = hes + v1 * ff + v2 * rf**2

          END IF ! End check within cutoff

       END DO ! End inner loop over atoms

    END DO ! End outer loop over atoms

  END FUNCTION hessian

END MODULE md_module
