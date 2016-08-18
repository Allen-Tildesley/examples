! md_lj_mts_module.f90
! Force routine for MD, LJ atoms, multiple timesteps
MODULE md_lj_mts_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, v, f
  PUBLIC :: allocate_arrays, deallocate_arrays, force

  INTEGER                                :: n      ! number of atoms
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: r      ! positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: v      ! velocities (3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: f      ! forces for each shell (3,n,k_max) 

CONTAINS

  SUBROUTINE allocate_arrays ( r_cut )
    REAL, DIMENSION(:), INTENT(in)  :: r_cut  ! shell cutoff distances

    INTEGER :: k_max
    k_max = SIZE(r_cut)
    ALLOCATE ( r(3,n), v(3,n), f(3,n,k_max) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, v, f )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( box, r_cut, lambda, k, pot, vir )
    REAL,                  INTENT(in)  :: box    ! box length
    REAL,    DIMENSION(:), INTENT(in)  :: r_cut  ! shell cutoff distances
    REAL,                  INTENT(in)  :: lambda ! switch function healing length
    INTEGER,               INTENT(in)  :: k      ! shell for this force evaluation
    REAL,                  INTENT(out) :: pot    ! total potential energy for this shell
    REAL,                  INTENT(out) :: vir    ! virial for this shell

    ! Calculates potential, virial and forces
    ! Each shell carries its own contribution to potential and forces
    ! The contribution includes a multiplicative switching function
    ! All quantities are calculated in units where sigma = 1 and epsilon = 1

    INTEGER            :: i, j, k_max
    REAL               :: rij_sq, rij_mag, sr2, sr6, sr12, s, ds, x
    REAL               :: pot_lj, vir_lj, potij, virij, pot_cut
    REAL, DIMENSION(3) :: rij, fij
    REAL               :: rk, rkm, rk_sq, rkm_sq         ! r_cut(k), r_cut(k-1) and squared values
    REAL               :: rk_l, rkm_l, rk_l_sq, rkm_l_sq ! same, but for r_cut(k)-lambda etc

    ! Variables to be returned
    f(:,:,k) = 0.0
    pot      = 0.0
    vir      = 0.0

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

    ! calculate shift in potential at outer cutoff
    sr2     = 1.0 / r_cut(k_max)**2
    sr6     = sr2 ** 3
    sr12    = sr6 ** 2
    pot_cut = 4.0* ( sr12 - sr6 )

    ! Double loop over atoms
    DO i = 1, n-1
       DO j = i+1, n

          rij(:) = r(:,i) - r(:,j)                     ! separation vector
          rij(:) = rij(:) - ANINT ( rij(:)/box ) * box ! periodic boundary conditions
          rij_sq = SUM ( rij**2 )                      ! squared separation

          IF ( rij_sq <= rk_sq .AND. rij_sq >= rkm_l_sq ) THEN ! test whether in shell

             sr2    = 1.0 / rij_sq
             sr6    = sr2 ** 3
             sr12   = sr6 ** 2
             pot_lj = 4.0 *  ( sr12 - sr6 ) - pot_cut    ! Lennard-Jones potential function
             vir_lj = 24.0 * ( 2.0*sr12 - sr6)           ! -rij_mag*derivative of pot_lj

             ! the following statements implement S(k,rij_mag) - S(k-1,rij_mag)
             ! it is assumed that r_cut(k-1) and r_cut(k) are at least lambda apart

             IF ( rij_sq < rkm_sq ) THEN               ! S(k,rij_mag)=1, S(k-1,rij_mag) varying
                rij_mag = SQRT ( rij_sq )              ! rij_mag lies between rkm-lambda and rkm
                x       = (rij_mag-rkm)/lambda         ! x lies between -1 and 0
                s       = 1.0 - (2.0*x+3.0)*x**2       ! 1 - S(k-1,rij_mag) lies between 0 and 1
                ds      = 6.0*(x+1.0)*x*rij_mag/lambda ! -rij_mag*derivative of (1-S(k-1,rij_mag))
                potij   = s * pot_lj                   ! potential includes switching function
                virij   = s * vir_lj + ds * pot_lj     ! virial also includes derivative

             ELSE IF ( rij_sq > rk_l_sq ) THEN ! S(k,rij_mag) varying, S(k-1,rij_mag)=0
                IF ( k == k_max ) THEN         ! no switch at outermost cutoff
                   potij = pot_lj              ! potential unchanged
                   virij = vir_lj              ! virial unchanged
                ELSE
                   rij_mag = SQRT ( rij_sq )               ! rij_mag lies between rk-lambda and rk
                   x       = (rij_mag-rk)/lambda           ! x lies between -1 and 0
                   s       = (2.0*x+3.0)*x**2              ! S(k,rij_mag) lies between 1 and 0
                   ds      = -6.0*(x+1.0)*x*rij_mag/lambda ! -rij_mag*derivative of S(k,rij_mag)
                   potij   = s * pot_lj                    ! potential includes switching function
                   virij   = s * vir_lj + ds * pot_lj      ! virial also includes derivative
                END IF
             ELSE               ! S(k,rij_mag)=1, S(k-1,rij_mag)=0, rij_mag lies between rkm and rk-lambda
                potij = pot_lj  ! potential unchanged
                virij = vir_lj  ! virial unchanged
             END IF

             pot      = pot + potij
             vir      = vir + virij
             virij    = virij / rij_sq
             fij      = rij * virij
             f(:,i,k) = f(:,i,k) + fij
             f(:,j,k) = f(:,j,k) - fij

          END IF ! end test whether in shell

       END DO ! End loop over pairs in this shell
    END DO

    ! multiply virial by numerical factor
    vir = vir / 3.0

  END SUBROUTINE force

END MODULE md_lj_mts_module
