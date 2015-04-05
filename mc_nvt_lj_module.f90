! mc_nvt_lj_module.f90 (used by mc_nvt_lj.f90)
! Monte Carlo simulation, constant-NVT ensemble, Lennard-Jones atoms
MODULE mc_nvt_lj_module

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, lt, ne, gt
  PUBLIC :: energy_1, energy, energy_lrc

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options

CONTAINS

  SUBROUTINE energy ( sigma, r_cut, pot, vir, overlap )
    REAL,    INTENT(in)  :: sigma, r_cut ! potential parameters
    REAL,    INTENT(out) :: pot, vir     ! potential and virial 
    LOGICAL, INTENT(out) :: overlap      ! shows if an overlap was detected

    ! Calculates potential and virial for whole system
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! The routine returns immediately on overlap, in which case
    ! the values of pot and vir should not be used

    REAL    :: pot_i, vir_i
    INTEGER :: i

    overlap = .FALSE.
    pot      = 0.0
    vir      = 0.0

    DO i = 1, n - 1

       CALL energy_1 ( r(:,i), i, gt, sigma, r_cut, pot_i, vir_i, overlap )

       IF ( overlap ) RETURN
       pot = pot + pot_i
       vir = vir + vir_i
    END DO

  END SUBROUTINE energy

  SUBROUTINE energy_1 ( ri, i, j_range, sigma, r_cut, pot, vir, overlap )

    REAL, DIMENSION(3), INTENT(in)  :: ri
    INTEGER,            INTENT(in)  :: i, j_range
    REAL,               INTENT(in)  :: r_cut, sigma
    REAL,               INTENT(out) :: pot, vir
    LOGICAL,            INTENT(out) :: overlap

    ! calculates potential energy and virial of atom i
    ! with j/=i, j>i, or j<i depending on j_range
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! The routine returns immediately on overlap, in which case
    ! the values of pot and vir should not be used

    INTEGER            :: j, j1, j2
    REAL               :: r_cut_sq, sigma_sq
    REAL               :: sr2, sr6, rij_sq
    REAL, DIMENSION(3) :: rij
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

    r_cut_sq = r_cut**2
    sigma_sq = sigma**2

    pot   = 0.0
    vir   = 0.0
    overlap = .FALSE.

    SELECT CASE ( j_range )
    CASE ( lt ) ! j < i
       j1 = 1
       j2 = i-1
    CASE ( gt ) ! j > i
       j1 = i+1
       j2 = n
    CASE ( ne ) ! j /= i
       j1 = 1
       j2 = n
    END SELECT

    DO j = j1, j2

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )

       IF ( rij_sq < r_cut_sq ) THEN

          sr2 = sigma_sq / rij_sq ! now in LJ units

          IF ( sr2 > sr2_overlap ) THEN
             overlap = .TRUE.
             RETURN
          END IF

          sr6 = sr2 * sr2 * sr2
          pot = pot + sr6 * ( sr6 - 1.0 )
          vir = vir + sr6 * ( sr6 - 0.5 )

       END IF

    END DO
    pot = 4.0 * pot        ! LJ units (sigma = 1)
    vir = 48.0 * vir / 3.0 ! LJ units (sigma = 1)

  END SUBROUTINE energy_1

  SUBROUTINE energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc ) ! Long-range corrections
    INTEGER, intent(in)  :: n                ! number of atoms
    REAL,    INTENT(in)  :: sigma, r_cut     ! LJ potential parameters
    REAL,    INTENT(out) :: pot_lrc, vir_lrc ! Results in LJ (sigma=1) units

    REAL            :: sr3, density
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3     = ( sigma / r_cut ) ** 3
    density = REAL(n)*sigma**3
    pot_lrc = REAL(n) * (8.0/9.0)  * pi * density * ( sr3**3 - 3.0*sr3 )
    vir_lrc = REAL(n) * (32.0/9.0) * pi * density * ( sr3**3 - 1.5*sr3 )

  END SUBROUTINE energy_lrc

END MODULE mc_nvt_lj_module
