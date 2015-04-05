! mc_npt_lj_module.f90 (used by mc_npt_lj.f90)
! Monte Carlo simulation, constant-NPT ensemble, Lennard-Jones atoms
MODULE mc_npt_lj_module

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, lt, ne, gt
  PUBLIC :: energy_1, energy, pot_lrc, vir_lrc

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options

CONTAINS

  SUBROUTINE energy ( sigma, r_cut, pot, vir, overlap )
    REAL,               INTENT(in)  :: sigma, r_cut ! potential parameters
    REAL, DIMENSION(2), INTENT(out) :: pot, vir     ! potential and virial 
    LOGICAL,            INTENT(out) :: overlap      ! shows if an overlap was detected

    ! Calculates potential and virial for whole system
    ! pot and vir are both returned as two components: LJ12 and LJ6 
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! The routine returns immediately on overlap, in which case
    ! the values of pot and vir should not be used

    REAL, dimension(2) :: pot_i, vir_i
    INTEGER            :: i

    overlap = .FALSE.
    pot     = 0.0
    vir     = 0.0

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
    REAL, DIMENSION(2), INTENT(out) :: pot, vir
    LOGICAL,            INTENT(out) :: overlap

    ! calculates potential energy and virial of atom i
    ! with j/=i, j>i, or j<i depending on j_range
    ! pot and vir are both returned as two components: LJ12 and LJ6 
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
             return
          END IF

          sr6    = sr2**3
          pot(1) = pot(1) + sr6**2
          pot(2) = pot(2) - sr6

       END IF

    END DO
    pot = 4.0 * pot              ! LJ units (sigma = 1)
    vir(1) = 12.0 * pot(1) / 3.0 ! LJ units (sigma = 1)
    vir(2) = 6.0  * pot(2) / 3.0 ! LJ units (sigma = 1)

  END SUBROUTINE energy_1

  FUNCTION pot_lrc ( sigma, r_cut, density ) ! Long-range correction to potential per atom
    REAL, dimension(2) :: pot_lrc            ! Function result in LJ (sigma=1) units
    REAL, INTENT(in)   :: sigma, r_cut, density

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3     = ( sigma / r_cut ) ** 3
    pot_lrc(1) =  (8.0/9.0) * pi * density * sr3**3
    pot_lrc(2) = -(8.0/3.0) * pi * density * sr3

  END FUNCTION pot_lrc

  FUNCTION vir_lrc ( sigma, r_cut, density ) ! Long-range correction to virial per atom
    REAL, dimension(2) :: vir_lrc            ! Function result in LJ (sigma=1) units
    REAL, INTENT(in)   :: sigma, r_cut, density

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3     = ( sigma / r_cut ) ** 3
    vir_lrc(1) =  (32.0/9.0) * pi * density * sr3**3
    vir_lrc(2) = -(32.0/6.0) * pi * density * sr3

  END FUNCTION vir_lrc

END MODULE mc_npt_lj_module
