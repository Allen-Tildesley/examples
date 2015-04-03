MODULE mc_nvt_lj_module

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, j_lt_i, j_ne_i, j_gt_i
  PUBLIC :: energy_i, energy, pot_lrc, vir_lrc

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions

  INTEGER, PARAMETER :: j_lt_i = -1, j_ne_i = 0, j_gt_i = 1 ! j-range options

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

       CALL energy_i ( r(:,i), i, j_gt_i, sigma, r_cut, pot_i, vir_i, overlap )

       IF ( overlap ) RETURN
       pot = pot + pot_i
       vir = vir + vir_i
    END DO

  END SUBROUTINE energy

  SUBROUTINE energy_i ( ri, i, j_range, sigma, r_cut, pot, vir, overlap )

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
    CASE ( j_lt_i )
       j1 = 1
       j2 = i-1
    CASE ( j_gt_i )
       j1 = i+1
       j2 = n
    CASE ( j_ne_i )
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

          sr6 = sr2 * sr2 * sr2
          pot = pot + sr6 * ( sr6 - 1.0 )
          vir = vir + sr6 * ( sr6 - 0.5 )

       END IF

    END DO
    pot = 4.0 * pot        ! LJ units (sigma = 1)
    vir = 48.0 * vir / 3.0 ! LJ units (sigma = 1)

  END SUBROUTINE energy_i

  FUNCTION pot_lrc ( sigma, r_cut, density ) ! Long-range correction to potential per atom
    REAL             :: pot_lrc              ! Function result in LJ (sigma=1) units
    REAL, INTENT(in) :: sigma, r_cut, density

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3     = ( sigma / r_cut ) ** 3
    pot_lrc = (8.0/9.0) * pi * density * ( sr3**3 - 3.0*sr3 )

  END FUNCTION pot_lrc

  FUNCTION vir_lrc ( sigma, r_cut, density ) ! Long-range correction to virial per atom
    REAL             :: vir_lrc              ! Function result in LJ (sigma=1) units
    REAL, INTENT(in) :: sigma, r_cut, density

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3     = ( sigma / r_cut ) ** 3
    vir_lrc = (32.0/9.0) * pi * density * ( sr3**3 - 1.5*sr3 )

  END FUNCTION vir_lrc

END MODULE mc_nvt_lj_module
