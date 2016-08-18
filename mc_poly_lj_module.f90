! mc_poly_lj_module.f90
! Routines for MC simulation, polyatomic molecule, LJ atoms
MODULE mc_poly_lj_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, na, r, e, d, lt, ne, gt
  PUBLIC :: allocate_arrays, deallocate_arrays, energy_1, energy, q_to_d
 
  INTEGER                                :: n  ! number of molecules
  INTEGER                                :: na ! number of atoms per molecule
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: r  ! centre of mass positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: e  ! quaternions (0:3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: d  ! bond vectors (3,na,n)

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options

CONTAINS

  SUBROUTINE allocate_arrays
    ALLOCATE ( r(3,n), e(0:3,n), d(3,na,n) )
  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, e, d )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE energy ( sigma, r_cut, rm_cut, overlap, pot, vir )
    REAL,    INTENT(in)  :: sigma, r_cut ! potential parameters
    REAL,    INTENT(in)  :: rm_cut       ! molecule-molecule cutoff distance
    LOGICAL, INTENT(out) :: overlap      ! shows if an overlap was detected
    REAL,    INTENT(out) :: pot, vir     ! potential and virial 

    ! Calculates potential and virial for whole system
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the values of pot and vir should not be used
    ! It is assumed that r, sigma, r_cut and rm_cut are in units where box = 1
    ! Results are in LJ units where sigma = 1, epsilon = 1

    REAL    :: pot_i, vir_i, pot_sum, vir_sum
    INTEGER :: i

    IF ( SIZE(r,dim=2)  /= n ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy'
    END IF
    
    overlap  = .FALSE.
    pot_sum  = 0.0
    vir_sum  = 0.0

    DO i = 1, n - 1
       CALL energy_1 ( r(:,i), d(:,:,i), i, gt, sigma, r_cut, rm_cut, overlap, pot_i, vir_i )
       IF ( overlap ) EXIT ! jump out of loop
       pot_sum  = pot_sum  + pot_i
       vir_sum  = vir_sum  + vir_i
    END DO

    pot  = pot_sum
    vir  = vir_sum
    
  END SUBROUTINE energy

  SUBROUTINE energy_1 ( ri, di, i, j_range, sigma, r_cut, rm_cut, overlap, pot, vir )

    REAL,    INTENT(in), DIMENSION(3)   :: ri           ! position of i, may differ from r(:,i)
    REAL,    INTENT(in), DIMENSION(:,:) :: di           ! bond vectors of i, may differ from d(:,:,i)
    INTEGER, INTENT(in)                 :: i, j_range   ! index, and partner index range
    REAL,    INTENT(in)                 :: r_cut, sigma ! LJ potential parameters
    REAL,    INTENT(in)                 :: rm_cut       ! molecule-molecule cutoff distance
    LOGICAL, INTENT(out)                :: overlap      ! shows if an overlap was detected
    REAL,    INTENT(out)                :: pot, vir     ! potential and virial

    ! Calculates potential energy and virial of molecule i
    ! with j/=i, j>i, or j<i depending on j_range
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the values of pot and vir should not be used
    ! It is assumed that r, sigma, r_cut and rm_cut are in units where box = 1
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: j, j1, j2, a, b, n_cut
    REAL               :: r_cut_sq, rm_cut_sq, sigma_sq
    REAL               :: sr2, sr6, sr12, rij_sq, rab_sq, potab, virab
    REAL, DIMENSION(3) :: rij, rab, fab
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

    IF ( SIZE(r,dim=2)  /= n ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy_1'
    END IF
    IF ( SIZE(di,dim=1) /= 3 ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Array bounds error for di', SIZE(di,dim=1)
       STOP 'Error in energy_1'
    END IF
    IF ( SIZE(di,dim=2) /= na ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for di', SIZE(di,dim=2), na
       STOP 'Error in energy_1'
    END IF

    r_cut_sq  = r_cut**2
    rm_cut_sq = rm_cut**2
    sigma_sq  = sigma**2

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

    pot = 0.0
    vir = 0.0
    n_cut = 0

    j_loop: DO j = j1, j2 ! loop over j-molecules

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )

       IF ( rij_sq < rm_cut_sq ) THEN ! test within molecular cutoff

          ! double loop over atoms on both molecules
          DO a = 1, na
             DO b = 1, na
                rab = rij + di(:,a) - d(:,b,j) ! atom-atom vector
                rab_sq = SUM ( rab**2 )

                IF ( rab_sq < r_cut_sq ) THEN ! test within potential cutoff 

                   sr2 = sigma_sq / rab_sq ! now dimensionless

                   IF ( sr2 > sr2_overlap ) THEN
                      overlap = .TRUE.
                      EXIT j_loop ! jump out of outer j loop
                   END IF

                   sr6   = sr2**3
                   sr12  = sr6**2
                   potab = sr12 - sr6
                   virab = potab + sr12
                   pot   = pot + potab
                   fab   = rab * virab / rab_sq
                   vir   = vir + DOT_PRODUCT ( rij, fab ) ! molecular virial
                   n_cut = n_cut + 1
                END IF ! end test within potential cutoff

             END DO
          END DO
          ! end double loop over atoms on both molecules

       END IF ! end test within molecular cutoff

    END DO j_loop ! end loop over j-molecules

    ! shift potentials at cutoff
    sr2 = sigma_sq / r_cut_sq
    sr6 = sr2**3
    sr12 = sr6**2
    potab  = sr12 - sr6
    pot = pot - potab * REAL(n_cut)

    ! scale potential and virial
    pot = pot * 4.0
    vir = vir * 24.0 / 3.0

  END SUBROUTINE energy_1

  FUNCTION q_to_d ( q, db ) RESULT ( ds ) ! convert quaternion to bond vectors in space-fixed frame
    USE utility_module, ONLY : q_to_a
    IMPLICIT NONE
    REAL, DIMENSION(0:3),                          INTENT(in) :: q  ! quaternion
    REAL, DIMENSION(:,:),                          INTENT(in) :: db ! body-fixed bonds
    REAL, DIMENSION(SIZE(db,dim=1),SIZE(db,dim=2))            :: ds ! space-fixed bonds

    REAL, DIMENSION(3,3) :: a ! rotation matrix
    INTEGER              :: i

    IF ( SIZE(db,dim=1) /= 3  ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Array mismatch', SIZE(db,dim=1)
       STOP 'Error in q_to_d'
    END IF
    IF ( SIZE(db,dim=2) /= na ) then ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array mismatch', SIZE(db,dim=2), na
       STOP 'Error in q_to_d'
    END IF
    a = q_to_a ( q )
    DO i = 1, na
       ds(:,i) = MATMUL ( db(:,i), a ) ! NB: using transpose of rotation matrix
    END DO

  END FUNCTION q_to_d
  
END MODULE mc_poly_lj_module
