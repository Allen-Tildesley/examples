! mc_poly_lj_module.f90
! Routines for MC simulation, polyatomic molecule, LJ atoms
MODULE mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, na, r, e, d, lt, ne, gt
  PUBLIC :: model_description, allocate_arrays, deallocate_arrays
  public :: energy_1, energy, q_to_d
 
  INTEGER                                :: n  ! number of molecules
  INTEGER                                :: na ! number of atoms per molecule
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: r  ! centre of mass positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: e  ! quaternions (0:3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: d  ! bond vectors (3,na,n)

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options
  REAL,    PARAMETER :: sigma = 1.0             ! Lennard-Jones diameter (unit of length)
  REAL,    PARAMETER :: epslj = 1.0             ! Lennard-Jones well depth (unit of energy)

CONTAINS

  SUBROUTINE model_description ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Shifted Lennard-Jones potential (no LRC)'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ',   sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epslj = ', epslj    
  END SUBROUTINE model_description

  SUBROUTINE allocate_arrays ( box, rm_cut )
    REAL, INTENT(in) :: box    ! simulation box length
    REAL, INTENT(in) :: rm_cut ! potential molecule cutoff distance

    REAL :: rm_cut_box
    
    ALLOCATE ( r(3,n), e(0:3,n), d(3,na,n) )

    rm_cut_box = rm_cut / box
    IF ( rm_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)') 'rm_cut/box too large ', rm_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    DEALLOCATE ( r, e, d )
  END SUBROUTINE deallocate_arrays

  SUBROUTINE energy ( box, r_cut, rm_cut, overlap, pot, vir )
    REAL,    INTENT(in)  :: box      ! simulation box length
    REAL,    INTENT(in)  :: r_cut    ! potential cutoff distance
    REAL,    INTENT(in)  :: rm_cut   ! molecule-molecule cutoff distance
    LOGICAL, INTENT(out) :: overlap  ! shows if an overlap was detected
    REAL,    INTENT(out) :: pot, vir ! potential and virial 

    ! Calculates potential and virial for whole system
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the values of pot and vir should not be used
    ! Actual calculation is performed by subroutine energy_1

    REAL    :: pot_i, vir_i
    INTEGER :: i

    IF ( SIZE(r,dim=2)  /= n ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy'
    END IF
    
    overlap = .FALSE.
    pot     = 0.0
    vir     = 0.0

    DO i = 1, n - 1
       CALL energy_1 ( r(:,i), d(:,:,i), i, gt, box, r_cut, rm_cut, overlap, pot_i, vir_i )
       IF ( overlap ) EXIT ! jump out of loop
       pot = pot + pot_i
       vir = vir + vir_i
    END DO
    
  END SUBROUTINE energy

  SUBROUTINE energy_1 ( ri, di, i, j_range, box, r_cut, rm_cut, overlap, pot, vir )

    REAL,    INTENT(in), DIMENSION(3)   :: ri         ! position of i, may differ from r(:,i)
    REAL,    INTENT(in), DIMENSION(:,:) :: di         ! bond vectors of i, may differ from d(:,:,i)
    INTEGER, INTENT(in)                 :: i, j_range ! index, and partner index range
    REAL,    INTENT(in)                 :: box        ! simulation box length
    REAL,    INTENT(in)                 :: r_cut      ! potential cutoff distance
    REAL,    INTENT(in)                 :: rm_cut     ! molecule-molecule cutoff distance
    LOGICAL, INTENT(out)                :: overlap    ! shows if an overlap was detected
    REAL,    INTENT(out)                :: pot, vir   ! potential and virial

    ! Calculates potential energy and virial of molecule i
    ! with j/=i, j>i, or j<i depending on j_range
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the values of pot and vir should not be used
    ! It is assumed that r and d have been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: j, j1, j2, a, b, n_cut
    REAL               :: r_cut_box, r_cut_box_sq, rm_cut_box, rm_cut_box_sq, box_sq
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

    r_cut_box     = r_cut / box
    r_cut_box_sq  = r_cut_box**2
    rm_cut_box    = rm_cut / box
    rm_cut_box_sq = rm_cut_box**2
    box_sq        = box**2

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

       IF ( rij_sq < rm_cut_box_sq ) THEN ! test within molecular cutoff

          ! double loop over atoms on both molecules
          DO a = 1, na
             DO b = 1, na
                rab = rij + di(:,a) - d(:,b,j) ! atom-atom vector, box=1 units
                rab_sq = SUM ( rab**2 )

                IF ( rab_sq < r_cut_box_sq ) THEN ! test within potential cutoff 

                   rab_sq = rab_sq * box_sq ! now in sigma=1 units
                   sr2    = 1.0 / rab_sq    ! (sigma/rab)**2

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
    sr2   = 1.0 / r_cut**2 ! in sigma=1 units
    sr6   = sr2**3
    sr12  = sr6**2
    potab = sr12 - sr6
    pot   = pot - potab * REAL(n_cut)

    ! scale potential and virial
    pot = pot * 4.0
    vir = vir * 24.0 / 3.0

  END SUBROUTINE energy_1

  FUNCTION q_to_d ( q, db ) RESULT ( ds ) ! convert quaternion to bond vectors in space-fixed frame
    USE maths_module, ONLY : q_to_a
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

    a = q_to_a ( q ) ! compute rotation matrix from quaternions
    DO i = 1, na
       ds(:,i) = MATMUL ( db(:,i), a ) ! NB: equivalent to ds = at*db, at=transpose of a
    END DO

  END FUNCTION q_to_d
  
END MODULE mc_module
