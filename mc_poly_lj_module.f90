! mc_poly_lj_module.f90
! Routines for MC simulation, polyatomic molecule, LJ atoms
MODULE mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, na, r, e, d
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: energy_1, energy, q_to_d
  PUBLIC :: potovr

  INTEGER                                :: n  ! number of molecules
  INTEGER                                :: na ! number of atoms per molecule
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: r  ! centre of mass positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: e  ! quaternions (0:3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: d  ! bond vectors (3,na,n)

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range
  REAL,    PARAMETER :: sigma = 1.0     ! Lennard-Jones diameter (unit of length)
  REAL,    PARAMETER :: epslj = 1.0     ! Lennard-Jones well depth (unit of energy)

  TYPE potovr ! A composite variable for interaction energies comprising
     REAL    :: pot ! the potential energy and
     REAL    :: vir ! the virial and
     LOGICAL :: ovr ! a flag indicating overlap (i.e. pot too high to use)
  END TYPE potovr

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_potovr
  END INTERFACE OPERATOR (+)

CONTAINS

  FUNCTION add_potovr ( a, b ) RESULT (c)
    TYPE(potovr)             :: c    ! Result is the sum of the two inputs
    TYPE(potovr), INTENT(in) :: a, b
    c%pot = a%pot +    b%pot
    c%vir = a%vir +    b%vir
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potovr

  SUBROUTINE introduction ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Shifted Lennard-Jones potential (no LRC)'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ',   sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epslj = ', epslj    
  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    INTEGER, INTENT(in) :: output_unit ! unit for standard output
    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'
  END SUBROUTINE conclusion

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

  FUNCTION energy ( box, r_cut, rm_cut )
    TYPE(potovr)     :: energy ! Returns a composite of pot, vir and ovr
    REAL, INTENT(in) :: box    ! Simulation box length
    REAL, INTENT(in) :: r_cut  ! Potential cutoff distance
    REAL, INTENT(in) :: rm_cut ! Molecule-molecule cutoff distance

    ! energy%pot is the nonbonded potential energy for whole system
    ! energy%vir is the corresponding virial for whole system
    ! energy%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of energy%pot, energy%vir should not be used
    ! Actual calculation is performed by function energy_1

    TYPE(potovr) :: energy_i
    INTEGER      :: i

    IF ( SIZE(r,dim=2)  /= n ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy'
    END IF

    energy = potovr ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO i = 1, n - 1
       energy_i = energy_1 ( r(:,i), d(:,:,i), i, box, r_cut, rm_cut, gt )
       IF ( energy_i%ovr ) THEN
          energy%ovr = .TRUE. ! Overlap detected
          RETURN              ! Return immediately
       END IF
       energy = energy + energy_i
    END DO

    energy%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION energy

  FUNCTION energy_1 ( ri, di, i, box, r_cut, rm_cut, j_range ) RESULT ( energy )
    TYPE(potovr)                        :: energy  ! Returns a composite of pot, vir and ovr
    REAL,    DIMENSION(3),   INTENT(in) :: ri      ! Coordinates of molecule of interest
    REAL,    DIMENSION(:,:), INTENT(in) :: di      ! Bond vectors of molecule of interest
    INTEGER,                 INTENT(in) :: i       ! Index of molecule of interest
    REAL,                    INTENT(in) :: box     ! Simulation box length
    REAL,                    INTENT(in) :: r_cut   ! Potential cutoff distance
    REAL,                    INTENT(in) :: rm_cut  ! Molecule-molecule cutoff distance
    INTEGER, OPTIONAL,       INTENT(in) :: j_range ! Optional partner index range

    ! energy%pot is the nonbonded potential energy of molecule ri/di with a set of other molecules
    ! energy%vir is the corresponding virial of molecule ri/di
    ! energy%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the value of energy%pot should not be used
    ! The coordinates in ri and di are not necessarily identical with those in r(:,i) and d(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that r and d have been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1
    ! Note that this is the shifted LJ potential

    INTEGER            :: j, j1, j2, a, b
    REAL               :: r_cut_box, r_cut_box_sq, rm_cut_box, rm_cut_box_sq, box_sq
    REAL               :: sr2, sr6, sr12, rij_sq, rab_sq, potab, virab, pot_cut
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

    IF ( PRESENT ( j_range ) ) THEN
       SELECT CASE ( j_range )
       CASE ( lt ) ! j < i
          j1 = 1
          j2 = i-1
       CASE ( gt ) ! j > i
          j1 = i+1
          j2 = n
       CASE default ! should never happen
          WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error ', j_range
          STOP 'Impossible error in energy_1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    r_cut_box     = r_cut / box
    r_cut_box_sq  = r_cut_box**2
    rm_cut_box    = rm_cut / box
    rm_cut_box_sq = rm_cut_box**2
    box_sq        = box**2

    sr2     = 1.0 / r_cut**2 ! in sigma=1 units
    sr6     = sr2**3
    sr12    = sr6**2
    pot_cut = sr12 - sr6 ! Potential at cutoff (without numerical factor 4)

    energy = potovr ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2 ! Loop over j-molecules

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )

       IF ( rij_sq < rm_cut_box_sq ) THEN ! Test within molecular cutoff

          ! Double loop over atoms on both molecules
          DO a = 1, na
             DO b = 1, na
                rab = rij + di(:,a) - d(:,b,j) ! Atom-atom vector, box=1 units
                rab_sq = SUM ( rab**2 )

                IF ( rab_sq < r_cut_box_sq ) THEN ! Test within potential cutoff 

                   rab_sq = rab_sq * box_sq ! now in sigma=1 units
                   sr2    = 1.0 / rab_sq    ! (sigma/rab)**2

                   IF ( sr2 > sr2_overlap ) THEN
                      energy%ovr = .TRUE. ! Overlap detected
                      RETURN              ! Return immediately
                   END IF

                   sr6        = sr2**3
                   sr12       = sr6**2
                   potab      = sr12 - sr6
                   energy%pot = energy%pot + potab - pot_cut ! accumulate shifted potential
                   virab      = potab + sr12
                   fab        = rab * virab / rab_sq
                   energy%vir = energy%vir + DOT_PRODUCT ( rij, fab ) ! molecular virial
                END IF ! End test within potential cutoff

             END DO
          END DO
          ! End double loop over atoms on both molecules

       END IF ! end test within molecular cutoff

    END DO ! end loop over j-molecules

    ! Include numerical factors
    energy%pot = 4.0 * energy%pot
    energy%vir = 24.0 * energy%vir
    energy%vir = energy%vir / 3.0

    energy%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION energy_1

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
    IF ( SIZE(db,dim=2) /= na ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array mismatch', SIZE(db,dim=2), na
       STOP 'Error in q_to_d'
    END IF

    a = q_to_a ( q ) ! compute rotation matrix from quaternions
    DO i = 1, na
       ds(:,i) = MATMUL ( db(:,i), a ) ! NB: equivalent to ds = at*db, at=transpose of a
    END DO

  END FUNCTION q_to_d

END MODULE mc_module
