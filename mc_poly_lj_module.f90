! mc_poly_lj_module.f90
! Routines for MC simulation, polyatomic molecule, LJ atoms
MODULE mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: potential_1, potential

  ! Public data
  INTEGER,                                PUBLIC :: n ! number of molecules
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: r ! centre of mass positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: e ! quaternions (0:3,n)

  ! Private data
  INTEGER, PARAMETER                  :: na = 3 ! Number of atoms per molecule
  REAL,    PARAMETER, DIMENSION(3,na) :: db = RESHAPE ( [ &
       &  0.0, 0.0,  1.0/SQRT(3.0), & 
       & -0.5, 0.0, -0.5/SQRT(3.0), &
       &  0.5, 0.0, -0.5/SQRT(3.0) ], [3,na] ) ! Bond vectors in body-fixed frame

  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range
  REAL,    PARAMETER :: sigma = 1.0     ! Lennard-Jones diameter (unit of length)
  REAL,    PARAMETER :: epslj = 1.0     ! Lennard-Jones well depth (unit of energy)

  REAL,    PARAMETER :: diameter = 2.0 * SQRT ( MAXVAL ( SUM(db**2,dim=1) ) ) ! Molecular diameter

  ! Public derived type
  PUBLIC :: potential_type
  TYPE potential_type   ! A composite variable for interactions comprising
     REAL    :: pot     ! the potential energy and
     REAL    :: vir     ! the virial and
     LOGICAL :: overlap ! a flag indicating overlap (i.e. pot too high to use)
  END TYPE potential_type

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_potential_type
  END INTERFACE OPERATOR (+)

CONTAINS

  FUNCTION add_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)             :: c    ! Result is the sum of the two inputs
    TYPE(potential_type), INTENT(in) :: a, b
    c%pot     = a%pot      +   b%pot
    c%vir     = a%vir      +   b%vir
    c%overlap = a%overlap .OR. b%overlap
  END FUNCTION add_potential_type

  SUBROUTINE introduction ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    INTEGER :: i

    WRITE ( unit=output_unit, fmt='(a)'           ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)'           ) 'Cut-and-shifted, with no long-range corrections'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, sigma = ',   sigma    
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Well depth, epslj = ', epslj    

    WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Number of atoms per molecule', na
    DO i = 1, na ! Loop over atoms
       WRITE ( unit=output_unit, fmt='(a,i1,t40,3f15.5)' ) 'Body-fixed atom vector ', i, db(:,i)
    END DO ! End loop over atoms
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Molecular diameter ', diameter

  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box, r_cut )
    IMPLICIT NONE
    REAL, INTENT(in) :: box   ! simulation box length
    REAL, INTENT(in) :: r_cut ! potential cutoff distance

    REAL :: rm_cut_box

    ALLOCATE ( r(3,n), e(0:3,n) )

    rm_cut_box = ( r_cut+diameter ) / box
    IF ( rm_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)') 'rm_cut/box too large ', rm_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, e )

  END SUBROUTINE deallocate_arrays

  FUNCTION potential ( box, r_cut ) RESULT ( system )
    IMPLICIT NONE
    TYPE(potential_type) :: system ! Returns a composite of pot, vir and overlap
    REAL, INTENT(in)     :: box    ! Simulation box length
    REAL, INTENT(in)     :: r_cut  ! Potential cutoff distance

    ! system%pot is the nonbonded potential system for whole system
    ! system%vir is the corresponding virial for whole system
    ! system%overlap is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of system%pot, system%vir should not be used
    ! Actual calculation is performed by function potential_1

    TYPE(potential_type) :: molecule
    INTEGER              :: i

    IF ( SIZE(r,dim=2)  /= n ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in potential'
    END IF

    system = potential_type ( pot=0.0, vir=0.0, overlap=.FALSE. ) ! Initialize

    DO i = 1, n - 1

       molecule = potential_1 ( r(:,i), e(:,i), i, box, r_cut, gt )

       IF ( molecule%overlap ) THEN
          system%overlap = .TRUE. ! Overlap detected
          RETURN                  ! Return immediately
       END IF

       system = system + molecule

    END DO

    system%overlap = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential

  FUNCTION potential_1 ( ri, ei, i, box, r_cut, j_range ) RESULT ( molecule )
    USE maths_module, ONLY : q_to_a
    IMPLICIT NONE
    TYPE(potential_type)                :: molecule ! Returns a composite of pot, vir and overlap
    REAL,    DIMENSION(3),   INTENT(in) :: ri       ! Coordinates of molecule of interest
    REAL,    DIMENSION(0:3), INTENT(in) :: ei       ! Quaternion of molecule of interest
    INTEGER,                 INTENT(in) :: i        ! Index of molecule of interest
    REAL,                    INTENT(in) :: box      ! Simulation box length
    REAL,                    INTENT(in) :: r_cut    ! Potential cutoff distance
    INTEGER, OPTIONAL,       INTENT(in) :: j_range  ! Optional partner index range

    ! molecule%pot is the nonbonded potential energy of molecule ri,ei with a set of other molecules
    ! molecule%vir is the corresponding virial of molecule ri,ei
    ! molecule%overlap is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the value of molecule%pot should not be used
    ! The coordinates in ri and ei are not necessarily identical with those in r(:,i) and e(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1
    ! Note that this is the shifted LJ potential

    INTEGER               :: j, j1, j2, a, b
    REAL                  :: rm_cut_box, rm_cut_box_sq, r_cut_sq
    REAL                  :: sr2, sr6, sr12, rij_sq, rab_sq, potab, virab, pot_cut
    REAL, DIMENSION(3)    :: rij, rab, fab
    REAL, DIMENSION(3,na) :: di, dj
    REAL, DIMENSION(3,3)  :: ai, aj
    REAL, PARAMETER       :: sr2_overlap = 1.8 ! overlap threshold

    IF ( SIZE(r,dim=2)  /= n ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in potential_1'
    END IF
    IF ( SIZE(e,dim=2)  /= n ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for e', n, SIZE(e,dim=2)
       STOP 'Error in potential_1'
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
          STOP 'Impossible error in potential_1'
       END SELECT
    ELSE
       j1 = 1
       j2 = n
    END IF

    rm_cut_box    = ( r_cut + diameter ) / box ! Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box**2              ! squared
    r_cut_sq      = r_cut**2                   ! Potential cutoff squared in sigma=1 units

    sr2     = 1.0 / r_cut_sq ! in sigma=1 units
    sr6     = sr2**3
    sr12    = sr6**2
    pot_cut = sr12 - sr6 ! Potential at cutoff (without numerical factor 4)

    molecule = potential_type ( pot=0.0, vir=0.0, overlap=.FALSE. ) ! Initialize

    ! Compute all space-fixed atom vectors for molecule i
    ai = q_to_a ( ei ) ! Rotation matrix for i
    DO a = 1, na ! Loop over all atoms
       di(:,a) = MATMUL ( db(:,a), ai ) ! NB: equivalent to ai_T*db, ai_T=transpose of ai
    END DO ! End loop over all atoms

    DO j = j1, j2 ! Loop over selected range of partner molecules

       IF ( i == j ) CYCLE ! Skip self

       rij(:) = ri(:) - r(:,j)            ! Centre-centre separation vector
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )            ! Squared centre-centre separation in box=1 units

       IF ( rij_sq < rm_cut_box_sq ) THEN ! Test within molecular cutoff

          rij = rij * box ! now in sigma = 1 units

          ! Compute all space-fixed atom vectors for molecule j
          aj = q_to_a ( e(:,j) ) ! Rotation matrix for j
          DO a = 1, na ! Loop over all atoms
             dj(:,a) = MATMUL ( db(:,a), aj ) ! NB: equivalent to aj_T*db, aj_T=transpose of aj
          END DO ! End loop over all atoms

          ! Double loop over atoms on both molecules
          DO a = 1, na      
             DO b = 1, na

                rab    = rij + di(:,a) - dj(:,b) ! Atom-atom vector, sigma=1 units
                rab_sq = SUM ( rab**2 )          ! Squared atom-atom separation, sigma=1 units

                IF ( rab_sq < r_cut_sq ) THEN ! Test within potential cutoff 

                   sr2 = 1.0 / rab_sq    ! (sigma/rab)**2

                   IF ( sr2 > sr2_overlap ) THEN
                      molecule%overlap = .TRUE. ! Overlap detected
                      RETURN                    ! Return immediately
                   END IF

                   sr6          = sr2**3
                   sr12         = sr6**2
                   potab        = sr12 - sr6
                   molecule%pot = molecule%pot + potab - pot_cut          ! LJ potential (cut-and-shifted) 
                   virab        = potab + sr12                            ! atom-atom virial
                   fab          = rab * virab / rab_sq                    ! atom-atom force
                   molecule%vir = molecule%vir + DOT_PRODUCT ( rij, fab ) ! LJ molecular virial

                END IF ! End test within potential cutoff

             END DO
          END DO
          ! End double loop over atoms on both molecules

       END IF ! End test within molecular cutoff

    END DO ! End loop over selected range of partner molecules

    ! Include numerical factors
    molecule%pot     = molecule%pot * 4.0
    molecule%vir     = molecule%vir * 24.0 / 3.0
    molecule%overlap = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential_1

END MODULE mc_module
