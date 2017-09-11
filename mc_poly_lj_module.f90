! mc_poly_lj_module.f90
! Routines for MC simulation, polyatomic molecule, LJ atoms
MODULE mc_module

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
  PUBLIC :: potential_1, potential

  ! Public data
  INTEGER,                                PUBLIC :: n ! number of molecules
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: r ! centre of mass positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: e ! quaternions (0:3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: d ! bond vectors (0:3,na,n)

  ! Bond vectors in body-fixed frame (na and db are public)
  ! Isosceles triangle, 3 sites, with unit bond length and bond angle alpha, which we set to 75 degrees here
  REAL,                  PARAMETER         :: pi = 4.0*ATAN(1.0)
  REAL,                  PARAMETER         :: alpha = 75.0 * pi / 180.0, alpha2 = alpha / 2.0
  INTEGER,               PARAMETER, PUBLIC :: na = 3
  REAL, DIMENSION(3,na), PARAMETER, PUBLIC :: db = RESHAPE ( [ &
       & -SIN(alpha2), 0.0,    -COS(alpha2)/3.0, & 
       &  0.0,         0.0, 2.0*COS(alpha2)/3.0, &
       &  SIN(alpha2), 0.0,    -COS(alpha2)/3.0 ], [3,na] )

  ! Cutoff distance and force-shift parameters (all private) chosen as per the reference:
  ! S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)
  REAL, PARAMETER :: r_cut   = 2.612 ! in sigma=1 units, where r_cut = 1.2616 nm, sigma = 0.483 nm
  REAL, PARAMETER :: sr_cut  = 1.0/r_cut, sr_cut6 = sr_cut**6, sr_cut12 = sr_cut6**2
  REAL, PARAMETER :: lambda1 = 4.0*(7.0*sr_cut6-13.0*sr_cut12)
  REAL, PARAMETER :: lambda2 = -24.0*(sr_cut6-2.0*sr_cut12)*sr_cut
  
  INTEGER, PARAMETER :: lt = -1, gt = 1 ! Options for j-range

  ! Public derived type
  TYPE, PUBLIC :: potential_type   ! A composite variable for interactions comprising
     REAL    :: pot ! the potential energy cut and shifted to 0 at r_cut, and
     REAL    :: vir ! the virial and
     LOGICAL :: ovr ! a flag indicating overlap (i.e. pot too high to use)
   CONTAINS
     PROCEDURE :: add_potential_type
     PROCEDURE :: subtract_potential_type
     GENERIC   :: OPERATOR(+) => add_potential_type
     GENERIC   :: OPERATOR(-) => subtract_potential_type
  END TYPE potential_type

CONTAINS

  FUNCTION add_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the sum of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot  +   b%pot
    c%vir = a%vir  +   b%vir
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  FUNCTION subtract_potential_type ( a, b ) RESULT (c)
    IMPLICIT NONE
    TYPE(potential_type)              :: c    ! Result is the difference of
    CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
    c%pot = a%pot  -   b%pot
    c%vir = a%vir  -   b%vir
    c%ovr = a%ovr .OR. b%ovr ! This is meaningless, but inconsequential
  END FUNCTION subtract_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    INTEGER :: i
    REAL    :: diameter

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut-and-force-shifted'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'    
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'  

    WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Number of atoms per molecule', na
    DO i = 1, na ! Loop over atoms
       WRITE ( unit=output_unit, fmt='(a,i1,t40,3f15.6)' ) 'Body-fixed atom vector ', i, db(:,i)
    END DO ! End loop over atoms

    diameter = 2.0 * SQRT ( MAXVAL ( SUM(db**2,dim=1) ) )
    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Molecular diameter', diameter
    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'r_cut', r_cut
    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Force-shift lambda1', lambda1
    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Force-shift lambda2', lambda2

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box )
    IMPLICIT NONE
    REAL, INTENT(in) :: box   ! simulation box length

    REAL :: rm_cut_box, diameter

    ALLOCATE ( r(3,n), e(0:3,n), d(3,na,n) )

    diameter   = 2.0 * SQRT ( MAXVAL ( SUM(db**2,dim=1) ) )
    rm_cut_box = ( r_cut+diameter ) / box
    IF ( rm_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)') 'rm_cut/box too large ', rm_cut_box
       STOP 'Error in allocate_arrays'
    END IF

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, e, d )

  END SUBROUTINE deallocate_arrays

  FUNCTION potential ( box ) RESULT ( total )
    IMPLICIT NONE
    TYPE(potential_type) :: total ! Returns a composite of pot, vir etc
    REAL, INTENT(in)     :: box   ! Simulation box length

    ! total%pot is the nonbonded cut-and-shifted potential energy for whole system
    ! total%vir is the corresponding virial for whole system
    ! total%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this flag is .true., the values of total%pot etc should not be used
    ! Actual calculation is performed by function potential_1

    TYPE(potential_type) :: partial ! Molecular contribution to total
    INTEGER              :: i

    IF ( ANY ( SHAPE(r)  /= [3,n] ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,3i15)' ) 'Array bounds error for r', n, SHAPE(r)
       STOP 'Error in potential'
    END IF

    IF ( ANY ( SHAPE(d)  /= [3,na,n] ) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,5i15)' ) 'Array bounds error for d', na, n, SHAPE(d)
       STOP 'Error in potential'
    END IF

    total = potential_type ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO i = 1, n - 1

       partial = potential_1 ( r(:,i), d(:,:,i), i, box, gt )

       IF ( partial%ovr ) THEN
          total%ovr = .TRUE. ! Overlap detected
          RETURN             ! Return immediately
       END IF

       total = total + partial

    END DO

    total%ovr = .FALSE. ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential

  FUNCTION potential_1 ( ri, di, i, box, j_range ) RESULT ( partial )
    IMPLICIT NONE
    TYPE(potential_type)                 :: partial ! Returns a composite of pot, vir etc for given molecule
    REAL,    DIMENSION(3),    INTENT(in) :: ri      ! Coordinates of molecule of interest
    REAL,    DIMENSION(3,na), INTENT(in) :: di      ! Bond vectors of molecule of interest
    INTEGER,                  INTENT(in) :: i       ! Index of molecule of interest
    REAL,                     INTENT(in) :: box     ! Simulation box length
    INTEGER, OPTIONAL,        INTENT(in) :: j_range ! Optional partner index range

    ! partial%pot is the nonbonded cut-and-shifted potential energy of molecule ri,ei with a set of other molecules
    ! partial%vir is the corresponding virial of molecule ri,ei
    ! partial%ovr is a flag indicating overlap (potential too high) to avoid overflow
    ! If this is .true., the values of partial%pot etc should not be used
    ! The coordinates in ri and ei are not necessarily identical with those in r(:,i) and e(:,i)
    ! The optional argument j_range restricts partner indices to j>i, or j<i

    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1
    ! Note that this is the force-shifted LJ potential with a linear smoothing term
    ! S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)

    INTEGER               :: j, j1, j2, a, b
    REAL                  :: diameter, rm_cut_box, rm_cut_box_sq, r_cut_sq
    REAL                  :: sr2, sr6, sr12, rij_sq, rab_sq, virab, rmag
    REAL, DIMENSION(3)    :: rij, rab, fab
    REAL, PARAMETER       :: sr2_ovr = 1.77 ! overlap threshold (pot > 100)
    TYPE(potential_type)  :: pair

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

    diameter      = 2.0 * SQRT ( MAXVAL ( SUM(db**2,dim=1) ) )
    rm_cut_box    = ( r_cut + diameter ) / box ! Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box**2              ! squared
    r_cut_sq      = r_cut**2                   ! Potential cutoff squared in sigma=1 units

    partial = potential_type ( pot=0.0, vir=0.0, ovr=.FALSE. ) ! Initialize

    DO j = j1, j2 ! Loop over selected range of partner molecules

       IF ( i == j ) CYCLE ! Skip self

       rij(:) = ri(:) - r(:,j)            ! Centre-centre separation vector
       rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )            ! Squared centre-centre separation in box=1 units

       IF ( rij_sq < rm_cut_box_sq ) THEN ! Test within molecular cutoff

          rij = rij * box ! Now in sigma = 1 units

          ! Double loop over atoms on both molecules
          DO a = 1, na      
             DO b = 1, na

                rab    = rij + di(:,a) - d(:,b,j) ! Atom-atom vector, sigma=1 units
                rab_sq = SUM ( rab**2 )           ! Squared atom-atom separation, sigma=1 units

                IF ( rab_sq < r_cut_sq ) THEN ! Test within potential cutoff 

                   sr2      = 1.0 / rab_sq  ! (sigma/rab)**2
                   pair%ovr = sr2 > sr2_ovr ! Overlap if too close

                   IF ( pair%ovr ) THEN
                      partial%ovr = .TRUE. ! Overlap detected
                      RETURN               ! Return immediately
                   END IF

                   rmag     = SQRT(rab_sq)
                   sr6      = sr2**3
                   sr12     = sr6**2
                   pair%pot = 4.0*(sr12-sr6) + lambda1 + lambda2*rmag ! LJ atom-atom pair potential (force-shifted)
                   virab    = 24.0*(2.0*sr12-sr6 ) - lambda2*rmag     ! LJ atom-atom pair virial
                   fab      = rab * virab * sr2                       ! LJ atom-atom pair force
                   pair%vir = DOT_PRODUCT ( rij, fab )                ! Contribution to molecular virial

                   partial = partial + pair

                END IF ! End test within potential cutoff

             END DO
          END DO
          ! End double loop over atoms on both molecules

       END IF ! End test within molecular cutoff

    END DO ! End loop over selected range of partner molecules

    ! Include numerical factors
    partial%vir = partial%vir / 3.0 ! Divide virial by 3
    partial%ovr = .FALSE.           ! No overlaps detected (redundant, but for clarity)

  END FUNCTION potential_1

END MODULE mc_module
