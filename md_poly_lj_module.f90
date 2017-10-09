! md_poly_lj_module.f90
! Force routine for MD simulation, polyatomic molecule, LJ atoms
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
  PUBLIC :: force

  ! Public data
  INTEGER,                                PUBLIC :: n   ! number of molecules
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: r   ! centre of mass positions (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: e   ! quaternions (0:3,n)
  REAL,    DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: d   ! bond vectors (0:3,na,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: v   ! centre of mass velocities (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: ell ! angular momenta (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: f   ! centre of mass forces (3,n)
  REAL,    DIMENSION(:,:),   ALLOCATABLE, PUBLIC :: tau ! torques (3,n)

  ! Bond vectors in body-fixed frame (na and db are public)
  ! Isosceles triangle, 3 sites, with unit bond length and bond angle alpha, which we set to 75 degrees here
  REAL,                  PARAMETER         :: pi = 4.0*ATAN(1.0)
  REAL,                  PARAMETER         :: alpha = 75.0 * pi / 180.0, alpha2 = alpha / 2.0
  INTEGER,               PARAMETER, PUBLIC :: na = 3
  REAL, DIMENSION(3,na), PARAMETER, PUBLIC :: db = RESHAPE ( [ &
       & -SIN(alpha2), 0.0,    -COS(alpha2)/3.0, & 
       &  0.0,         0.0, 2.0*COS(alpha2)/3.0, &
       &  SIN(alpha2), 0.0,    -COS(alpha2)/3.0 ], [3,na] )

  ! Atomic masses and (public) moments of inertia
  REAL, DIMENSION(na), PARAMETER        :: m = [ 1.0/3.0, 1.0/3.0, 1.0/3.0 ] ! Masses add up to 1.0
  REAL, DIMENSION(3),            PUBLIC :: inertia
  
  ! Cutoff distance and force-shift parameters (all private) chosen as per the reference:
  ! S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)
  REAL, PARAMETER :: r_cut   = 2.612 ! in sigma=1 units, where r_cut = 1.2616 nm, sigma = 0.483 nm
  REAL, PARAMETER :: sr_cut  = 1.0/r_cut, sr_cut6 = sr_cut**6, sr_cut12 = sr_cut6**2
  REAL, PARAMETER :: lambda1 = 4.0*(7.0*sr_cut6-13.0*sr_cut12)
  REAL, PARAMETER :: lambda2 = -24.0*(sr_cut6-2.0*sr_cut12)*sr_cut

  ! Public derived type
  TYPE, PUBLIC :: potential_type   ! A composite variable for interactions comprising
     REAL    :: pot ! the potential energy
     REAL    :: vir ! the virial and
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
    c%pot = a%pot  +   b%pot
    c%vir = a%vir  +   b%vir
    c%ovr = a%ovr .OR. b%ovr
  END FUNCTION add_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    INTEGER               :: i
    REAL                  :: diameter
    REAL, DIMENSION(3)    :: com
    REAL, PARAMETER       :: tol = 1.0e-9

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

    ! The following section sets the diagonal moments of inertia "realistically"
    ! based on the values of atomic masses and bond vectors (above), with some checking.
    ! However, there is nothing to stop the user replacing this section with a statement setting
    ! the values of inertia(1:3). The masses m are not passed back to the calling program.
    ! It might be advantageous, for instance, to artificially increase the values in inertia.
    
    ! Ensure that the db bonds, xyz molecular axes, and masses are chosen such that
    ! the total mass is 1 and the centre-of-mass is at the origin
    IF ( ABS ( SUM(m) - 1.0 ) > tol ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)') 'Molecular mass is not 1.0 ', SUM(m)
       STOP 'Error in introduction'
    END IF
    com = SUM ( SPREAD(m,dim=1,ncopies=3)*db, dim = 2 )
    IF ( ANY ( ABS(com) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.6)') 'Molecular centre-of-mass error', com
       STOP 'Error in introduction'
    END IF

    ! Ensure that the db bonds, xyz molecular axes, and masses are chosen such that
    ! the off-diagonal elements of the inertia tensor are zero
    inertia(1) = -SUM ( m*db(2,:)*db(3,:) )
    inertia(2) = -SUM ( m*db(3,:)*db(1,:) )
    inertia(3) = -SUM ( m*db(1,:)*db(2,:) )
    IF ( ANY ( ABS(inertia) > tol ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,3f15.6)') 'Molecular off-diagonal inertia error', inertia
       STOP 'Error in introduction'
    END IF

    ! Calculate the diagonal elements of the inertia tensor
    inertia(1) = SUM ( m*db(2,:)**2 ) + SUM ( m*db(3,:)**2 ) 
    inertia(2) = SUM ( m*db(3,:)**2 ) + SUM ( m*db(1,:)**2 ) 
    inertia(3) = SUM ( m*db(1,:)**2 ) + SUM ( m*db(2,:)**2 )
    WRITE ( unit=output_unit, fmt='(a,t40,3f15.6)' ) 'Inertia Ixx, Iyy, Izz', inertia

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
    ALLOCATE ( v(3,n), ell(3,n), f(3,n), tau(3,n) )

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
    DEALLOCATE ( v, ell, f, tau )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE force ( box, total )
    USE maths_module, ONLY : cross_product
    IMPLICIT NONE
    REAL,                 INTENT(in)  :: box   ! Simulation box length
    TYPE(potential_type), INTENT(out) :: total ! Composite of pot, vir, etc

    ! total%pot is the nonbonded cut-and-shifted potential energy for whole system
    ! total%vir is the corresponding virial
    ! total%ovr is a warning flag that there is an overlap
    ! This routine also calculates forces and torques
    ! and stores them in the arrays f, tau

    ! It is assumed that r has been divided by box
    ! The d array of space-fixed bond vectors for each molecule should be computed already
    ! Results are in LJ units where sigma = 1, epsilon = 1
    ! Note that this is the force-shifted LJ potential with a linear smoothing term
    ! S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)

    INTEGER              :: i, j, a, b
    REAL                 :: diameter, rm_cut_box, rm_cut_box_sq, r_cut_sq
    REAL                 :: sr2, sr6, sr12, rij_sq, rab_sq, virab, rmag
    REAL, DIMENSION(3)   :: rij, rab, fab
    REAL, PARAMETER      :: sr2_ovr = 1.77 ! overlap threshold (pot > 100)
    TYPE(potential_type) :: pair

    diameter      = 2.0 * SQRT ( MAXVAL ( SUM(db**2,dim=1) ) )
    rm_cut_box    = ( r_cut + diameter ) / box ! Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box**2              ! squared
    r_cut_sq      = r_cut**2                   ! Potential cutoff squared in sigma=1 units

    ! Initialize
    f     = 0.0
    tau   = 0.0
    total = potential_type ( pot=0.0, vir=0.0, ovr=.FALSE. )

    DO i = 1, n - 1 ! Begin outer loop over molecules

       DO j = i + 1, n ! Begin inner loop over molecules

          rij(:) = r(:,i) - r(:,j)           ! Centre-centre separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! Periodic boundaries in box=1 units
          rij_sq = SUM ( rij**2 )            ! Squared centre-centre separation in box=1 units

          IF ( rij_sq < rm_cut_box_sq ) THEN ! Test within molecular cutoff

             rij = rij * box ! Now in sigma = 1 units

             ! Double loop over atoms on both molecules
             DO a = 1, na      
                DO b = 1, na

                   rab    = rij + d(:,a,i) - d(:,b,j) ! Atom-atom vector, sigma=1 units
                   rab_sq = SUM ( rab**2 )            ! Squared atom-atom separation, sigma=1 units

                   IF ( rab_sq < r_cut_sq ) THEN ! Test within potential cutoff 

                      sr2      = 1.0 / rab_sq  ! (sigma/rab)**2
                      pair%ovr = sr2 > sr2_ovr ! Overlap if too close

                      rmag     = SQRT(rab_sq)
                      sr6      = sr2**3
                      sr12     = sr6**2
                      pair%pot = 4.0*(sr12-sr6) + lambda1 + lambda2*rmag ! LJ atom-atom pair potential (force-shifted)
                      virab    = 24.0*(2.0*sr12-sr6 ) - lambda2*rmag     ! LJ atom-atom pair virial
                      fab      = rab * virab * sr2                       ! LJ atom-atom pair force
                      pair%vir = DOT_PRODUCT ( rij, fab )                ! Contribution to molecular virial

                      total    = total + pair
                      f(:,i)   = f(:,i) + fab
                      f(:,j)   = f(:,j) - fab
                      tau(:,i) = tau(:,i) + cross_product ( d(:,a,i), fab )
                      tau(:,j) = tau(:,j) - cross_product ( d(:,b,j), fab )

                   END IF ! End test within potential cutoff

                END DO
             END DO
             ! End double loop over atoms on both molecules

          END IF ! End test within molecular cutoff

       END DO ! End inner loop over molecules
    END DO ! End outer loop over molecules

    ! Include numerical factors
    total%vir = total%vir / 3.0 ! Divide virial by 3

  END SUBROUTINE force

END MODULE md_module
