! t_tensor.f90
! Electrostatic interactions: T-tensors compared with angles
PROGRAM t_tensor

  ! The dipole moment of molecule 1 is aligned along the axial vector e1
  ! The quadrupole tensor, quad1, is diagonal and traceless with
  ! quad1_xx = -0.5*quad1_mag, quad1_yy = -0.5*quad1_mag, and quad1_zz = quad1_mag in the molecule-fixed system.
  ! Similarly for molecule 2
  ! The vector r12 = r1-r2 points from 2 to 1.

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit

  USE maths_module, ONLY : init_random_seed, random_vector, outer_product

  IMPLICIT NONE

  REAL, DIMENSION(3)      :: r12, r12_hat, e1, e2, econ, mu1, mu2, f12t, f12e
  REAL, DIMENSION(3,3)    :: tt2, quad1, quad2, emat
  REAL, DIMENSION(3,3,3)  :: tt3
  REAL                    :: r12_mag, r12_sq, r12_cu, r12_fo
  REAL                    :: vddt, vdde, vdqt, vdqe, vqdt, vqde, c1, c2, c12
  INTEGER                 :: i

  REAL, PARAMETER      :: d_min = 0.5, d_max = 1.5 ! desired range of separations
  REAL, PARAMETER      :: mu1_mag    = 1.0         ! dipole moment of molecule 1
  REAL, PARAMETER      :: mu2_mag    = 1.0         ! dipole moment of molecule 2
  REAL, PARAMETER      :: quad1_mag  = 1.0         ! quadrupole moment of molecule 1
  REAL, PARAMETER      :: quad2_mag  = 1.0         ! quadrupole moment of molecule 2

  WRITE ( unit=output_unit, fmt='(a)' ) 'T-tensor'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Calculation of electrostatic interactions between linear molecules' 
  WRITE ( unit=output_unit, fmt='(a)' ) 'using T-tensors and Euler angles'

  ! Initialize random number generator                          
  CALL init_random_seed

  ! Choose orientations at random
  e1 = random_vector ( )
  e2 = random_vector ( )

  ! Place atom 2 at origin and atom 1 in a random direction within desired distance range
  r12_hat = random_vector ( ) ! unit vector
  CALL RANDOM_NUMBER ( r12_mag )
  r12_mag = d_min + (d_max-d_min)*r12_mag ! magnitude of r12
  r12     = r12_hat * r12_mag             ! within desired range of origin
  r12_sq  = r12_mag**2                    ! squared distance
  r12_cu  = r12_sq * r12_mag              ! cubed distance
  r12_fo  = r12_sq * r12_sq               ! fourth power of distance

  c1     = DOT_PRODUCT ( e1, r12_hat )    ! cosine of angle between e1 and r12
  c2     = DOT_PRODUCT ( e2, r12_hat )    ! cosine of angle between e2 and r12
  c12    = DOT_PRODUCT ( e1, e2   )       ! cosine of angle between e1 and e2

  WRITE ( unit=output_unit, fmt='(a,t30,3f10.6)' ) 'Displacement r12 = ', r12
  WRITE ( unit=output_unit, fmt='(a,t30,3f10.6)' ) 'Orientation  e1  = ', e1
  WRITE ( unit=output_unit, fmt='(a,t30,3f10.6)' ) 'Orientation  e2  = ', e2

  ! dipole vectors in space-fixed frame
  mu1 = mu1_mag * e1                             
  mu2 = mu2_mag * e2

  ! quadrupole tensors in space-fixed frame (traceless)
  quad1 = 1.5 * outer_product ( e1, e1 )        
  FORALL (i=1:3) quad1(i,i) = quad1(i,i) - 0.5
  quad1 = quad1_mag * quad1
  quad2 = 1.5 * outer_product ( e2, e2 )
  FORALL (i=1:3) quad2(i,i) = quad2(i,i) - 0.5
  quad2 = quad2_mag * quad2

  tt2 = t2_tensor ( r12_hat, r12_cu ) ! Calculate the second rank tensor, T2
  tt3 = t3_tensor ( r12_hat, r12_fo ) ! Calculate the third rank tensor, T3

  ! Headings
  WRITE ( unit=output_unit, fmt='(/,t30,a30,t70,a30)' ) '.....Result from T tensor', '.....Result from Euler angles' 

  ! Calculate the dipole-dipole energy
  vddt = -DOT_PRODUCT ( mu1, MATMUL ( tt2, mu2 ) ) ! Contract both dipoles with T2
  vdde = mu1_mag * mu2_mag * ( c12 - 3.0 * c1 * c2 ) / r12_cu
  WRITE ( unit=output_unit, fmt='(a,t30,f10.6,t70,f10.6)' ) 'Dipole-dipole energy =', vddt, vdde
  
  ! Calculate the dipole-quadrupole energy
  econ = contract_32 ( tt3, quad2 ) ! Contract T3 with the quadrupole tensor for molecule 2
  vdqt = -(1.0/3.0) * DOT_PRODUCT ( mu1, econ ) ! and contract with dipole for molecule 1
  vdqe = 1.5 * mu1_mag * quad2_mag * ( c1 * (1.0 - 5.0 * c2 * c2) + 2.0 * c2 * c12 ) / r12_fo
  WRITE ( unit=output_unit, fmt='(a,t30,f10.6,t70,f10.6)' ) 'Dipole-quadrupole energy =', vdqe, vdqt

  ! Calculate the quadrupole-dipole energy
  econ = contract_32 ( tt3, quad1 )  ! Contract T3 with the quadrupole tensor for molecule 1
  vqdt  =  (1.0/3.0) * DOT_PRODUCT ( econ, mu2 ) ! and contract with dipole for molecule 2
  vqde  = - 1.5 * quad1_mag * mu2_mag * ( c2 * (1.0 - 5.0 * c1 * c1 ) + 2.0 * c1 * c12 ) / r12_fo
  WRITE ( unit=output_unit, fmt='(a,t30,f10.6,t70,f10.6)' ) 'Quadrupole-dipole energy =', vqde, vqdt

  ! Calculate the dipole-dipole force
  emat = outer_product ( mu1, mu2 ) ! dyadic of the dipole vectors
  econ = contract_32 ( tt3, emat )  ! Contract with T3
  f12t = - econ
  f12e = (3.0 * mu1_mag * mu2_mag / r12_fo ) * ((c12 - 5.0 *c1 * c2) * r12_hat + c2 * e1 + c1 * e2 )
  WRITE ( unit=output_unit, fmt='(a,t30,3f10.6,t70,3f10.6)' ) 'Dipole-dipole force  =', f12t, f12e

CONTAINS

  FUNCTION t2_tensor ( r, r3 ) RESULT ( t2 ) ! second-rank 3x3 interaction tensor

    IMPLICIT NONE

    REAL, DIMENSION(3),   INTENT(in)  :: r  ! unit vector from 2 to 1
    REAL,                 INTENT(in)  :: r3 ! third power of the modulus of r12
    REAL, DIMENSION(3,3)              :: t2 ! function result

    INTEGER :: i

    t2 = 3.0 * outer_product ( r, r )
    FORALL (i=1:3) t2(i,i) = t2(i,i) - 1.0 ! make traceless

    t2 = t2 / r3 ! Scale by third power of distance

  END FUNCTION t2_tensor

  FUNCTION t3_tensor ( r, r4 ) RESULT ( t3 ) ! third-rank 3x3x3 interaction tensor (note positive sign)

    IMPLICIT NONE

    REAL, DIMENSION(3),     INTENT(in)  :: r  ! unit vector from 2 to 1
    REAL,                   INTENT(in)  :: r4 ! fourth power of the modulus of r12
    REAL, DIMENSION(3,3,3)              :: t3 ! function result

    INTEGER :: i, j

    t3 = 15.0 * outer_product ( r, r, r )

    DO i = 1, 3
       t3(i,i,i) = t3(i,i,i) - 9.0 * r(i) ! correction for all indices the same

       DO j = 1, 3
          IF ( j == i ) CYCLE
          t3(i,i,j) = t3(i,i,j) - 3.0 * r(j) ! correction for two indices the same
          t3(i,j,i) = t3(i,j,i) - 3.0 * r(j) ! correction for two indices the same
          t3(j,i,i) = t3(j,i,i) - 3.0 * r(j) ! correction for two indices the same
       END DO
    END DO

    t3 = t3 / r4 ! Scale by fourth power of distance

  END FUNCTION t3_tensor

  FUNCTION contract_32 ( mat3, mat2 ) RESULT ( mat1 ) ! Contracts a rank-3 and a rank-2 tensor to a rank-1 vector

    IMPLICIT NONE

    REAL, DIMENSION(3,3,3), INTENT(in)  :: mat3       ! third-rank  tensor
    REAL, DIMENSION(3,3),   INTENT(in)  :: mat2       ! second-rank tensor
    REAL, DIMENSION(3)                  :: mat1       ! function result

    INTEGER :: i

    DO i = 1, 3
       mat1(i) = SUM ( mat3(i,:,:) * mat2(:,:) ) 
    END DO

  END FUNCTION contract_32

END PROGRAM t_tensor
