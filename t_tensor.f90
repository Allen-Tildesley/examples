! t_tensor.f90
! Electrostatic interactions: T-tensors compared with angles
PROGRAM t_tensor

  ! The dipole moment of molecule 1 is aligned along the axial vector e1.
  ! The quadrupole tensor, theta1, is diagonal and traceless with
  ! theta1xx = -0.5*thetam1, theta1yy = -0.5*thetam1, and theta1zz = thetam1 in the molecule-fixed system.
  ! The vector r12=r1-r2 points from 2 to 1.

  ! TODO MPA tidy up this program

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit
  USE               utility_module,  ONLY : random_orientation_vector
  IMPLICIT NONE

  REAL, DIMENSION(3)      :: r1, r2
  REAL, DIMENSION(3)      :: e1, e2, econ
  REAL, DIMENSION(3)      :: r12, r12u
  REAL, DIMENSION(3,3)    :: tt2, theta1, theta2, emat
  REAL, DIMENSION(3,3,3)  :: tt3
  REAL                    :: mu1, mu2, thetam1, thetam2
  REAL                    :: modr12, r12sq, r12cu, r12fo
  REAL                    :: vddt, vdde, vdqt, vdqe, vqdt, vqde, c1, c2, c12
  REAL                    :: f12tx, f12ty, f12tz, f12ex, f12ey, f12ez

  WRITE ( unit=output_unit, fmt='(a)' ) 'T-tensor '
  WRITE ( unit=output_unit, fmt='(a)' ) 'Calculation of electrostatic interactions beteeen linear molecules' 
  WRITE ( unit=output_unit, fmt='(a)' ) 'using T-tensors and Euler angles'

  ! Set up parameters for the pair interaction

  mu1      = 1.0                 ! dipole moment of molecule 1
  mu2      = 1.0                 ! dipole moment of molecule 2
  thetam1  = 1.0                 ! quadrupole moment of molecule 1
  thetam2  = 1.0                 ! quadrupole moment of molecule 2
  r1       = [ 2.0,  1.0,  1.0 ] ! set position of molecule 1
  r2       = [ 0.0,  0.0,  0.0 ] ! molecule is at the origin
  e1       = [ 0.0,  0.0,  1.0 ] ! axial vector of molecule 1
  e2       = [ 0.0,  1.0,  0.0 ] ! axial vector of molecule 2

  r12    = r1 - r2
  r12sq  = DOT_PRODUCT(r12, r12)
  modr12 = SQRT(r12sq)
  r12u   = r12 / modr12
  r12cu  = r12sq * modr12
  r12fo  = r12sq * r12sq

  ! Choose orientations at random

  CALL random_orientation_vector ( e1 )
  CALL random_orientation_vector ( e2 )

  WRITE ( unit=output_unit, fmt='(a,3f10.6)' ) 'Position of molecule    1', r1
  WRITE ( unit=output_unit, fmt='(a,3f10.6)' ) 'Orientation of molecule 1', e1
  WRITE ( unit=output_unit, fmt='(a,3f10.6)' ) 'Position of molecule    2', r2
  WRITE ( unit=output_unit, fmt='(a,3f10.6)' ) 'Orientation of molecule 2', e2

  ! calculate the elements of the quadrupole tensors in the space-fixed system

  theta1(1,1) = 1.5 * e1(1) * e1(1) - 0.5
  theta1(1,2) = 1.5 * e1(1) * e1(2)
  theta1(1,3) = 1.5 * e1(1) * e1(3)
  theta1(2,1) = 1.5 * e1(2) * e1(1)
  theta1(2,2) = 1.5 * e1(2) * e1(2) - 0.5
  theta1(2,3) = 1.5 * e1(2) * e1(3)
  theta1(3,1) = 1.5 * e1(3) * e1(1)
  theta1(3,2) = 1.5 * e1(3) * e1(2)
  theta1(3,3) = 1.5 * e1(3) * e1(3) - 0.5

  theta2(1,1) = 1.5 * e2(1) * e2(1) - 0.5
  theta2(1,2) = 1.5 * e2(1) * e2(2)
  theta2(1,3) = 1.5 * e2(1) * e2(3)
  theta2(2,1) = 1.5 * e2(2) * e2(1)
  theta2(2,2) = 1.5 * e2(2) * e2(2) - 0.5
  theta2(2,3) = 1.5 * e2(2) * e2(3)
  theta2(3,1) = 1.5 * e2(3) * e2(1)
  theta2(3,2) = 1.5 * e2(3) * e2(2)
  theta2(3,3) = 1.5 * e2(3) * e2(3) - 0.5

  ! calculate the second rank tensor, T2

  CALL t2_tensor(r12u, r12cu, tt2 )

  ! Calculate the dipole-dipole energy

  vddt    = - mu1 * mu2 * DOT_PRODUCT( e1, MATMUL(tt2, e2))

  c1      = DOT_PRODUCT(e1, r12u)
  c2      = DOT_PRODUCT(e2, r12u)
  c12     = DOT_PRODUCT(e1, e2)
  vdde    = mu1 * mu2 * ( c12 - 3.0 * c1 * c2 ) / r12cu

  WRITE ( unit=output_unit, fmt='(a,f10.6)' ) 'Dipole-dipole energy from T tensor     =', vddt
  WRITE ( unit=output_unit, fmt='(a,f10.6)' ) 'Dipole-dipole energy from Euler angles =', vdde

  ! Calculate the third rank tensor T3

  CALL t3_tensor(r12u, r12fo, tt3 )

  ! Contract T3 with the quadrupole tensor for molecule 2

  CALL contract_32(tt3, theta2, econ)

  ! Calculate the dipole-quadrupole energy

  vdqt = -(1.0/3.0) * mu1 * thetam2 * DOT_PRODUCT( e1, econ)

  vdqe =   1.5 * mu1 * thetam2 * ( c1 * (1.0 - 5.0 * c2 * c2) + 2.0 * c2 * c12 ) / r12fo

  WRITE ( unit=output_unit, fmt='(a,f10.6)' ) 'Dipole-quadrupole energy from T-tensors    =', vdqe
  WRITE ( unit=output_unit, fmt='(a,f10.6)' ) 'Dipole-quadrupole energy from Euler angles =', vdqt

  ! Calculate the quadrupole-dipole energy

  ! Contract T3 with the quadrupole tensor for molecule 1

  CALL contract_32(tt3, theta1, econ)

  vqdt  =  (1.0/3.0) * thetam1 * mu2 * DOT_PRODUCT( econ, e2 )
  vqde  = - 1.5 * thetam1 * mu2 * ( c2 * (1.0 - 5.0 * c1 * c1 ) + 2.0 * c1 * c12 ) / r12fo

  WRITE ( unit=output_unit, fmt='(a,f10.6)' ) 'Quadrupole-dipole energy from T-tensors    =', vqde
  WRITE ( unit=output_unit, fmt='(a,f10.6)' ) 'Quadrupole-dipole energy from Euler angles =', vqdt

  ! Calculate of the dipole-dipole force

  emat(1,1) = e1(1) * e2(1)
  emat(1,2) = e1(1) * e2(2)
  emat(1,3) = e1(1) * e2(3)
  emat(2,1) = e1(2) * e2(1)
  emat(2,2) = e1(2) * e2(2)
  emat(2,3) = e1(2) * e2(3)
  emat(3,1) = e1(3) * e2(1)
  emat(3,2) = e1(3) * e2(2)
  emat(3,3) = e1(3) * e2(3)

  ! Contract T3 with the dyadic of the axial vectors

  CALL contract_32(tt3, emat, econ)

  f12tx = - mu1 * mu2 * econ(1)
  f12ty = - mu1 * mu2 * econ(2)
  f12tz = - mu1 * mu2 * econ(3)

  f12ex = (3 * mu1 * mu2 / r12fo ) * ((c12 - 5.0 *c1 * c2) * r12u(1) + c2 * e1(1) + c1 * e2(1) )
  f12ey = (3 * mu1 * mu2 / r12fo ) * ((c12 - 5.0 *c1 * c2) * r12u(2) + c2 * e1(2) + c1 * e2(2) )
  f12ez = (3 * mu1 * mu2 / r12fo ) * ((c12 - 5.0 *c1 * c2) * r12u(3) + c2 * e1(3) + c1 * e2(3) )

  WRITE ( unit=output_unit, fmt='(a,3f10.6)' ) 'Dipole-dipole force from T-tensors    =', f12tx, f12ty, f12tz
  WRITE ( unit=output_unit, fmt='(a,3f10.6)' ) 'Dipole-dipole force from Euler angles =', f12ex, f12ey, f12ez


CONTAINS
  
  SUBROUTINE t2_tensor(r, r3, t2 )

    IMPLICIT NONE

    REAL, DIMENSION(:),   INTENT(in)  :: r  ! unit vector from j to i
    REAL,                 INTENT(in)  :: r3 ! third power of the modulus of rij
    REAL, DIMENSION(:,:), INTENT(out) :: t2 ! second-rank interaction tensor

    REAL :: txx, tyy, tzz
    REAL :: txy, txz, tyz

    ! calculate the 6 independent second rank elements

    txx = (3.0 * r(1) * r(1) - 1.0 ) / r3
    tyy = (3.0 * r(2) * r(2) - 1.0 ) / r3
    tzz = (3.0 * r(3) * r(3) - 1.0 ) / r3
    txy =  3.0 * r(1) * r(2)         / r3
    txz =  3.0 * r(1) * r(3)         / r3
    tyz =  3.0 * r(2) * r(3)         / r3

    ! assign the 9 elements of the second rank tensor

    t2(1, 1) = txx
    t2(2, 2) = tyy
    t2(3, 3) = tzz
    t2(1, 2) = txy
    t2(2, 1) = txy
    t2(1, 3) = txz
    t2(3, 1) = txz
    t2(2, 3) = tyz
    t2(3, 2) = tyz

  END SUBROUTINE t2_tensor

  SUBROUTINE t3_tensor(r, r4, t3 )

    ! calculates the components of a rank 3 interaction tensor

    IMPLICIT NONE

    REAL, DIMENSION(:),     INTENT(in)  :: r  ! unit vector from 2 to 1
    REAL,                   INTENT(in)  :: r4 ! fourth power of the modulus of r12
    REAL, DIMENSION(:,:,:), INTENT(out) :: t3 ! third-rank interaction tensor (note positive sign)


    REAL :: txxx, tyyy, tzzz, txyz
    REAL :: txxz, txxy, txyy, txzz, tyyz, tyzz

    ! calculate the unique third rank elements

    txxx = (15.0 * r(1) * r(1) * r(1) - 9.0 * r(1) ) / r4
    tyyy = (15.0 * r(2) * r(2) * r(2) - 9.0 * r(2) ) / r4
    tzzz = (15.0 * r(3) * r(3) * r(3) - 9.0 * r(3) ) / r4
    txyz =  15.0 * r(1) * r(2) * r(3)                / r4
    txxz =   3.0 * (5.0 * r(1) * r(1) * r(3) - r(3)) / r4
    txxy =   3.0 * (5.0 * r(1) * r(1) * r(2) - r(2)) / r4
    txyy =   3.0 * (5.0 * r(1) * r(2) * r(2) - r(1)) / r4
    txzz =   3.0 * (5.0 * r(1) * r(3) * r(3) - r(1)) / r4
    tyyz =   3.0 * (5.0 * r(2) * r(2) * r(3) - r(3)) / r4
    tyzz =   3.0 * (5.0 * r(2) * r(3) * r(3) - r(2)) / r4

    ! assign the 27 elements of the third rank tensor

    t3(1, 1, 1) = txxx
    t3(2, 2, 2) = tyyy
    t3(3, 3, 3) = tzzz
    t3(3, 2, 1) = txyz
    t3(2, 3, 1) = txyz
    t3(3, 1, 2) = txyz
    t3(1, 3, 2) = txyz
    t3(2, 1, 3) = txyz
    t3(1, 2, 3) = txyz
    t3(1, 1, 3) = txxz
    t3(1, 3, 1) = txxz
    t3(3, 1, 1) = txxz
    t3(1, 1, 2) = txxy
    t3(1, 2, 1) = txxy
    t3(2, 1, 1) = txxy
    t3(1, 2, 2) = txyy
    t3(2, 2, 1) = txyy
    t3(2, 1, 2) = txyy
    t3(1, 3, 3) = txzz
    t3(3, 3, 1) = txzz
    t3(3, 1, 3) = txzz
    t3(2, 2, 3) = tyyz
    t3(2, 3, 2) = tyyz
    t3(3, 2, 2) = tyyz
    t3(2, 3, 3) = tyzz
    t3(3, 2, 3) = tyzz
    t3(3, 3, 2) = tyzz

  END SUBROUTINE t3_tensor

  SUBROUTINE contract_32(mat3, mat2, mat1)

    ! Contracts a rank 3 and a rank 2 tensor to a rank 1 vector

    IMPLICIT NONE

    REAL, DIMENSION(:,:,:), INTENT(in)  :: mat3       ! third-rank  tensor
    REAL, DIMENSION(:,:),   INTENT(in)  :: mat2       ! second-rank tensor
    REAL, DIMENSION(:),     INTENT(out) :: mat1       ! first-rank contraction

    INTEGER :: i, j, k
    REAL    :: sum

    ! loop over first component of mat1

    DO i = 1, 3

       sum = 0.0

       DO j = 1, 3

          DO k = 1, 3

             sum = sum + mat3(i, j, k) * mat2(j, k)

          ENDDO

       ENDDO

       mat1(i) = sum

    ENDDO

  END SUBROUTINE contract_32

END PROGRAM t_tensor
