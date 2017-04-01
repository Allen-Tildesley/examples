! t_tensor_module.f90
! Routines to define T-tensors and carry out index contractions
MODULE t_tensor_module

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <dominic.tildesley@epfl.ch> ("the authors"),                          !
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

  USE maths_module, ONLY : outer_product

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: t2_tensor, t3_tensor, t4_tensor, t5_tensor, contract, skew_contract

  ! The contraction routines always exhaust the indices of the second argument
  ! Also, the first argument is always completely symmetric with respect to its indices
  INTERFACE contract
     MODULE PROCEDURE contract_11_0 ! Contract rank-1 with rank-1 to give rank-0
     MODULE PROCEDURE contract_21_1 ! Contract rank-2 with rank-1 to give rank-1
     MODULE PROCEDURE contract_22_0 ! Contract rank-2 with rank-2 to give rank-0
     MODULE PROCEDURE contract_32_1 ! Contract rank-3 with rank-2 to give rank-1
     MODULE PROCEDURE contract_31_2 ! Contract rank-3 with rank-1 to give rank-2
     MODULE PROCEDURE contract_42_2 ! Contract rank-4 with rank-2 to give rank-2
     MODULE PROCEDURE contract_52_3 ! Contract rank-5 with rank-2 to give rank-3
  END INTERFACE contract

CONTAINS

  FUNCTION t2_tensor ( r, r3 ) RESULT ( t2 )
    IMPLICIT NONE
    REAL, DIMENSION(3,3)              :: t2 ! Returns second-rank 3x3 interaction tensor
    REAL, DIMENSION(3),   INTENT(in)  :: r  ! Unit vector from 2 to 1
    REAL,                 INTENT(in)  :: r3 ! Third power of the modulus of r12

    INTEGER :: i

    t2 = 3.0 * outer_product ( r, r ) ! Starting point

    FORALL (i=1:3) t2(i,i) = t2(i,i) - 1.0 ! Make traceless

    t2 = t2 / r3 ! Scale by third power of distance

  END FUNCTION t2_tensor

  FUNCTION t3_tensor ( r, r4 ) RESULT ( t3 )
    IMPLICIT NONE
    REAL, DIMENSION(3,3,3)              :: t3 ! Returns third-rank 3x3x3 interaction tensor (note positive sign)
    REAL, DIMENSION(3),     INTENT(in)  :: r  ! Unit vector from 2 to 1
    REAL,                   INTENT(in)  :: r4 ! Fourth power of the modulus of r12

    INTEGER :: i, j

    t3 = 15.0 * outer_product ( r, r, r ) ! Starting point

    DO i = 1, 3

       t3(i,i,i) = t3(i,i,i) - 9.0 * r(i) ! Correction for all indices the same

       DO j = 1, 3
          IF ( j == i ) CYCLE
          t3(i,i,j) = t3(i,i,j) - 3.0 * r(j) ! Correction for two indices the same
          t3(i,j,i) = t3(i,j,i) - 3.0 * r(j) ! Correction for two indices the same
          t3(j,i,i) = t3(j,i,i) - 3.0 * r(j) ! Correction for two indices the same
       END DO

    END DO

    t3 = t3 / r4 ! Scale by fourth power of distance

  END FUNCTION t3_tensor

  FUNCTION t4_tensor ( r, r5 ) RESULT ( t4 )
    IMPLICIT NONE
    REAL, DIMENSION(3,3,3,3)            :: t4 ! Returns fourth-rank 3x3x3x3 interaction tensor
    REAL, DIMENSION(3),     INTENT(in)  :: r  ! Unit vector from 2 to 1
    REAL,                   INTENT(in)  :: r5 ! Fifth power of the modulus of r12

    INTEGER :: i, j, k, l

    ! Define 3x3 unit matrix or Kronecker delta
    REAL, DIMENSION(3,3), PARAMETER :: u = RESHAPE([ 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0 ],[3,3])

    t4 = 105.0 * outer_product ( r, r, r, r ) ! Starting point

    DO i = 1, 3
       DO j = 1, 3
          DO k = 1, 3
             DO l = 1, 3
                t4(i,j,k,l) = t4(i,j,k,l) - 15.0 * ( &
                     &   r(i) * r(j) * u(k,l) + r(i) * r(k) * u(j,l) &
                     & + r(i) * r(l) * u(j,k) + r(j) * r(k) * u(i,l) &
                     & + r(j) * r(l) * u(i,k) + r(k) * r(l) * u(i,j) ) &
                     & + 3.0 * ( u(i,j) * u(k,l) + u(i,k) * u(j,l) + u(i,l) * u(j,k) )
             END DO
          END DO
       END DO
    END DO

    t4 = t4 / r5 ! Scale by fifth power of distance

  END FUNCTION t4_tensor

  FUNCTION t5_tensor ( r, r6 ) RESULT ( t5 )
    IMPLICIT NONE
    REAL, DIMENSION(3,3,3,3,3)          :: t5 ! Returns fifth-rank 3x3x3x3X3 interaction tensor
    REAL, DIMENSION(3),     INTENT(in)  :: r  ! Unit vector from 2 to 1
    REAL,                   INTENT(in)  :: r6 ! Sixth power of the modulus of r12

    INTEGER :: i, j, k, l, m

    ! Define 3x3 unit matrix or Kronecker delta
    REAL, DIMENSION(3,3), PARAMETER :: u = RESHAPE([ 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0 ],[3,3])

    t5 = 945.0 * outer_product ( r, r, r, r, r ) ! Starting point

    DO i = 1, 3
       DO j = 1, 3
          DO k = 1, 3
             DO l = 1, 3
                DO m = 1, 3
                   t5(i,j,k,l,m) = t5(i,j,k,l,m) - 105.0 * ( &
                        &   r(i) * r(j) * r(k) * u(l,m) + r(i) * r(j) * r(l) * u(k,m)     &
                        & + r(i) * r(j) * r(m) * u(k,l) + r(i) * r(k) * r(l) * u(j,m)     &
                        & + r(i) * r(k) * r(m) * u(j,l) + r(i) * r(l) * r(m) * u(j,k)     &
                        & + r(j) * r(k) * r(l) * u(i,m) + r(j) * r(k) * r(m) * u(i,l)     &
                        & + r(j) * r(l) * r(m) * u(i,k) + r(k) * r(l) * r(m) * u(i,j) )   &
                        & + 15.0 * ( &
                        &   r(i) * ( u(j,k) * u(l,m) + u(j,l) * u(k,m) + u(j,m) * u(k,l) ) &
                        & + r(j) * ( u(i,k) * u(l,m) + u(i,l) * u(k,m) + u(i,m) * u(k,l) ) &
                        & + r(k) * ( u(i,j) * u(l,m) + u(i,l) * u(j,m) + u(i,m) * u(j,l) ) &
                        & + r(l) * ( u(i,j) * u(k,m) + u(i,k) * u(j,m) + u(i,m) * u(j,k) ) &
                        & + r(m) * ( u(i,j) * u(k,l) + u(i,k) * u(j,l) + u(i,l) * u(j,k) ) )
                END DO
             END DO
          END DO
       END DO
    END DO

    t5 = t5 / r6 ! Scale by sixth power of distance

  END FUNCTION t5_tensor

  FUNCTION contract_11_0 ( a, b ) RESULT (c)
    IMPLICIT NONE
    REAL                           :: c ! Returns a zero-rank contraction
    REAL, DIMENSION(3), INTENT(in) :: a ! of a first-rank tensor
    REAL, DIMENSION(3), INTENT(in) :: b ! with a first-rank tensor

    c = dot_PRODUCT ( a, b )

  END FUNCTION contract_11_0

  FUNCTION contract_21_1 ( a, b ) RESULT (c)
    IMPLICIT NONE
    REAL, DIMENSION(3)                :: c ! Returns a first-rank contraction
    REAL, DIMENSION(3,3), INTENT(in)  :: a ! of a second-rank tensor
    REAL, DIMENSION(3),   INTENT(in)  :: b ! with a first-rank tensor

    c = MATMUL ( a, b )

  END FUNCTION contract_21_1

  FUNCTION contract_22_0 ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL                              :: c ! Returns a zero-rank contraction
    REAL, DIMENSION(3,3), INTENT(in)  :: a ! of a second-rank tensor
    REAL, DIMENSION(3,3), INTENT(in)  :: b ! with another second-rank tensor

    c = SUM ( a * b )

  END FUNCTION contract_22_0

  FUNCTION contract_31_2 ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3,3)                :: c ! Returns a second-rank contraction of 
    REAL, DIMENSION(3,3,3), INTENT(in)  :: a ! a third-rank tensor
    REAL, DIMENSION(3),     INTENT(in)  :: b ! and a first-rank tensor

    INTEGER :: i, j

    DO i = 1, 3
       DO j = 1, 3
          c(i,j) = SUM ( a(i,j,:) * b(:) )
       END DO
    END DO

  END FUNCTION contract_31_2

  FUNCTION contract_32_1 ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3)                  :: c ! Returns a first-rank contraction of 
    REAL, DIMENSION(3,3,3), INTENT(in)  :: a ! a third-rank tensor
    REAL, DIMENSION(3,3),   INTENT(in)  :: b ! and a second-rank tensor

    INTEGER :: i

    DO i = 1, 3
       c(i) = SUM ( a(i,:,:) * b(:,:) ) 
    END DO

  END FUNCTION contract_32_1

  FUNCTION contract_42_2 ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3,3)                 :: c ! Returns a second-rank contraction of 
    REAL, DIMENSION(3,3,3,3), INTENT(in) :: a ! a fourth-rank tensor
    REAL, DIMENSION(3,3),     INTENT(in) :: b ! and a second-rank tensor

    INTEGER :: i, j 

    DO i = 1, 3
       DO j = 1, 3
          c(i,j) = SUM ( a(i,j,:,:) * b(:,:) ) 
       END DO
    END DO

  END FUNCTION contract_42_2

  FUNCTION contract_52_3 ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3,3,3)                 :: c ! Returns a third-rank contraction of 
    REAL, DIMENSION(3,3,3,3,3), INTENT(in) :: a ! a fifth-rank tensor
    REAL, DIMENSION(3,3),   INTENT(in)     :: b ! and a second-rank tensor

    INTEGER :: i, j, k

    DO i = 1, 3
       DO j = 1, 3
          DO k = 1, 3
             c(i,j,k) = SUM ( a(i,j,k,:,:) * b(:,:) )
          END DO
       END DO
    END DO

  END FUNCTION contract_52_3

  FUNCTION skew_contract ( a, b ) RESULT ( c )

    IMPLICIT NONE
    REAL, DIMENSION(3)   :: c ! Returns a first-rank tensor by doubly contracting the Levi-Civita symbol with the matrix product of
    REAL, DIMENSION(3,3) :: a ! a second-rank tensor
    REAL, DIMENSION(3,3) :: b ! and a second-rank tensor

    REAL, DIMENSION(3,3) :: x

    ! This function evaluates epsilon_abc A_bd B_cd (both A and B are symmetric)
    ! as needed for the torque calculation

    x    = MATMUL ( a, b )
    c(1) = x(2,3) - x(3,2)
    c(2) = x(3,1) - x(1,3)
    c(3) = x(1,2) - x(2,1)

  END FUNCTION skew_contract

END MODULE t_tensor_module
