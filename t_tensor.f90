! t_tensor.f90
! Electrostatic interactions: T-tensors compared with angles
PROGRAM t_tensor

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

  ! The dipole moment of molecule 1 is aligned along the axial vector e1
  ! The quadrupole tensor, quad1, is diagonal and traceless with
  ! quad1_xx = -0.5*quad1_mag, quad1_yy = -0.5*quad1_mag, and quad1_zz = quad1_mag
  ! in the molecule-fixed system. Similarly for molecule 2.
  ! The vector r12 = r1-r2 points from 2 to 1.

  ! Forces are calculated by differentiating the T-tensor, giving the next higher rank T-tensor
  ! Torques are calculated from the angular dependence of dipole, quadrupole etc.
  ! potential V = mu_i g_i => torque tau_i = -epsilon_ijk mu_j g_k (=-cross_product)
  ! potential V = Q_ij G_ij => torque tau_l = -2 epsilon_lij Q_ik G_jk
  ! where ijkl are Cartesian indices and epsilon is the Levi-Civita symbol
  ! It is just necessary to identify the constants g_i, G_ij, in terms of the T tensor and the
  ! multipole on the other molecule.

  ! NB in the text, eqn (1.15), the signs of the odd-rank terms in the energy are wrong.
  ! See https://github.com/Allen-Tildesley/corrections. The correct formulae are used here.

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE maths_module,    ONLY : init_random_seed, random_vector, outer_product, cross_product

  IMPLICIT NONE

  REAL, DIMENSION(3)         :: r12, r12_hat, e1, e2, mu1, mu2, f12t, f12e, t1t, t2t, t1e, t2e, g
  REAL, DIMENSION(3,3)       :: tt2, quad1, quad2, gg
  REAL, DIMENSION(3,3,3)     :: tt3, ggg
  REAL, DIMENSION(3,3,3,3)   :: tt4
  REAL, DIMENSION(3,3,3,3,3) :: tt5

  REAL    :: r12_mag, c1, c2, c12, v12t, v12e
  INTEGER :: i, ioerr
  REAL    :: d_min, d_max, mu1_mag, mu2_mag, quad1_mag, quad2_mag

  NAMELIST /nml/ d_min, d_max, mu1_mag, mu2_mag, quad1_mag, quad2_mag

  WRITE ( unit=output_unit, fmt='(a)' ) 'T-tensor'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Calculation of electrostatic interactions between linear molecules' 
  WRITE ( unit=output_unit, fmt='(a)' ) 'using T-tensors and Euler angles'

  ! Initialize random number generator                          
  CALL init_random_seed

  ! Default parameters
  d_min     = 0.5 ! Minimum separation
  d_max     = 1.5 ! Maximum separation
  mu1_mag   = 1.0 ! Dipole moment of molecule 1    
  mu2_mag   = 1.0 ! Dipole moment of molecule 2    
  quad1_mag = 1.0 ! Quadrupole moment of molecule 1
  quad2_mag = 1.0 ! Quadrupole moment of molecule 2

  !Read parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in t_tensor'
  END IF

  ! Write out parameters
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'  ) 'Min separation d_min',            d_min
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'  ) 'Max separation d_max',            d_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'  ) 'Dipole moment of molecule 1',     mu1_mag   
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'  ) 'Dipole moment of molecule 2',     mu2_mag   
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'  ) 'Quadrupole moment of molecule 1', quad1_mag 
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'  ) 'Quadrupole moment of molecule 2', quad2_mag 

  ! Choose orientations at random
  e1 = random_vector ( )
  e2 = random_vector ( )

  ! Place atom 2 at origin and atom 1 in a random direction within desired distance range
  r12_hat = random_vector ( ) ! unit vector
  CALL RANDOM_NUMBER ( r12_mag )
  r12_mag = d_min + (d_max-d_min)*r12_mag ! Magnitude of r12
  r12     = r12_hat * r12_mag             ! Within desired range of origin

  c1  = DOT_PRODUCT ( e1, r12_hat ) ! Cosine of angle between e1 and r12
  c2  = DOT_PRODUCT ( e2, r12_hat ) ! Cosine of angle between e2 and r12
  c12 = DOT_PRODUCT ( e1, e2   )    ! Cosine of angle between e1 and e2

  WRITE ( unit=output_unit, fmt='(a,t40,3f12.6)' ) 'Displacement r12 = ', r12
  WRITE ( unit=output_unit, fmt='(a,t40,3f12.6)' ) 'Orientation  e1  = ', e1
  WRITE ( unit=output_unit, fmt='(a,t40,3f12.6)' ) 'Orientation  e2  = ', e2

  ! Dipole vectors in space-fixed frame
  mu1 = mu1_mag * e1                             
  mu2 = mu2_mag * e2

  ! Quadrupole tensors in space-fixed frame (traceless)
  quad1 = 1.5 * outer_product ( e1, e1 )        
  FORALL (i=1:3) quad1(i,i) = quad1(i,i) - 0.5
  quad1 = quad1_mag * quad1
  quad2 = 1.5 * outer_product ( e2, e2 )
  FORALL (i=1:3) quad2(i,i) = quad2(i,i) - 0.5
  quad2 = quad2_mag * quad2

  ! The T tensors of each rank: T2, T3, T4, T5
  tt2 = t2_tensor ( r12_hat, r12_mag**3 )
  tt3 = t3_tensor ( r12_hat, r12_mag**4 )
  tt4 = t4_tensor ( r12_hat, r12_mag**5 )
  tt5 = t5_tensor ( r12_hat, r12_mag**6 )

  ! Headings
  WRITE ( unit=output_unit, fmt='(/,t30,a36,t70,a36,t110,a30)' )  &
       & '.....Result from T tensor', '.....Result from Euler angles', '.........Difference'

  WRITE ( unit=output_unit, fmt='(/,a)') 'Dipole-dipole'

  ! Calculate the dipole-dipole energy
  g    =  contract_ij_j ( tt2, mu2 ) ! Contract T2 with dipole 2
  v12t = -contract_i_i  ( mu1, g   ) ! Contract result with dipole 1
  v12e = (mu1_mag*mu2_mag/r12_mag**3) * ( c12 - 3.0 * c1 * c2 ) ! Compare result from angles
  WRITE ( unit=output_unit, fmt='(a,t30,f12.6,t70,f12.6,t110,es10.2)' ) 'Energy =', v12t, v12e, v12t-v12e

  ! Calculate the dipole-dipole force
  gg   =  contract_ijk_k ( tt3, mu2 ) ! Contract T3 with dipole 2
  f12t = -contract_ij_j  ( gg,  mu1 ) ! Contract result with dipole 1
  f12e = (3.0*mu1_mag*mu2_mag/r12_mag**4) * ( (c12-5.0*c1*c2)*r12_hat + c2*e1 + c1*e2 ) ! Compare result from angles
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Force  =', f12t, f12e, f12t-f12e

  ! Calculate the dipole-dipole torques
  g   = -contract_ij_j ( tt2, mu2 ) ! Contract T2 with dipole 2
  t1t = -cross_product ( mu1, g   ) ! Cross-product result with dipole 1
  g   = e2 - 3.0*c2*r12_hat ! Compare result from angles
  t1e = -(mu1_mag*mu2_mag/r12_mag**3) * cross_product ( e1, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 1  =', t1t, t1e, t1t-t1e
  g   = -contract_ij_j ( tt2, mu1 ) ! Contract T2 with dipole 1
  t2t = -cross_product ( mu2, g )   ! Cross-product result with dipole 2
  g   = e1 - 3.0*c1 * r12_hat ! Compare result from angles
  t2e = -(mu1_mag*mu2_mag/r12_mag**3) * cross_product ( e2, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 2  =', t2t, t2e, t2t-t2e

  WRITE ( unit=output_unit, fmt='(/,a)') 'Dipole-quadrupole'

  ! Calculate the dipole-quadrupole energy
  g    = contract_ijk_jk ( tt3, quad2 )           ! Contract T3 with quadrupole 2
  v12t = -(1.0/3.0) * contract_i_i ( mu1, g     ) ! Contract result with dipole 1
  v12e = (1.5*mu1_mag*quad2_mag/r12_mag**4) * ( c1*(1.0-5.0*c2*c2) + 2.0*c2*c12 )
  WRITE ( unit=output_unit, fmt='(a,t30,f12.6,t70,f12.6,t110,es10.2)' ) 'Energy =', v12t, v12e, v12t-v12e

  ! Calculate the dipole-quadrupole force
  gg   = contract_ijkl_kl ( tt4, quad2 )        ! Contract T4 with quadrupole 2
  f12t = -(1.0/3.0) * contract_ij_j ( gg, mu1 ) ! Contract result with dipole 1
  f12e = -(1.5*mu1_mag*quad2_mag/r12_mag**5) * ( & ! Compare result from angles
       &     ( 35.0*c1*c2**2 - 10.0*c2*c12 - 5.0*c1 ) * r12_hat  &
       &   + ( 1.0 - 5.0*c2**2 ) * e1 &
       &   + ( 2.0*c12 - 10.0*c1*c2 ) * e2 )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Force  =', f12t, f12e, f12t-f12e

  ! Calculate the dipole-quadrupole torques
  g   = -(1.0/3.0)*contract_ijk_jk ( tt3, quad2 ) ! Contract T3 with quadrupole 2
  t1t = -cross_product ( mu1, g )                 ! Cross-product result with dipole 1
  g   =  (1.0-5.0*c2**2) * r12_hat + 2.0*c2 * e2 ! Compare result from angles
  t1e = -(1.5*mu1_mag*quad2_mag/r12_mag**4) * cross_product ( e1, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 1  =', t1t, t1e, t1t-t1e
  gg  = -(1.0/3.0)*contract_ijk_k ( tt3, mu1 ) ! Contract T3 with dipole 1
  gg = contract_ik_jk ( quad2, gg )            ! Contract result with quadrupole 2
  t2t = -2.0*skew ( gg )                       ! Contract with Levi-Civita symbol
  g   =  (c12-5.0*c1*c2) * r12_hat + c2 * e1 ! Compare result from angles
  t2e = -(3.0*mu1_mag*quad2_mag/r12_mag**4) * cross_product ( e2, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 2  =', t2t, t2e, t2t-t2e

  WRITE ( unit=output_unit, fmt='(/,a)') 'Quadrupole-dipole'

  ! Calculate the quadrupole-dipole energy
  g    = contract_ijk_jk ( tt3, quad1 )      ! Contract T3 with quadrupole 1
  v12t = (1.0/3.0) * contract_i_i ( g, mu2 ) ! Contract result with dipole 2
  v12e = -(1.5*quad1_mag*mu2_mag/r12_mag**4) * ( c2*(1.0-5.0*c1**2) + 2.0*c1*c12 ) ! Compare result from angles
  WRITE ( unit=output_unit, fmt='(a,t30,f12.6,t70,f12.6,t110,es10.2)' ) 'Energy =', v12t, v12e, v12t-v12e

  ! Calculate the quadrupole-dipole force
  gg   = contract_ijkl_kl ( tt4, quad1 )       ! Contract T4 with quadrupole 1 
  f12t = (1.0/3.0) * contract_ij_j ( gg, mu2 ) ! Contract result with dipole 2
  f12e = (1.5*quad1_mag*mu2_mag/r12_mag**5) * ( & ! Compare result from angles
       &     ( 35.0*c2*c1**2 - 10.0*c1*c12 - 5.0*c2 ) * r12_hat  &
       &   + ( 1.0-5.0*c1**2 ) * e2 &
       &   + ( 2.0*c12 - 10.0*c1*c2 ) * e1 )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Force  =', f12t, f12e, f12t-f12e

  ! Calculate the quadrupole-dipole torques
  gg  = (1.0/3.0)*contract_ijk_k ( tt3, mu2 ) ! Contract T3 with dipole 2
  gg  = contract_ik_jk ( quad1, gg )          ! Contract result with quadrupole 1
  t1t = -2.0*skew ( gg )                      ! Contract with Levi-Civita symbol
  g   = (c12-5.0*c1*c2) * r12_hat + c1 * e2 ! Compare result from angles
  t1e = (3.0*quad1_mag*mu2_mag/r12_mag**4) * cross_product ( e1, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 1  =', t1t, t1e, t1t-t1e
  g   = (1.0/3.0)*contract_ijk_jk ( tt3, quad1 ) ! Contract T3 with quadrupole 1
  t2t = -cross_product ( mu2, g )         ! Cross product result with dipole 2
  g   = (1.0-5.0*c1**2) * r12_hat + 2.0*c1 * e1 ! Compare result from angles
  t2e = (1.5*quad1_mag*mu2_mag/r12_mag**4) * cross_product ( e2, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 2  =', t2t, t2e, t2t-t2e

  WRITE ( unit=output_unit, fmt='(/,a)') 'Quadrupole-quadrupole'

  ! Calculate the quadrupole-quadrupole energy
  gg   = contract_ijkl_kl ( tt4, quad2 )          ! Contract T4 with quadrupole 2
  v12t = (1.0/9.0) * contract_ij_ij ( quad1, gg ) ! Contract result with quadrupole 1
  v12e = (0.75*quad1_mag*quad2_mag/r12_mag**5) * ( & ! Compare result from angles
       &    1.0 - 5.0*c1**2 - 5.0*c2**2 + 2.0*c12**2 + 35.0*(c1*c2)**2 - 20.0*c1*c2*c12 )
  WRITE ( unit=output_unit, fmt='(a,t30,f12.6,t70,f12.6,t110,es10.2)' ) 'Energy =', v12t, v12e, v12t-v12e

  ! Calculate the quadrupole-quadrupole force
  ggg  = contract_ijklm_lm ( tt5, quad2 )           ! Contract T5 with quadrupole 2
  f12t = (1.0/9.0) * contract_ijk_jk ( ggg, quad1 ) ! Contract result with quadrupole 1
  f12e = (0.75*quad1_mag*quad2_mag/r12_mag**6) * ( & ! Compare result from angles
       &    ( 5.0 - 35.0*c1**2 - 35.0*c2**2 + 10.0*c12**2 + 315.0*(c1*c2)**2 - 140.0*c1*c2*c12 ) * r12_hat  &
       &  + ( 10.0*c1 - 70.0*c1*c2**2 + 20.0*c2*c12 ) * e1 &
       &  + ( 10.0*c2 - 70.0*c2*c1**2 + 20.0*c1*c12 ) * e2   )          
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Force  =', f12t, f12e, f12t-f12e

  ! Calculate the quadrupole-quadrupole torques
  gg  = (1.0/9.0)*contract_ijkl_kl ( tt4, quad2 ) ! Contract T4 with quadrupole 2
  gg  = contract_ik_jk(quad1, gg)                 ! Contract result with quadrupole 1
  t1t = -2.0*skew ( gg )                          ! Contract with Levi-Civita symbol
  g   = 2.5*(c1*(7.0*c2**2-1.0)-2.0*c2*c12) * r12_hat - (5.0*c1*c2-c12) * e2 ! Compare result from angles
  t1e = -(3.0*quad1_mag*quad2_mag/r12_mag**5) * cross_product ( e1, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 1  =', t1t, t1e, t1t-t1e
  gg  = (1.0/9.0)*contract_ijkl_kl ( tt4, quad1 ) ! Contract T4 with quadrupole 1
  gg  = contract_ik_jk(quad2, gg)                 ! Contract result with quadrupole 2
  t2t = -2.0*skew ( gg )                          ! Contract with Levi-Civita symbol
  g   = 2.5*(c2*(7.0*c1**2-1.0)-2.0*c1*c12) * r12_hat -(5.0*c1*c2-c12) * e1 ! Compare result from angles
  t2e = -(3.0*quad1_mag*quad2_mag/r12_mag**5) * cross_product ( e2, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 2  =', t2t, t2e, t2t-t2e

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

  FUNCTION contract_i_i ( a, b ) RESULT (c)
    IMPLICIT NONE
    REAL                           :: c ! Returns a zero-rank contraction
    REAL, DIMENSION(3), INTENT(in) :: a ! of a first-rank tensor
    REAL, DIMENSION(3), INTENT(in) :: b ! with a first-rank tensor

    ! Sum over index i, equivalent to built-in dot_product
    c = DOT_PRODUCT ( a, b )

  END FUNCTION contract_i_i

  FUNCTION contract_ij_j ( a, b ) RESULT (c)
    IMPLICIT NONE
    REAL, DIMENSION(3)                :: c ! Returns a first-rank contraction
    REAL, DIMENSION(3,3), INTENT(in)  :: a ! of a second-rank tensor
    REAL, DIMENSION(3),   INTENT(in)  :: b ! with a first-rank tensor

    ! Free index i, sum over index j, equivalent to built-in mat_mul
    c = MATMUL ( a, b )

  END FUNCTION contract_ij_j

  FUNCTION contract_ij_ij ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL                              :: c ! Returns a zero-rank contraction
    REAL, DIMENSION(3,3), INTENT(in)  :: a ! of a second-rank tensor
    REAL, DIMENSION(3,3), INTENT(in)  :: b ! with another second-rank tensor

    ! Sum over indices i, j
    c = SUM ( a * b )

  END FUNCTION contract_ij_ij

  FUNCTION contract_ik_jk ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3,3)             :: c ! Returns a second-rank contraction
    REAL, DIMENSION(3,3), INTENT(in) :: a ! of a second-rank tensor
    REAL, DIMENSION(3,3), INTENT(in) :: b ! with another second-rank tensor
    INTEGER :: i, j

    ! Free indices i, j, sum over index k (if b symmetric, would be equivalent to built-in matmal)
    DO i = 1, 3
       DO j = 1, 3
          c(i,j) = SUM ( a(i,:) * b(j,:) )
       END DO
    END DO

  END FUNCTION contract_ik_jk

  FUNCTION contract_ijk_k ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3,3)                :: c ! Returns a second-rank contraction of 
    REAL, DIMENSION(3,3,3), INTENT(in)  :: a ! a third-rank tensor
    REAL, DIMENSION(3),     INTENT(in)  :: b ! and a first-rank tensor

    INTEGER :: i, j

    ! Free indices i, j, sum over index k.
    DO i = 1, 3
       DO j = 1, 3
          c(i,j) = SUM ( a(i,j,:) * b(:) )
       END DO
    END DO

  END FUNCTION contract_ijk_k

  FUNCTION contract_ijk_jk ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3)                  :: c ! Returns a first-rank contraction of 
    REAL, DIMENSION(3,3,3), INTENT(in)  :: a ! a third-rank tensor
    REAL, DIMENSION(3,3),   INTENT(in)  :: b ! and a second-rank tensor

    INTEGER :: i

    ! Free index i, sum over indices j, k.
    DO i = 1, 3
       c(i) = SUM ( a(i,:,:) * b(:,:) ) 
    END DO

  END FUNCTION contract_ijk_jk

  FUNCTION contract_ijkl_kl ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3,3)                 :: c ! Returns a second-rank contraction of 
    REAL, DIMENSION(3,3,3,3), INTENT(in) :: a ! a fourth-rank tensor
    REAL, DIMENSION(3,3),     INTENT(in) :: b ! and a second-rank tensor

    INTEGER :: i, j 

    ! Free indices i,j, sum over indices k, l
    DO i = 1, 3
       DO j = 1, 3
          c(i,j) = SUM ( a(i,j,:,:) * b(:,:) ) 
       END DO
    END DO

  END FUNCTION contract_ijkl_kl

  FUNCTION contract_ijklm_lm ( a, b ) RESULT ( c )
    IMPLICIT NONE
    REAL, DIMENSION(3,3,3)                 :: c ! Returns a third-rank contraction of 
    REAL, DIMENSION(3,3,3,3,3), INTENT(in) :: a ! a fifth-rank tensor
    REAL, DIMENSION(3,3),   INTENT(in)     :: b ! and a second-rank tensor

    INTEGER :: i, j, k

    ! Free indices i,j,k, sum over indices l, m
    DO i = 1, 3
       DO j = 1, 3
          DO k = 1, 3
             c(i,j,k) = SUM ( a(i,j,k,:,:) * b(:,:) )
          END DO
       END DO
    END DO

  END FUNCTION contract_ijklm_lm

  FUNCTION skew ( a ) RESULT ( b )

    IMPLICIT NONE
    REAL, DIMENSION(3)   :: b ! Returns a first-rank tensor by contracting the Levi-Civita symbol with
    REAL, DIMENSION(3,3) :: a ! a second-rank tensor

    b(1) = a(2,3) - a(3,2)
    b(2) = a(3,1) - a(1,3)
    b(3) = a(1,2) - a(2,1)

  END FUNCTION skew

END PROGRAM t_tensor
