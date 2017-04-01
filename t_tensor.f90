! t_tensor.f90
! Electrostatic interactions: T-tensors compared with angles
PROGRAM t_tensor

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

  ! The dipole moment of molecule 1 is aligned along the axial vector e1
  ! The quadrupole tensor, quad1, is diagonal and traceless with
  ! quad1_xx = -0.5*quad1_mag, quad1_yy = -0.5*quad1_mag, and quad1_zz = quad1_mag
  ! in the molecule-fixed system. Similarly for molecule 2.
  ! The vector r12 = r1-r2 points from 2 to 1.

  ! Forces are calculated by differentiating the T-tensor, giving the next higher rank T-tensor
  ! Torques are calculated from the angular dependence of dipole, quadrupole etc.
  ! potential V = mu_a g_a => torque tau_a = -epsilon_abc mu_b g_c
  ! potential V = Q_ab G_ab => torque tau_a = -2 epsilon_abc Q_bd G_cd
  ! where abcd are Cartesian indices and epsilon is the Levi-Civita symbol
  ! It is just necessary to identify the constants g_a, G_ab, in terms of the T tensor and the
  ! multipole on the other molecule.

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE maths_module,    ONLY : init_random_seed, random_vector, outer_product, cross_product
  use t_tensor_module, only : t2_tensor, t3_tensor, t4_tensor, t5_tensor, contract, skew_contract

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
  g    =  contract ( tt2, mu2 ) ! Contract T2 with dipole 2
  v12t = -contract ( mu1, g   ) ! Contract result with dipole 1
  v12e = (mu1_mag*mu2_mag/r12_mag**3) * ( c12 - 3.0 * c1 * c2 ) ! Compare result from angles
  WRITE ( unit=output_unit, fmt='(a,t30,f12.6,t70,f12.6,t110,es10.2)' ) 'Energy =', v12t, v12e, v12t-v12e

  ! Calculate the dipole-dipole force
  gg   =  contract ( tt3, mu2 ) ! Contract T3 with dipole 2
  f12t = -contract ( gg,  mu1 ) ! Contract result with dipole 1
  f12e = (3.0*mu1_mag*mu2_mag/r12_mag**4) * ( (c12-5.0*c1*c2)*r12_hat + c2*e1 + c1*e2 ) ! Compare result from angles
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Force  =', f12t, f12e, f12t-f12e

  ! Calculate the dipole-dipole torques
  g   = -contract ( tt2, mu2 )    ! Contract T2 with dipole 2
  t1t = -cross_product ( mu1, g ) ! Cross-product result with dipole 1
  g   = e2 - 3.0*c2*r12_hat ! Compare result from angles
  t1e = -(mu1_mag*mu2_mag/r12_mag**3) * cross_product ( e1, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 1  =', t1t, t1e, t1t-t1e
  g   = -contract ( tt2, mu1 )    ! Contract T2 with dipole 1
  t2t = -cross_product ( mu2, g ) ! Cross-product result with dipole 2
  g   = e1 - 3.0*c1 * r12_hat ! Compare result from angles
  t2e = -(mu1_mag*mu2_mag/r12_mag**3) * cross_product ( e2, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 2  =', t2t, t2e, t2t-t2e

  WRITE ( unit=output_unit, fmt='(/,a)') 'Dipole-quadrupole'

  ! Calculate the dipole-quadrupole energy
  g    = contract ( tt3, quad2 )           ! Contract T3 with quadrupole 2
  v12t = -(1.0/3.0) * contract ( mu1, g ) ! Contract result with dipole 1
  v12e = (1.5*mu1_mag*quad2_mag/r12_mag**4) * ( c1*(1.0-5.0*c2*c2) + 2.0*c2*c12 )
  WRITE ( unit=output_unit, fmt='(a,t30,f12.6,t70,f12.6,t110,es10.2)' ) 'Energy =', v12t, v12e, v12t-v12e

  ! Calculate the dipole-quadrupole force
  gg   = contract ( tt4, quad2 )           ! Contract T4 with quadrupole 2
  f12t = -(1.0/3.0) * contract ( gg, mu1 ) ! Contract result with dipole 1
  f12e = -(1.5*mu1_mag*quad2_mag/r12_mag**5) * ( & ! Compare result from angles
       &     ( 35.0*c1*c2**2 - 10.0*c2*c12 - 5.0*c1 ) * r12_hat  &
       &   + ( 1.0 - 5.0*c2**2 ) * e1 &
       &   + ( 2.0*c12 - 10.0*c1*c2 ) * e2 )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Force  =', f12t, f12e, f12t-f12e

  ! Calculate the dipole-quadrupole torques
  g   = -(1.0/3.0)*contract ( tt3, quad2 ) ! Contract T3 with quadrupole 2
  t1t = -cross_product ( mu1, g )          ! Cross-product result with dipole 1
  g   =  (1.0-5.0*c2**2) * r12_hat + 2.0*c2 * e2 ! Compare result from angles
  t1e = -(1.5*mu1_mag*quad2_mag/r12_mag**4) * cross_product ( e1, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 1  =', t1t, t1e, t1t-t1e
  gg  = -(1.0/3.0)*contract ( tt3, mu1 ) ! Contract T3 with dipole 1
  t2t = -2.0*skew_contract ( quad2, gg ) ! Skew-contract result with quadrupole 2
  g   =  (c12-5.0*c1*c2) * r12_hat + c2 * e1 ! Compare result from angles
  t2e = -(3.0*mu1_mag*quad2_mag/r12_mag**4) * cross_product ( e2, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 2  =', t2t, t2e, t2t-t2e

  WRITE ( unit=output_unit, fmt='(/,a)') 'Quadrupole-dipole'

  ! Calculate the quadrupole-dipole energy
  g    = contract ( tt3, quad1 )           ! Contract T3 with quadrupole 1
  v12t =  (1.0/3.0) * contract ( g, mu2 ) ! Contract result with dipole 2
  v12e = -(1.5*quad1_mag*mu2_mag/r12_mag**4) * ( c2*(1.0-5.0*c1**2) + 2.0*c1*c12 ) ! Compare result from angles
  WRITE ( unit=output_unit, fmt='(a,t30,f12.6,t70,f12.6,t110,es10.2)' ) 'Energy =', v12t, v12e, v12t-v12e

  ! Calculate the quadrupole-dipole force
  gg   = contract ( tt4, quad1 )          ! Contract T4 with quadrupole 1 
  f12t = (1.0/3.0) * contract ( gg, mu2 ) ! Contract result with dipole 2
  f12e = (1.5*quad1_mag*mu2_mag/r12_mag**5) * ( & ! Compare result from angles
       &     ( 35.0*c2*c1**2 - 10.0*c1*c12 - 5.0*c2 ) * r12_hat  &
       &   + ( 1.0-5.0*c1**2 ) * e2 &
       &   + ( 2.0*c12 - 10.0*c1*c2 ) * e1 )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Force  =', f12t, f12e, f12t-f12e

  ! Calculate the quadrupole-dipole torques
  gg  = (1.0/3.0)*contract ( tt3, mu2 )  ! Contract T3 with dipole 2
  t1t = -2.0*skew_contract ( quad1, gg ) ! Skew-contract result with quadrupole 1
  g   = (c12-5.0*c1*c2) * r12_hat + c1 * e2 ! Compare result from angles
  t1e = (3.0*quad1_mag*mu2_mag/r12_mag**4) * cross_product ( e1, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 1  =', t1t, t1e, t1t-t1e
  g   = (1.0/3.0)*contract ( tt3, quad1 ) ! Contract T3 with quadrupole 1
  t2t = -cross_product ( mu2, g )         ! Cross product result with dipole 2
  g   = (1.0-5.0*c1**2) * r12_hat + 2.0*c1 * e1 ! Compare result from angles
  t2e = (1.5*quad1_mag*mu2_mag/r12_mag**4) * cross_product ( e2, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 2  =', t2t, t2e, t2t-t2e

  WRITE ( unit=output_unit, fmt='(/,a)') 'Quadrupole-quadrupole'

  ! Calculate the quadrupole-quadrupole energy
  gg   = contract ( tt4, quad2 )            ! Contract T4 with quadrupole 2
  v12t = (1.0/9.0) * contract ( quad1, gg ) ! Contract result with quadrupole 1
  v12e = (0.75*quad1_mag*quad2_mag/r12_mag**5) * ( & ! Compare result from angles
       &    1.0 - 5.0*c1**2 - 5.0*c2**2 + 2.0*c12**2 + 35.0*(c1*c2)**2 - 20.0*c1*c2*c12 )
  WRITE ( unit=output_unit, fmt='(a,t30,f12.6,t70,f12.6,t110,es10.2)' ) 'Energy =', v12t, v12e, v12t-v12e

  ! Calculate the quadrupole-quadrupole force
  ggg  = contract ( tt5, quad2 )             ! Contract T5 with quadrupole 2
  f12t = (1.0/9.0) * contract ( ggg, quad1 ) ! Contract result with quadrupole 1
  f12e = (0.75*quad1_mag*quad2_mag/r12_mag**6) * ( & ! Compare result from angles
       &    ( 5.0 - 35.0*c1**2 - 35.0*c2**2 + 10.0*c12**2 + 315.0*(c1*c2)**2 - 140.0*c1*c2*c12 ) * r12_hat  &
       &  + ( 10.0*c1 - 70.0*c1*c2**2 + 20.0*c2*c12 ) * e1 &
       &  + ( 10.0*c2 - 70.0*c2*c1**2 + 20.0*c1*c12 ) * e2   )          
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Force  =', f12t, f12e, f12t-f12e

  ! Calculate the quadrupole-quadrupole torques
  gg  = (1.0/9.0)*contract ( tt4, quad2 ) ! Contract T4 with quadrupole 2
  t1t = -2.0*skew_contract(quad1, gg)     ! Skew-contract result with quadrupole 1
  g   = 2.5*(c1*(7.0*c2**2-1.0)-2.0*c2*c12) * r12_hat - (5.0*c1*c2-c12) * e2 ! Compare result from angles
  t1e = -(3.0*quad1_mag*quad2_mag/r12_mag**5) * cross_product ( e1, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 1  =', t1t, t1e, t1t-t1e
  gg  = (1.0/9.0)*contract ( tt4, quad1 ) ! Contract T4 with quadrupole 1
  t2t = -2.0*skew_contract(quad2, gg)     ! Skew-contract result with quadrupole 2
  g   = 2.5*(c2*(7.0*c1**2-1.0)-2.0*c1*c12) * r12_hat -(5.0*c1*c2-c12) * e1 ! Compare result from angles
  t2e = -(3.0*quad1_mag*quad2_mag/r12_mag**5) * cross_product ( e2, g )
  WRITE ( unit=output_unit, fmt='(a,t30,3f12.6,t70,3f12.6,t110,3es10.2)' ) 'Torque on 2  =', t2t, t2e, t2t-t2e

END PROGRAM t_tensor
