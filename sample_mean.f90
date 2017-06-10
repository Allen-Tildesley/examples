! sample_mean.f90
! Estimates volume of polyhedron by sample mean integration
PROGRAM sample_mean

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

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit

  IMPLICIT NONE

  REAL                          :: v, f
  REAL, DIMENSION(2)            :: r, zeta
  REAL, DIMENSION(2), PARAMETER :: r_0 = [1.0, 2.0]
  REAL,               PARAMETER :: a_0 = PRODUCT(r_0)
  INTEGER                       :: tau, tau_max

  CALL RANDOM_SEED()
  tau_max = 1000000

  f = 0.0
  DO tau = 1, tau_max
     CALL RANDOM_NUMBER ( zeta ) ! uniform in (0,1)
     r = zeta * r_0 ! uniform in xy rectangle
     IF ( r(2) < 2.0-2.0*r(1) ) THEN ! in xy triangle
        f = f + ( 1.0 + r(2) )       ! value of z
     END IF
  END DO
  v = a_0 * f / REAL ( tau_max ) 
  WRITE ( unit=output_unit, fmt='(a,f10.5)' ) 'Estimate = ', v

END PROGRAM sample_mean
