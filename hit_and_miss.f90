! hit_and_miss.f90
! Estimates volume of polyhedron by simple MC
PROGRAM hit_and_miss

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

  REAL                          :: v
  REAL, DIMENSION(3)            :: r, zeta
  REAL, DIMENSION(3), PARAMETER :: r_0 = [1.0, 2.0, 3.0]
  REAL,               PARAMETER :: v_0 = PRODUCT(r_0)
  INTEGER                       :: tau, tau_shot, tau_hit

  CALL RANDOM_SEED()
  tau_hit  = 0
  tau_shot = 1000000

  DO tau = 1, tau_shot
     CALL RANDOM_NUMBER ( zeta(:) ) ! uniform in range (0,1)
     r = zeta * r_0                 ! uniform in v_0
     IF (   r(2) < ( 2.0 - 2.0*r(1) ) .AND. &
          & r(3) < ( 1.0 + r(2) )   ) THEN ! in polyhedron
        tau_hit = tau_hit + 1
     END IF
  END DO
  v = v_0 * REAL ( tau_hit ) / REAL ( tau_shot ) 
  WRITE ( unit=output_unit, fmt='(a,f10.5)' ) 'Estimate = ', v

END PROGRAM hit_and_miss
