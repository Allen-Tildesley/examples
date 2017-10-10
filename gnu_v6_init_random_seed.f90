! gnu_v6_init_random_seed.f90
! Specialized random number initialization for old version of GNU compiler

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

! This routine, and the next one, are taken from the online GNU documentation
! https://gcc.gnu.org/onlinedocs/gcc-6.2.0/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
! They are specific to the gfortran compiler v6 (i.e. the version that comes with gcc6)
! For this version, calling RANDOM_SEED() initializes the random number generator
! with the same random seed to a default state, which may result in the same sequence
! being generated every time. The routines below are intended to generate different
! sequences on different calls. If you have this compiler version, you may like to replace the
! routine init_random_seed in maths_module.f90 with the following two routines.

! In gfortran v7 the random number generator was changed.
! Now, calling RANDOM_SEED() initializes the random number generator with random data
! retrieved from the operating system. Various other compilers behave the same way.
! In this case, the following two routines are redundant.

! YOU SHOULD INVESTIGATE THE BEHAVIOUR FOR YOUR OWN COMPILER AND MACHINE IMPLEMENTATION 

SUBROUTINE init_random_seed
  USE iso_fortran_env, ONLY: int64
  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: seed(:)
  INTEGER :: i, n, un, istat, dt(8), pid
  INTEGER(int64) :: t

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  ! First try if the OS provides a random number generator
  OPEN(newunit=un, file='/dev/urandom', access='stream', &
       form='unformatted', action='read', status='old', iostat=istat)
  IF (istat == 0) THEN
     READ(un) seed
     CLOSE(un)
  ELSE
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     CALL SYSTEM_CLOCK(t)
     IF (t == 0) THEN
        CALL DATE_AND_TIME(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     END IF
     pid = getpid() ! apparently getpid is declared somewhere else in GCC
     t = IEOR(t, INT(pid, KIND(t)))
     DO i = 1, n
        seed(i) = lcg(t)
     END DO
  END IF
  CALL RANDOM_SEED(put=seed)
END SUBROUTINE init_random_seed

! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
FUNCTION lcg(s)
  USE iso_fortran_env, ONLY: int64
  IMPLICIT NONE
  INTEGER :: lcg
  INTEGER(int64) :: s
  IF (s == 0) THEN
     s = 104729
  ELSE
     s = MOD(s, 4294967296_int64)
  END IF
  s = MOD(s * 279470273_int64, 4294967291_int64)
  lcg = INT(MOD(s, INT(HUGE(0), int64)), KIND(0))
END FUNCTION lcg
