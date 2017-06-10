! wl_hist.f90
! Wang-Landau histogram post-processing
PROGRAM wl_hist

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

  ! Program to post-process histograms produced by mc_chain_wl_sw.f90
  ! Reads the following from standard input, until end-of-file:
  !    histogram filename (20 characters max)
  !    number of atoms in chain (solely to make PE and Cv intensive)
  !    one temperature per line (as many lines as you like)

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  IMPLICIT NONE

  REAL, DIMENSION(:), ALLOCATABLE :: e      ! Array of energies
  REAL, DIMENSION(:), ALLOCATABLE :: h      ! Histogram of probabilities (should be "flat")
  REAL, DIMENSION(:), ALLOCATABLE :: g      ! Histogram of average radii of gyration, resolved by energy
  REAL, DIMENSION(:), ALLOCATABLE :: s      ! Entropy values (log of DOS)
  REAL, DIMENSION(:), ALLOCATABLE :: boltz  ! Array of Boltzmann factors (including DOS)

  INTEGER           :: q, q_max, qs_max, ioerr, hist_unit, n
  REAL              :: t, s_max, norm, g_avg, e_avg, e_msd
  CHARACTER(len=20) :: hist_file

  ! Read histogram filename
  READ  ( unit=input_unit,  fmt='(a)'   ) hist_file
  WRITE ( unit=output_unit, fmt='(a,a)' ) 'Reading histograms from ', hist_file

  OPEN ( newunit=hist_unit, file=hist_file, action='read', status='old' )

  ! First pass through histogram file determining number of energy states
  ! For simplicity we number them 0:q_max, which is usually the case
  q_max = -1
  DO
     READ ( unit=hist_unit, fmt=*, iostat=ioerr )
     IF ( ioerr == iostat_end ) EXIT
     IF ( ioerr /= 0 ) STOP 'Unexpected error'
     q_max = q_max + 1
  END DO
  WRITE ( unit=output_unit, fmt='(a,t10,i15)' ) 'q_max = ', q_max

  ! Allocate arrays to exactly the right size
  ALLOCATE ( e(0:q_max), h(0:q_max), g(0:q_max), s(0:q_max), boltz(0:q_max) )

  ! Second pass through histogram file reading in data
  REWIND ( unit=hist_unit )
  DO q = 0, q_max
     READ ( unit=hist_unit, fmt=*, iostat=ioerr ) e(q), h(q), g(q), s(q)
     IF ( ioerr /= 0 ) STOP 'Error in histogram data'
  END DO

  CLOSE ( unit=hist_unit )

  READ  ( unit=input_unit, fmt=* ) n
  WRITE ( unit=output_unit, fmt='(a,t10,i15)' ) 'N = ', n

  ! Write headings
  WRITE ( unit=output_unit, fmt='(4a15)') 'T', 'Rg', 'PE/N', 'Cv(ex)/N'

  ! Read and process temperatures until end of input file
  DO
     READ  ( unit=input_unit, fmt=*, iostat=ioerr ) t
     IF ( ioerr == iostat_end ) EXIT
     IF ( ioerr /= 0 ) STOP 'Unexpected error'

     ! Locate maximum Boltzmann factor (helps avoid overflows)
     qs_max = MAXLOC ( s(0:q_max) - e(0:q_max) / t, dim=1 )
     qs_max = qs_max - 1 ! Correct for the slice numbering starting at zero
     s_max = s(qs_max)

     ! Compute Boltzmann weights including density of states
     boltz = s(0:q_max) - s_max - e(0:q_max) / t
     boltz = EXP ( boltz )

     ! Calculate averages
     norm  = SUM ( boltz )
     g_avg = SUM ( boltz * g(0:q_max) ) / norm
     e_avg = SUM ( boltz * e(0:q_max) ) / norm
     e_msd = SUM ( boltz * (e(0:q_max) - e_avg )**2  ) / norm
     e_avg = e_avg / REAL(n)            ! Energy per atom
     e_msd = e_msd / ( REAL(n) * t**2 ) ! Heat capacity per atom

     WRITE ( unit=output_unit, fmt='(4f15.6)' ) t, g_avg, e_avg, e_msd

  END DO

  DEALLOCATE ( e, h, g, s, boltz )

END PROGRAM wl_hist
