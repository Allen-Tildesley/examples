! pair_distribution.f90
! Calculates pair distribution function g(r)
PROGRAM pair_distribution

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

  ! Reads a trajectory from a sequence of configuration files
  ! For illustration and simplicity, we adopt a scheme of formatted files of the same kind
  ! as those that are saved at the end of each block of our MD simulation examples
  ! We assume that the initial configuration of a run has been copied to cnf.000
  ! and subsequent configurations are called cnf.001 cnf.002 etc., up to (at most) cnf.999
  ! Obviously, in a practical application, a binary trajectory file would fulfil this role.

  ! Cubic periodic boundary conditions are assumed
  ! r and box are assumed to be in the same units (e.g. LJ sigma)
  ! box is assumed to be unchanged throughout
  ! dr, the desired resolution of g(r), should be provided in the same units
  ! The entire calculation is performed in box=1 units, but the grid points are
  ! converted back to sigma=1 units before output.

  ! The value of dr is read from standard input using a namelist nml
  ! Leave namelist empty to accept supplied default

  ! Results are written to a file 'pair_distribution.out' with diagnostics to standard output

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms

  IMPLICIT NONE

  INTEGER :: n   ! number of atoms
  REAL    :: box ! box length (assumed constant throughout)
  REAL    :: dr  ! spacing in g(r)

  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,n)
  integer, DIMENSION(:),   ALLOCATABLE :: h ! histogram of separations (nk)
  REAL,    DIMENSION(:),   ALLOCATABLE :: g ! pair distribution function (nk)

  INTEGER                     :: i, j, k, nk, nstep
  REAL, DIMENSION(3)          :: rij
  REAL                        :: rij_sq, r_hi, r_lo, h_id, const, rho
  REAL,             PARAMETER :: pi = 4.0*ATAN(1.0)

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3)            :: sav_tag
  INTEGER                     :: unit, ioerr
  LOGICAL                     :: exists

  NAMELIST /nml/ dr

  ! Example default values
  dr = 0.02

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in pair_distribution'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'g(r) spacing dr = ', dr

  WRITE ( sav_tag, fmt='(i3.3)' ) 0 ! Use initial configuration to get n and box
  INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! Check the file exists
  IF ( .NOT. exists ) THEN
     WRITE ( unit=error_unit, fmt='(a,a)') 'File does not exist: ', cnf_prefix//sav_tag
     STOP 'Error in pair_distribution'
  END IF
  CALL read_cnf_atoms ( cnf_prefix//sav_tag, n, box ) ! Read n and box

  dr = dr / box         ! Convert to box=1 units
  nk = FLOOR ( 0.5/dr ) ! Accumulate out to half box length
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'number of bins = ', nk

  ALLOCATE ( r(3,n), h(nk), g(nk) )

  h(:) = 0 ! Initialize to zero
  nstep = 0

  DO ! Single sweep through data until end

     IF ( nstep >= 1000 ) EXIT ! Our naming scheme only goes up to cnf.999

     WRITE ( sav_tag, fmt='(i3.3)' ) nstep                  ! Number of configuration
     INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! Check the file exists
     IF ( .NOT. exists ) EXIT                               ! We have come to the end of the data

     CALL read_cnf_atoms ( cnf_prefix//sav_tag, n, box, r ) ! Read configuration
     WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Processing ', nstep
     r = r / box ! Convert to box=1 units

     ! The following is carried out for nstep configurations

     DO i = 1, n-1
        DO j = i+1, n
           rij(:) = r(:,i) - r(:,j)
           rij(:) = rij(:) - ANINT ( rij(:) )
           rij_sq = SUM ( rij**2 )
           k      = FLOOR ( SQRT ( rij_sq ) / dr ) + 1
           IF ( k <= nk ) h(k) = h(k) + 2
        END DO
     END DO

     nstep = nstep + 1 ! Increment step counter ready for next time
     
  END DO ! End single sweep through data until end

  rho = REAL(n) ! Our calculation is done in box=1 units
  const = 4.0 * pi * rho / 3.0
  DO k = 1, nk
     g(k) = REAL ( h(k) ) / REAL ( n * nstep ) ! Average number
     r_lo = REAL ( k - 1 ) * dr
     r_hi = r_lo + dr
     h_id = const * ( r_hi ** 3 - r_lo ** 3 ) ! Ideal number
     g(k) = g(k) / h_id
  END DO
  WRITE ( unit=output_unit, fmt='(a)' ) 'Output to pair_distribution.out'

  OPEN ( newunit=unit, file='pair_distribution.out', status='replace', iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error opening file', ioerr
     STOP 'Error in pair_distribution'
  END IF

  dr = dr*box ! Grid spacing in sigma=1 units; g(r) is dimensionless
  DO k = 1, nk
     WRITE ( unit=unit, fmt='(2f15.8)' ) (REAL(k)-0.5)*dr, g(k)
  END DO
  CLOSE(unit=unit)

END PROGRAM pair_distribution
