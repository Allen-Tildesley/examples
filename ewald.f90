! ewald.f90
! Test routines provided in ewald_module.f90
PROGRAM ewald

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

  ! Reads an atomic configuration with periodic boundary conditions from cnf.inp
  ! Calculates r-space and k-space contributions to potential energy
  ! for given screening parameter kappa and number of wave vectors determined by nk
  ! Adds the surface term.
  ! Compares with brute force summation in real space over shells of periodic boxes
  ! Leave namelist empty to accept supplied defaults

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms
  USE               ewald_module,     ONLY : pot_r_ewald, pot_k_ewald, pot_k_pm_ewald

  IMPLICIT NONE

  INTEGER                              :: n         ! Number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r         ! Positions (3,n)
  REAL,    DIMENSION(:),   ALLOCATABLE :: q         ! Charges (n)
  REAL,    DIMENSION(:),   ALLOCATABLE :: pot_shell ! Potential for each shell (0:nbox**2)

  CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp'
  REAL, DIMENSION(3)          :: rij, dipole, rbox_vec
  REAL                        :: box, kappa, pot, pot_r, pot_k, pot_s
  INTEGER                     :: ioerr
  INTEGER                     :: i, j, xbox, ybox, zbox, rbox, rbox_sq, nbox, nbox_sq, nk
  REAL,             PARAMETER :: pi = 4.0*ATAN(1.0), twopi = 2.0*pi

  NAMELIST /nml/ kappa, nk, nbox

  ! Sensible default values
  kappa = 6.0
  nk    = 8
  nbox  = 8

  READ ( unit=input_unit, nml=nml, iostat=ioerr ) ! Namelist input
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in ewald'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'kappa*box',                           kappa
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Wave-vector limit in each direction', nk
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Brute force box limit',               nbox
  nbox_sq = nbox**2

  CALL read_cnf_atoms ( filename, n, box ) ! First call to obtain n
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'  ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Box (in sigma units)', box

  ALLOCATE ( r(3,n), q(n), pot_shell(0:nbox_sq) )

  CALL read_cnf_atoms ( filename, n, box, r ) ! Second call to read in configuration
  r = r / box         ! Work throughout in unit box
  r = r - ANINT ( r ) ! Apply periodic boundaries

  ! Assign +/- charges ensuring net charge is zero
  IF ( MOD(n,2) /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Error, we require even N for charge balance'
     STOP 'Error in ewald'
  END IF
  q(1:n:2) =  1.0
  q(2:n:2) = -1.0
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Net charge', SUM(q)

  ! Compute potential energy by Ewald method
  pot_r  = pot_r_ewald ( n, r, q, kappa )                ! Real-space term involving screened Coulomb potential
  pot_k  = pot_k_ewald ( nk, n, r, q, kappa )            ! Reciprocal space term
  dipole = SUM ( SPREAD(q(:),dim=1,ncopies=3)*r, dim=2 ) ! Calculate overall box dipole
  pot_s  = ( twopi / 3.0 ) * SUM ( dipole**2 )           ! Surface term
  pot    = pot_r + pot_k + pot_s                         ! Total potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'r-space potential energy', pot_r
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'k-space potential energy', pot_k
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'surface potential energy', pot_s
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'total   potential energy', pot

  ! Compare with an illustrative (simplified) particle-mesh method
  pot_k  = pot_k_pm_ewald ( nk, n, r, q, kappa ) ! Reciprocal space term
  pot    = pot_r + pot_k + pot_s                 ! Total potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'k-space potential energy (PME)', pot_k
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'total   potential energy (PME)', pot

  ! Compare with brute force calculation
  ! Big multiple loop over all pairs and surrounding periodic boxes
  ! For clarity we count all pairs twice, as ij and ji
  ! Potentials are stored according to squared distance of neighbouring box
  WRITE ( unit=output_unit, fmt='(a)' ) 'Brute force method'

  pot_shell = 0.0 ! Zero array of potentials (not all the elements of this array will be filled)

  DO xbox = -nbox, nbox
     DO ybox = -nbox, nbox
        DO zbox = -nbox, nbox
           rbox_sq = xbox**2 + ybox**2 + zbox**2
           IF ( rbox_sq > nbox_sq ) CYCLE ! Skip if outside maximum sphere of boxes
           rbox_vec = REAL ( [xbox,ybox,zbox] )
           DO i = 1, n
              DO j = 1, n
                 IF ( i==j .AND. rbox_sq==0 ) CYCLE ! Skip only for i=j in central box
                 ! Bare Coulomb term, including box vector, no periodic box correction
                 rij = r(:,i) - r(:,j) - rbox_vec
                 pot_shell(rbox_sq) = pot_shell(rbox_sq) + q(i)*q(j) / SQRT ( SUM ( rij**2 ) )
              END DO
           END DO
        END DO
     END DO
  END DO

  ! Correct for double counting
  pot_shell = pot_shell / 2.0

  ! Convert to cumulative sum
  DO rbox_sq = 1, nbox_sq
     pot_shell(rbox_sq) = pot_shell(rbox_sq) + pot_shell(rbox_sq-1)
  END DO

  ! Write out results for increasing spherical cutoff
  WRITE ( unit=output_unit, fmt='(a)' ) 'Shell      Potential'
  DO rbox = 0, nbox
     WRITE ( unit=output_unit, fmt='(i5,f15.6)' ) rbox, pot_shell(rbox**2)
  END DO

  DEALLOCATE ( r, q, pot_shell )

END PROGRAM ewald



