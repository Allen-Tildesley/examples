! grint.f90
! ! g(z,c,r) in a planar interface
PROGRAM grint
  !
  ! TODO DJT to complete the code

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

  ! Reads a trajectory from a sequence of configuration files
  ! Calculates pair distribution function for a planar interface in the xy plane,
  ! including dependence on z and symmetry breaking with respect to z direction

  ! Single-particle density profile in box coordinates is written to a file 'den.out'
  ! The combined profile for both interfaces, relative to interface position, is written to rho1.out
  ! More work is needed to check the two-particle density calculation, its conversion to the appropriate
  ! distribution function, and output to a file rho2.out. DJT to complete this!

  ! For illustration and simplicity, we adopt a scheme of formatted files of the same kind
  ! as those that are saved at the end of each block of our MD simulation examples
  ! We assume that the initial configuration of a run has been copied to cnf.000
  ! and subsequent configurations are called cnf.001 cnf.002 etc., up to (at most) cnf.999
  ! Obviously, in a practical application, a binary trajectory file would fulfil this role.

  ! Cubic periodic boundary conditions are assumed
  ! r and box are assumed to be in the same units (e.g. LJ sigma)
  ! box is assumed to be unchanged throughout

  ! Values of basic parameters are read from standard input using a namelist nml

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms
  USE               grint_module,     ONLY : fit, nterms

  IMPLICIT NONE

  ! Most important variables
  INTEGER :: n      ! Number of atoms
  INTEGER :: t      ! Time (in blocks, equivalent to number of file)
  REAL    :: box    ! Box length
  REAL    :: dz     ! Spacing in z-direction
  REAL    :: dz_box ! Spacing in z-direction adjusted to fit box
  REAL    :: dr     ! Spacing in r_ij
  REAL    :: dc     ! Spacing in cos(theta_ij)
  REAL    :: z_mid  ! Rough estimate of mid-point of liquid slab
  INTEGER :: nz     ! Number of z points relative to interface location
  INTEGER :: nr     ! Number of r points in pair density
  INTEGER :: nc     ! Number of cos(theta) points in pair density
  INTEGER :: nk     ! Number of z points in simulation box

  REAL, DIMENSION(:,:),   ALLOCATABLE :: r    ! Positions (3,n)
  REAL, DIMENSION(:),     ALLOCATABLE :: rz   ! Saved z-coordinates
  REAL, DIMENSION(:),     ALLOCATABLE :: d    ! Single-particle density profile snapshot (nk)
  REAL, DIMENSION(:),     ALLOCATABLE :: dens ! Single-particle density profile average (nk)
  REAL, DIMENSION(:),     ALLOCATABLE :: z    ! Positions for density profile (nk)
  REAL, DIMENSION(:),     ALLOCATABLE :: rho1 ! One-particle density relative to interface location (-nz:nz)
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: rho2 ! Two-particle density relative to interface location (-nz:nz,nc,nr)
  REAL, DIMENSION(nterms)             :: cfit ! Coefficients for fit

  CHARACTER(len=4), PARAMETER :: cnf_prefix  = 'cnf.'
  CHARACTER(len=4), PARAMETER :: den_prefix  = 'den.'
  CHARACTER(len=5), PARAMETER :: rho1_prefix = 'rho1.'
  CHARACTER(len=5), PARAMETER :: rho2_prefix = 'rho2.'
  CHARACTER(len=3)            :: sav_tag
  INTEGER                     :: i, k, unit, ioerr
  LOGICAL                     :: exists, fail
  REAL                        :: norm, area
  REAL, PARAMETER             :: pi = 4.0*ATAN(1.0), twopi = 2.0*pi

  NAMELIST /nml/ dz, dr, z_mid, nz, nc, nr

  ! Example default values
  dz    = 0.1 ! Spacing in z-direction
  nz    = 20  ! Number of z-points relative to interface location
  nc    = 20  ! Number of cos(theta) points covering range (-1,1)
  nr    = 20  ! Number of r points
  dr    = 0.1 ! r-spacing
  z_mid = 0.0 ! First estimate of location of liquid slab midpoint

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in grint'
  END IF
  dc = 2.0 / REAL(nc) ! Cosine spacing to cover the range (-1,1)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Spacing in z-direction',      dz
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Spacing in cos(theta)',       dc
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Spacing in r',                dr
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of z points',          nz
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of cos(theta) points', nc
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of r points',          nr
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Liquid slab midpoint',        z_mid

  sav_tag = '000' ! Use initial configuration to get N
  INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! Check the file exists
  IF ( .NOT. exists ) THEN
     WRITE ( unit=error_unit, fmt='(a,a)') 'File does not exist: ', cnf_prefix//sav_tag
     STOP 'Error in grint'
  END IF
  CALL read_cnf_atoms ( cnf_prefix//sav_tag, n, box ) ! Read n and box

  ! We define dz_box to fit the box exactly
  nk     = NINT(box/dz)
  dz_box = box / REAL(nk)
  area   = box**2

  ALLOCATE ( r(3,n), rz(n), d(nk), dens(nk), z(nk), rho1(-nz:nz), rho2(-nz:nz,nc,nr) )

  z(:) = [ ( -0.5*box + (REAL(k)-0.5)*dz_box, k = 1, nk ) ]

  dens = 0.0
  rho1 = 0.0
  rho2 = 0.0
  norm = 0.0
  t    = 0

  ! Initial guesses at slab fit parameters
  ! These will be passed on at each step, assuming that changes are small
  cfit(1) = -0.25*box ! Gas-liquid interface position
  cfit(2) =  0.25*box ! Liquid-gas interface position
  cfit(3) = 1.0       ! Width of interface
  cfit(4) = 0.0       ! Gas density
  cfit(5) = 0.8       ! Liquid density

  DO ! Loop through configuration files

     IF ( t > 999 ) EXIT ! Our naming scheme only goes up to cnf.999

     WRITE ( sav_tag, fmt='(i3.3)' ) t                      ! Number of this configuration
     INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! Check the file exists
     IF ( .NOT. exists ) EXIT                               ! We have come to the end of the data

     CALL read_cnf_atoms ( cnf_prefix//sav_tag, n, box, r ) ! Read configuration
     r(3,:) = r(3,:) - z_mid              ! Place liquid slab approximately in centre of box
     r      = r - ANINT ( r / box ) * box ! Apply periodic boundary conditions

     d(:) = 0.0

     DO i = 1, n
        k = 1 + FLOOR ( ( r(3,i) + box/2.0 ) / dz_box )
        k = MAX ( 1,  k ) ! Guard against roundoff
        k = MIN ( nk, k ) ! Guard against roundoff
        d(k)    = d(k)    + 1.0
        dens(k) = dens(k) + 1.0
     END DO

     d(:) = d(:) / ( area * dz_box )

     ! Fit profile in order to obtain interface positions in cfit(1) and cfit(2)
     CALL fit ( z, d, cfit, fail )
     IF ( fail ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Failed in fit routine'
        STOP 'Error in grint'
     END IF
     
     rz(:) = r(3,:) ! Save z-coordinates of all atoms

     ! Process gas-liquid interface
     r(3,:) = rz(:) - cfit(1)                       ! Shift gas-liquid interface to origin
     r(3,:) = r(3,:) - ANINT ( r(3,:) / box ) * box ! Apply PBC
     CALL calculate

     ! Process liquid-gas interface
     r(3,:) = cfit(2) - rz(:)                       ! Shift liquid-gas interface to origin and reflect
     r(3,:) = r(3,:) - ANINT ( r(3,:) / box ) * box ! Apply PBC
     CALL calculate
     
     norm  = norm + 1.0
     t     = t + 1 ! Ready for next file
     z_mid = z_mid + 0.5*(cfit(1)+cfit(2)) ! Refine estimate of slab midpoint for next time

  END DO ! End loop through configuration files

  ! Normalize
  dens = dens / ( norm * area * dz_box )
  rho1 = rho1 / ( 2.0 * norm * area * dz )
  rho2 = rho2 / ( twopi * norm * area * dr * dc * dz )

  ! Output results
  
  OPEN ( newunit=unit, file=den_prefix//'out', status='replace', iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error opening file', ioerr
     STOP 'Error in grint'
  END IF
  DO k = 1, nk
     WRITE ( unit=unit, fmt='(2f15.6)' ) z(k), dens(k)
  END DO
  CLOSE(unit=unit)

  OPEN ( newunit=unit, file=rho1_prefix//'out', status='replace', iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error opening file', ioerr
     STOP 'Error in grint'
  END IF
  DO k = -nz, nz
     WRITE ( unit=unit, fmt='(2f15.6)' ) REAL(k)*dz, rho1(k)
  END DO
  CLOSE(unit=unit)

  DEALLOCATE ( r, rz, d, dens, rho1, rho2, z )

CONTAINS
  
  SUBROUTINE calculate

    ! This routine carries out the histogramming for rho1 and rho2

    IMPLICIT NONE
    INTEGER            :: i, j, ic, ir, iz
    REAL               :: rij_sq, rij_mag, c, r_cut, r_cut_sq
    REAL, DIMENSION(3) :: rij

    r_cut = REAL(nr)*dr
    r_cut_sq = r_cut ** 2

    DO i = 1, n ! Loop over atoms

       iz = NINT ( r(3,i) / dz )
       IF ( ABS(iz) > nz ) CYCLE ! No need to consider i-atoms outside range

       rho1(iz) = rho1(iz) + 1.0

       DO j = 1, n ! Loop over all other atoms

          IF ( i==j ) CYCLE ! Skip self

          rij    = r(:,i) - r(:,j)
          rij    = rij - ANINT ( rij / box ) * box ! Apply PBC
          rij_sq = SUM(rij**2) ! Magnitude of separation

          IF ( rij_sq > r_cut_sq ) CYCLE

          rij_mag = SQRT ( rij_sq )
          ir = 1 + FLOOR ( rij_mag/dr )
          ir = MAX(1, ir) ! Guard against roundoff
          ir = MIN(nr,ir) ! Guard against roundoff

          c  = rij(3) / rij_mag ! Cosine of angle
          ic = 1 + FLOOR ( (1.0+c)/dc )
          ic = MAX(1, ic) ! Guard against roundoff
          ic = MIN(nc,ic) ! Guard against roundoff

          rho2(iz,ic,ir) = rho2(iz,ic,ir) + 1.0 / rij_sq ! Incorporating r**2 factor here

       END DO ! End loop over all other atoms

    END DO ! End loop over atoms

  END SUBROUTINE calculate
  
END PROGRAM grint
