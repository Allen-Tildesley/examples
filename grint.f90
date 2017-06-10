! grint.f90
! g(z,c,r) in a planar interface
PROGRAM grint

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
  ! Calculates pair distribution function for a planar interface in the xy plane,
  ! including dependence on z and symmetry breaking with respect to z direction

  ! Single-particle density profile in box coordinates is written to a file 'den.out'
  ! The combined profile for both interfaces, relative to interface position, is written to rho.out
  ! Slices through the pair distribution function g2(z,c,r) where z=z1 is the z-coordinate of atom 1
  ! and c=cos(theta) is the angle of the r12 vector, are written out to files as a function of r.
  ! Assuming that the liquid phase is more-or-less central in the box, the interfaces are combined
  ! in the analysis and oriented so that z<0 is towards the gas and z>0 is towards the liquid.
  ! The cosine is defined so that c<0 corresponds to z1<z2 and c>0 to z1>z2.

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
  USE               grint_module,     ONLY : fit

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

  REAL, DIMENSION(:,:),   ALLOCATABLE :: r     ! Positions (3,n)
  REAL, DIMENSION(:),     ALLOCATABLE :: rz    ! Saved z-coordinates (n)
  REAL, DIMENSION(:),     ALLOCATABLE :: d     ! Single-particle density profile snapshot (nk)
  REAL, DIMENSION(:),     ALLOCATABLE :: dens  ! Single-particle density profile average (nk)
  REAL, DIMENSION(:),     ALLOCATABLE :: z_box ! Positions for density profile (nk)
  REAL, DIMENSION(:),     ALLOCATABLE :: rho1  ! One-particle density relative to interface location (-nz:nz)
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: rho2  ! Two-particle density relative to interface location (-nz:nz,-nc:nc,nr)
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: g2    ! Pair distribution relative to interface location (-nz:nz,-nc:nc,nr)
  REAL, DIMENSION(:),     ALLOCATABLE :: z     ! Positions for rho1, rho2, g2 (-nz:nz)

  REAL, DIMENSION(4) :: c1tanh ! Coefficients for 1-tanh fit
  REAL, DIMENSION(5) :: c2tanh ! Coefficients for 2-tanh fit

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: den_prefix = 'den'
  CHARACTER(len=3), PARAMETER :: rho_prefix = 'rho'
  CHARACTER(len=3), PARAMETER :: g2_prefix  = 'g2_'
  CHARACTER(len=4), PARAMETER :: out_tag = '.out' ! NB with dot
  CHARACTER(len=3)            :: sav_tag, z_tag, c_tag
  INTEGER                     :: i, k, unit, ioerr, zskip, cskip, iz, ir, ic, ic_max, iz_max
  LOGICAL                     :: exists, fail
  REAL                        :: norm, area, rij_mag, c, z1, z2, rho1_z1, rho1_z2
  REAL, PARAMETER             :: pi = 4.0*ATAN(1.0), twopi = 2.0*pi

  NAMELIST /nml/ dz, dr, z_mid, nz, nc, nr, zskip, cskip, iz_max

  ! Example default values
  dz    = 0.2  ! Spacing in z-direction
  nz    = 15   ! Number of z-points relative to interface location (-nz:+nz)
  nc    = 6    ! Number of cos(theta) points covering range (-1,1) (-nc:+nc)
  dr    = 0.02 ! r-spacing
  z_mid = 0.0  ! First estimate of location of liquid slab midpoint
  iz_max = 10  ! Limit on z values for g2 output
  zskip = 5    ! With iz_max=10 will give output files at iz = -10, -5, 0, +5, +10
  cskip = 3    ! With nc=6  will give output files at ic = -6,  -3, 0, +3, +6
  nr    = 200  ! Number of r-points

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in grint'
  END IF
  dc = 2.0 / REAL(2*nc+1) ! Cosine spacing to cover the range (-1,1)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Spacing in z-direction',       dz
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of z points',           nz
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) '+/- zmax',                     nz*dz
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Spacing in cos(theta)',        dc
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of cos(theta) points',  nc
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Spacing in r',                 dr
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of r points',           nr
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'rmax',                         nr*dr
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Output z skip',                zskip
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Output z max',                 iz_max
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Output c skip',                cskip
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Liquid slab midpoint (guess)', z_mid

  sav_tag = '000' ! Use initial configuration to get N
  INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! Check the file exists
  IF ( .NOT. exists ) THEN
     WRITE ( unit=error_unit, fmt='(a,a)') 'File does not exist: ', cnf_prefix//sav_tag
     STOP 'Error in grint'
  END IF
  CALL read_cnf_atoms ( cnf_prefix//sav_tag, n, box ) ! Read n and box

  ! We must remember that artefacts are expected whenever z approaches the "other" interface
  ! This depends on widths of the two phases, on nz*dz, and on nr*dr
  ! It's up to you if you ignore this warning
  IF ( ( nz*dz + nr*dr ) > 0.25*box ) THEN
     WRITE ( unit=output_unit, fmt='(a,t40,2f15.6)' ) 'Warning: max z > box/4 = ', (nz*dz+nr*dr), 0.25*box
  END IF

  ! We define dz_box to fit the box exactly
  nk     = NINT(box/dz)
  dz_box = box / REAL(nk)
  area   = box**2

  ALLOCATE ( r(3,n), rz(n), d(nk), dens(nk), z_box(nk) )
  ALLOCATE ( rho1(-nz:nz), rho2(-nz:nz,-nc:nc,nr), g2(-nz:nz,-nc:nc,nr), z(-nz:nz) )

  z_box(:) = [ ( -0.5*box + (REAL(k)-0.5)*dz_box, k = 1, nk ) ] ! Coordinates in simulation box
  z(:)     = [ ( REAL(k)*dz, k = -nz, nz ) ]                    ! Coordinates around interface

  dens = 0.0
  rho1 = 0.0
  rho2 = 0.0
  norm = 0.0
  t    = 0

  ! Initial guesses at slab fit parameters
  ! These will be passed on at each step, assuming that changes are small
  c2tanh(1) = -0.25*box ! Gas-liquid interface position
  c2tanh(2) =  0.25*box ! Liquid-gas interface position
  c2tanh(3) = 1.0       ! Width of interface
  c2tanh(4) = 0.0       ! Gas density
  c2tanh(5) = 0.8       ! Liquid density

  DO ! Loop through configuration files

     IF ( t > 999 ) EXIT ! Our naming scheme only goes up to cnf.999

     WRITE ( sav_tag, fmt='(i3.3)' ) t                      ! Number of this configuration
     INQUIRE ( file = cnf_prefix//sav_tag, exist = exists ) ! Check the file exists
     IF ( .NOT. exists ) EXIT                               ! We have come to the end of the data

     IF ( MOD(t,10)==0 ) WRITE ( unit=output_unit, fmt='(a)' ) 'Processing file '//cnf_prefix//sav_tag

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

     ! Fit profile in order to obtain interface positions in c2tanh(1) and c2tanh(2)
     CALL fit ( z_box, d, c2tanh, f2tanh, d2tanh, fail )
     IF ( fail ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Failed in fit routine'
        STOP 'Error in grint'
     END IF

     rz(:) = r(3,:) ! Save z-coordinates of all atoms

     ! Process gas-liquid interface
     r(3,:) = rz(:) - c2tanh(1)                       ! Shift gas-liquid interface to origin
     r(3,:) = r(3,:) - ANINT ( r(3,:) / box ) * box ! Apply PBC
     CALL calculate

     ! Process liquid-gas interface
     r(3,:) = c2tanh(2) - rz(:)                       ! Shift liquid-gas interface to origin and reflect
     r(3,:) = r(3,:) - ANINT ( r(3,:) / box ) * box ! Apply PBC
     CALL calculate

     norm  = norm + 1.0
     t     = t + 1 ! Ready for next file
     z_mid = z_mid + 0.5*(c2tanh(1)+c2tanh(2)) ! Refine estimate of slab midpoint for next time

  END DO ! End loop through configuration files

  ! Normalize (including factor for 2 interfaces)
  dens = dens / ( norm * area * dz_box )
  rho1 = rho1 / ( 2.0 * norm * area * dz )
  rho2 = rho2 / ( 2.0 * twopi * norm * area * dr * dc * dz )

  ! Fit the averaged density profile
  CALL fit ( z_box, dens, c2tanh, f2tanh, d2tanh, fail )

  ! Fit the single particle density
  c1tanh(1) = 0.0       ! Position of interface
  c1tanh(2) = c2tanh(3) ! Width of interface
  c1tanh(3) = c2tanh(4) ! Gas density
  c1tanh(4) = c2tanh(5) ! Liquid density
  CALL fit ( z, rho1, c1tanh, f1tanh, d1tanh, fail )

  ! Convert rho2 to g2, normalizing by fitted single-particle densities at z1 and z2

  DO iz = -nz, nz ! Loop over z coordinate

     z1  = z(iz)
     rho1_z1 = f1tanh ( z1, c1tanh ) ! Use fitted single-particle density (an approximation)

     DO  ic = -nc, nc ! Loop over cos(theta)

        c = REAL(ic) * dc

        DO ir = 1, nr ! Loop over radial distance

           rij_mag = ( REAL( ir ) - 0.5 ) * dr
           z2      = z1 - c * rij_mag
           rho1_z2 = f1tanh ( z2, c1tanh ) ! Use fitted single-particle density (an approximation)

           g2(iz,ic,ir) = rho2(iz,ic,ir) / ( rho1_z1 * rho1_z2  )

        END DO ! End loop over radial distance

     END DO ! End loop over cos(theta)

  END DO ! End loop over z coordinate

  ! Output results

  WRITE ( unit=output_unit, fmt='(a)' ) 'Box average density profile output to ' // den_prefix//out_tag

  OPEN ( newunit=unit, file=den_prefix//out_tag, status='replace', iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error opening dens file', ioerr
     STOP 'Error in grint'
  END IF
  DO k = 1, nk
     WRITE ( unit=unit, fmt='(3f15.6)' ) z_box(k), dens(k), f2tanh(z_box(k),c2tanh)
  END DO
  CLOSE(unit=unit)

  WRITE ( unit=output_unit, fmt='(a)' ) 'Single-particle density profile rho1 output to ' // rho_prefix//out_tag

  OPEN ( newunit=unit, file=rho_prefix//out_tag, status='replace', iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error opening rho file', ioerr
     STOP 'Error in grint'
  END IF
  DO k = -nz, nz
     WRITE ( unit=unit, fmt='(3f15.6)' ) z(k), rho1(k), f1tanh(z(k),c1tanh)
  END DO
  CLOSE(unit=unit)

  WRITE ( unit=output_unit, fmt='(a)' ) 'Pair distribution function g2 output in selected slices'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Each slice has fixed z=z1 and c=cos(theta)'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Filenames have the form g2_ZZZ_CCC'//out_tag

  ! Slices for selected values of iz, ic

  iz_max = MIN ( iz_max, nz )
  iz_max = iz_max / zskip
  iz_max = iz_max * zskip
  IF ( iz_max>99 ) STOP 'The output filename format will only cope with iz_max<100'
  WRITE ( unit=output_unit, fmt='(a)' ) 'ZZZ             z1'
  DO iz = -iz_max, iz_max, zskip
     WRITE ( unit=output_unit, fmt='(sp,i3.2,f15.5)' ) iz, z(iz)
  END DO
  WRITE ( unit=output_unit, fmt='(a)' ) '-ve sign means z1 on gas side, +ve sign means z1 on liquid side'

  ic_max = nc / cskip
  ic_max = ic_max * cskip
  IF ( ic_max>99 ) STOP 'The output filename format will only cope with ic_max<100'
  WRITE ( unit=output_unit, fmt='(a)' ) 'CCC     cos(theta)'
  DO ic = -ic_max, ic_max, cskip
     WRITE ( unit=output_unit, fmt='(sp,i3.2,f15.5)' ) ic, REAL(ic) * dc
  END DO
  WRITE ( unit=output_unit, fmt='(a)' ) '-ve sign means z1<z2, +ve sign means z1>z2'

  DO iz = -iz_max, iz_max, zskip
     WRITE ( z_tag, fmt='(sp,i3.2)' ) iz ! encoded with +/- sign
     DO ic = -ic_max, ic_max, cskip
        WRITE ( c_tag, fmt='(sp,i3.2)' ) ic ! encode with +/- sign
        OPEN ( newunit=unit, file=g2_prefix//z_tag//'_'//c_tag//out_tag, status='replace', iostat=ioerr )
        IF ( ioerr /= 0 ) THEN
           WRITE ( unit=error_unit, fmt='(a,i15)') 'Error opening g2 file', ioerr
           STOP 'Error in grint'
        END IF
        DO ir = 1, nr ! Loop over radial distance
           rij_mag = ( REAL( ir ) - 0.5 ) * dr
           WRITE ( unit=unit, fmt='(2f15.6)' ) rij_mag, g2(iz,ic,ir)
        END DO
        CLOSE(unit=unit)
     END DO
  END DO

  DEALLOCATE ( r, rz, d, dens, rho1, rho2, g2, z_box, z )

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
          rij_sq = SUM(rij**2) ! Squared magnitude of separation

          IF ( rij_sq > r_cut_sq ) CYCLE ! No need to consider separations out of range

          rij_mag = SQRT ( rij_sq ) ! Magnitude of separation
          ir = 1 + FLOOR ( rij_mag/dr )
          ir = MAX (  1, ir ) ! Guard against roundoff
          ir = MIN ( nr, ir ) ! Guard against roundoff

          c  = rij(3) / rij_mag ! Cosine of angle
          ic = NINT ( c / dc )
          ic = MAX ( -nc, ic ) ! Guard against roundoff
          ic = MIN (  nc, ic ) ! Guard against roundoff

          rho2(iz,ic,ir) = rho2(iz,ic,ir) + 1.0 / rij_sq ! Incorporating r**2 factor here

       END DO ! End loop over all other atoms

    END DO ! End loop over atoms

  END SUBROUTINE calculate

  FUNCTION f1tanh ( x, c ) RESULT ( f )
    IMPLICIT NONE
    REAL                           :: f ! Returns fitting function
    REAL,               INTENT(in) :: x ! Abscissa
    REAL, DIMENSION(:), INTENT(in) :: c ! Coefficients

    REAL :: t

    ! 1-tanh fit function (function value)
    IF ( SIZE(c) /= 4 ) STOP 'Error in function f1tanh'

    t = TANH ( ( x - c(1) ) / c(2) )

    f = 0.5 * ( c(4) + c(3) ) + 0.5 * ( c(4) - c(3) ) * t

  END FUNCTION f1tanh

  FUNCTION d1tanh ( x, c ) RESULT ( d )
    IMPLICIT NONE
    REAL,                    INTENT(in) :: x ! Abscissa
    REAL, DIMENSION(:),      INTENT(in) :: c ! Coefficients
    REAL, DIMENSION(SIZE(c))            :: d ! Returns fitting function derivatives

    REAL :: t

    ! 1-tanh fit function (derivatives)
    IF ( SIZE(c) /= 4 ) STOP 'Error in function d1tanh'

    t = TANH ( ( x - c(1) ) / c(2) )

    d(1) = 0.5 * ( c(4) - c(3) ) * ( t**2 - 1.0 ) / c(2)
    d(2) = 0.5 * ( c(4) - c(3) ) * ( x - c(1) ) * ( t**2 - 1.0 ) / c(2)**2
    d(3) = 0.5 - 0.5 * t
    d(4) = 0.5 + 0.5 * t

  END FUNCTION d1tanh

  FUNCTION f2tanh ( x, c ) RESULT ( f )
    IMPLICIT NONE
    REAL                           :: f ! Returns fitting function
    REAL,               INTENT(in) :: x ! Abscissa
    REAL, DIMENSION(:), INTENT(in) :: c ! Coefficients

    REAL :: t1, t2

    ! 2-tanh fit function (function value)
    IF ( SIZE(c) /= 5 ) STOP 'Error in function f2tanh'

    t1 = TANH ( ( x - c(1) ) / c(3) )
    t2 = TANH ( ( x - c(2) ) / c(3) )

    f = c(4) + 0.5 * ( c(5) - c(4) ) * ( t1 - t2 )

  END FUNCTION f2tanh

  FUNCTION d2tanh ( x, c ) RESULT ( d )
    IMPLICIT NONE
    REAL,                    INTENT(in) :: x ! Abscissa
    REAL, DIMENSION(:),      INTENT(in) :: c ! Coefficients
    REAL, DIMENSION(SIZE(c))            :: d ! Returns fitting function derivatives

    REAL :: t1, t2

    ! 2-tanh fit function (derivatives)
    IF ( SIZE(c) /= 5 ) STOP 'Error in function d2tanh'

    t1 = TANH ( ( x - c(1) ) / c(3) )
    t2 = TANH ( ( x - c(2) ) / c(3) )

    d(1) = 0.5 * ( c(5) - c(4) ) * ( t1**2 - 1.0 ) / c(3)
    d(2) = 0.5 * ( c(5) - c(4) ) * ( 1.0 - t2**2 ) / c(3)
    d(3) = 0.5 * ( c(5) - c(4) ) * ( ( x - c(1) ) * ( t1**2 - 1.0 )  + ( x - c(2) ) * ( 1.0 - t2**2 ) ) / c(3)**2
    d(4) = 1.0 - 0.5 * ( t1 - t2 )
    d(5) = 0.5 * ( t1 - t2 )

  END FUNCTION d2tanh

END PROGRAM grint
