! fft3dwrap.f90
! 3D fast Fourier transform applied to a Gaussian function
PROGRAM fft3dwrap

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE, INTRINSIC :: iso_c_binding

  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

  ! The program calls use the C subroutine library FFTW to perform the finite Fourier transforms
  ! The details of this library are available at http://www.fftw.org/
  ! We assume that compiler flags are set such that real and integer fortran variables
  ! have the appropriate precision to match their C counterparts

  ! In this example the box lengths and numbers of grid points are the same in each dimension
  INTEGER :: sc2 ! half number of grid points
  REAL    :: box ! periodic repeat distance
  REAL    :: dr  ! grid spacing in real space
  REAL    :: dk  ! grid spacing in reciprocal space

  INTEGER            :: ix, iy, iz, wx, wy, wz, ioerr
  REAL, DIMENSION(3) :: r, k
  REAL               :: r_sq, k_sq, k_mag, g, ghat

  INTEGER(C_INT)                                           :: sc      ! number of points for FFT
  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:,:), ALLOCATABLE :: fft_inp ! Data to be transformed (0:sc-1,0:sc-1,0:sc-1)
  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:,:), ALLOCATABLE :: fft_out ! Output data (0:sc-1,0:sc-1,0:sc-1)
  TYPE(C_PTR)                                              :: fft_plan! Plan needed for FFTW

  REAL,    PARAMETER :: pi = 4.0 * ATAN( 1.0 ) ! decay parameter of Gaussian is set to pi
  integer, parameter :: out_max = 15

  NAMELIST /nml/ sc2, box

  ! Typical default values
  sc2 = 2**6 ! Not essential to be a power of 2, but usually more efficient
  box = 6.0  ! Large enough to accommodate the chosen 3D Gaussian

  ! Namelist from standard input
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in fft3dwrap'
  END IF

  sc = sc2 * 2
  dr = box / REAL (sc)
  dk = (2.0 * pi) / dr / REAL(sc) ! interval in reciprocal space
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of grid points in each dimension, sc = ', sc
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Periodic repeat length (box) = ',                box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Grid spacing in real space (dr) = ',             dr
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Grid spacing in reciprocal space (dk) = ',       dk

  ALLOCATE ( fft_inp(0:sc-1,0:sc-1,0:sc-1), fft_out(0:sc-1,0:sc-1,0:sc-1) )

  WRITE ( unit=output_unit, fmt='(a)' ) 'Initial Gaussian'

  DO ix = 0, sc-1
     DO iy = 0, sc-1
        DO iz = 0, sc-1

           wx   = wraparound ( ix )
           wy   = wraparound ( iy )
           wz   = wraparound ( iz )
           r    = REAL ( [ wx, wy, wz ] ) * dr
           r_sq = SUM ( r**2 )
           g    = EXP ( - pi * r_sq ) ! Setup 3d Gaussian

           fft_inp(ix,iy,iz) = CMPLX ( g, 0.0 ) ! Feed into complex array for FFT

           ! Write some elements of data
           IF ( ix**2 + iy**2 + iz**2 <= out_max ) THEN
              WRITE ( unit=output_unit, fmt='(3i5,2f15.5)' ) ix, iy, iz, fft_inp(ix,iy,iz)
           END IF

        END DO
     END DO
  END DO

  WRITE(*,'(/)')

  ! Forward FFT

  fft_plan = fftw_plan_dft_3d ( sc, sc, sc, fft_inp, fft_out, FFTW_FORWARD, FFTW_ESTIMATE) ! Set up plan for the Fourier transform
  CALL fftw_execute_dft ( fft_plan, fft_inp, fft_out )                                     ! Execute FFT
  CALL fftw_destroy_plan ( fft_plan )                                                      ! Release plan

  WRITE ( unit=output_unit, fmt='(a)' ) 'Reciprocal-space transform'

  DO ix = 0, sc-1
     DO iy = 0, sc-1
        DO iz = 0, sc-1

           wx    = wraparound ( ix )
           wy    = wraparound ( iy )
           wz    = wraparound ( iz )
           k     = REAL ( [ wx, wy, wz ] ) * dk
           k_sq  = SUM ( k**2 )
           k_mag = SQRT( k_sq )            ! modulus of wave vector
           ghat  = EXP(- k_sq / 4.0 / pi ) ! analytical transform of the Gaussian

           ! Write some elements of data in reciprocal space including factor of dr**3
           ! Compare with the analytical expression for the transform of the Gaussian test function

           IF ( ix**2 + iy**2 + iz**2 <= out_max ) THEN
              WRITE ( unit=output_unit, fmt='(3i5,4f15.5)') ix, iy, iz, k_mag, fft_out(ix,iy,iz)*dr**3, ghat
           END IF

        END DO
     END DO
  END DO

  WRITE ( unit=output_unit, fmt='(/)')

  ! Backward Fourier transform

  fft_plan = fftw_plan_dft_3d ( sc, sc, sc, fft_out, fft_inp, FFTW_BACKWARD, FFTW_ESTIMATE) ! set up plan for the Fourier transform
  CALL fftw_execute_dft ( fft_plan, fft_out, fft_inp )                                      ! execute FFT
  CALL fftw_destroy_plan ( fft_plan )                                                       ! release plan

  ! Write some elements of data in real space after the back transform
  ! including the normalising factor 1/sc**3
  ! These should be compared with the input data

  WRITE ( unit=output_unit, fmt='(a)' ) 'Back Transform to real space'
  DO ix = 0, sc-1
     DO iy = 0, sc-1
        DO iz = 0, sc-1
           IF ( ix**2 + iy**2 + iz**2 <= out_max ) THEN
              WRITE ( unit=output_unit, fmt='(3i4,2f10.4)' ) ix, iy, iz, fft_inp(ix,iy,iz)/REAL(sc**3)
           END IF
        END DO
     END DO
  END DO

  DEALLOCATE ( fft_inp, fft_out )

CONTAINS

  FUNCTION wraparound ( m ) RESULT ( w ) ! implements wrapping of indices
    INTEGER, INTENT(in) :: m ! index to be wrapped
    INTEGER             :: w ! result

    IF ( m < 0 .OR. m >= sc ) STOP 'This should never happen'

    IF ( m < sc2 ) THEN
       w = m
    ELSE
       w = m - sc
    END IF

  END FUNCTION wraparound

END PROGRAM fft3dwrap

