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
  INTEGER :: sc2 ! half the number of grid points
  REAL    :: box ! periodic repeat distance
  REAL    :: dr  ! grid spacing in real space
  REAL    :: dk  ! grid spacing in reciprocal space

  INTEGER            :: ix, iy, iz, ioerr
  REAL, DIMENSION(3) :: r, k
  REAL               :: r_sq, k_sq, g

  INTEGER(C_INT)                                           :: sc      ! number of points for FFT
  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:,:), ALLOCATABLE :: fft_inp ! Data to be transformed (0:sc-1,0:sc-1,0:sc-1)
  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:,:,:), ALLOCATABLE :: fft_out ! Output data (0:sc-1,0:sc-1,0:sc-1)
  TYPE(C_PTR)                                              :: fft_plan! Plan needed for FFTW

  REAL,    PARAMETER :: pi = 4.0 * ATAN( 1.0 )
  INTEGER, PARAMETER :: out_max = 15

  NAMELIST /nml/ sc2, box

  ! Typical default values
  sc2 = 2**6 ! Not essential to be a power of 2, but usually more efficient
  box = 6.0  ! Large enough to accommodate the chosen 3D Gaussian, for good comparison with analytical result

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

  WRITE ( unit=output_unit, fmt='(a)'    ) 'Initial real-space Gaussian'
  WRITE ( unit=output_unit, fmt='(5a15)' ) '   ix   iy   iz', '|r|', 'Gaussian(r)', 'FFT (real)', 'FFT (imag)'

  ! Triple loop over xyz grid points (uses wraparound indexing)
  DO ix = 0, sc-1
     r(1) = REAL ( wraparound ( ix ) ) * dr
     DO iy = 0, sc-1
        r(2) = REAL ( wraparound ( iy ) ) * dr
        DO iz = 0, sc-1
           r(3) = REAL ( wraparound ( iz ) ) * dr

           r_sq = SUM ( r**2 )        ! squared distance from origin
           g    = EXP ( - pi * r_sq ) ! Setup 3d Gaussian (decay parameter chosen to be pi)

           fft_inp(ix,iy,iz) = CMPLX ( g, 0.0 ) ! Feed into complex array for FFT

           ! Write some elements of data in same form as later output
           IF ( ix**2 + iy**2 + iz**2 <= out_max ) THEN
              WRITE ( unit=output_unit, fmt='(3i5,4f15.5)' ) ix, iy, iz, SQRT(r_sq), g, fft_inp(ix,iy,iz)
           END IF

        END DO
     END DO
  END DO
  ! End triple loop over xyz grid points

  WRITE(*,'(/)')

  ! Forward FFT

  fft_plan = fftw_plan_dft_3d ( sc, sc, sc, fft_inp, fft_out, FFTW_FORWARD, FFTW_ESTIMATE) ! Set up plan for the FFT
  CALL fftw_execute_dft  ( fft_plan, fft_inp, fft_out )                                    ! Execute FFT
  CALL fftw_destroy_plan ( fft_plan )                                                      ! Release plan

  WRITE ( unit=output_unit, fmt='(a)' ) 'Reciprocal-space transform'
  WRITE ( unit=output_unit, fmt='(5a15)' ) '   ix   iy   iz', '|k|', 'Gaussian(k)', 'FFT (real)', 'FFT (imag)'

  ! Triple loop over xyz grid points (uses wraparound indexing)
  DO ix = 0, sc-1
     k(1) = REAL ( wraparound ( ix ) ) * dk
     DO iy = 0, sc-1
        k(2) = REAL ( wraparound ( iy ) ) * dk
        DO iz = 0, sc-1
           k(3) = REAL ( wraparound ( iz ) ) * dk

           ! Write some elements of data in reciprocal space including factor of dr**3
           ! Compare with the (real) analytical expression for the transform of the Gaussian test function

           IF ( ix**2 + iy**2 + iz**2 <= out_max ) THEN
              k_sq = SUM ( k**2 )             ! squared magnitude of wave vector
              g    = EXP ( -k_sq / 4.0 / pi ) ! analytical transform of the Gaussian
              WRITE ( unit=output_unit, fmt='(3i5,4f15.5)') ix, iy, iz, SQRT(k_sq), g, fft_out(ix,iy,iz)*dr**3
           END IF

        END DO
     END DO
  END DO
  ! End triple loop over xyz grid points

  WRITE ( unit=output_unit, fmt='(/)')

  ! Backward Fourier transform

  fft_plan = fftw_plan_dft_3d ( sc, sc, sc, fft_out, fft_inp, FFTW_BACKWARD, FFTW_ESTIMATE) ! Set up plan for the FFT
  CALL fftw_execute_dft  ( fft_plan, fft_out, fft_inp )                                     ! Execute FFT
  CALL fftw_destroy_plan ( fft_plan )                                                       ! Release plan

  ! Write some elements of data in real space after the back transform including the normalising factor 1/sc**3
  ! Compare with the (real) input data

  WRITE ( unit=output_unit, fmt='(a)'    ) 'Back Transform to real space'
  WRITE ( unit=output_unit, fmt='(5a15)' ) '   ix   iy   iz', '|r|', 'Gaussian(r)', 'FFT (real)', 'FFT (imag)'

  ! Triple loop over xyz grid points (uses wraparound indexing)
  DO ix = 0, sc-1
     r(1) = REAL ( wraparound ( ix ) ) * dr
     DO iy = 0, sc-1
        r(2) = REAL ( wraparound ( iy ) ) * dr
        DO iz = 0, sc-1
           r(3) = REAL ( wraparound ( iz ) ) * dr

           IF ( ix**2 + iy**2 + iz**2 <= out_max ) THEN
              r_sq = SUM ( r**2 )        ! squared distance from origin
              g    = EXP ( - pi * r_sq ) ! Original 3d Gaussian
              WRITE ( unit=output_unit, fmt='(3i5,4f15.5)' ) ix, iy, iz, SQRT(r_sq), g, fft_inp(ix,iy,iz)/REAL(sc**3)
           END IF

        END DO
     END DO
  END DO
  ! End triple loop over xyz grid points

  DEALLOCATE ( fft_inp, fft_out )

CONTAINS

  FUNCTION wraparound ( i ) RESULT ( w ) ! implements wrapping of indices
    INTEGER, INTENT(in) :: i ! index to be wrapped
    INTEGER             :: w ! result

    IF ( i < 0 .OR. i >= sc ) STOP 'This should never happen'

    IF ( i < sc2 ) THEN
       w = i
    ELSE
       w = i - sc
    END IF

  END FUNCTION wraparound

END PROGRAM fft3dwrap

