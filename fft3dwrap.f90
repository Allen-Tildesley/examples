! fft3dwrap.f90
! 3D fast Fourier transform applied to a Gaussian function
PROGRAM fft3dwrap
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

  ! The program calls use the C subroutine library FFTW to perform the finite Fourier transforms
  ! The details of this library are available at http://www.fftw.org/
  ! We assume that compiler flags are set such that real and integer fortran variables
  ! have the appropriate precision to match their C counterparts


  INTEGER(C_INT), PARAMETER :: sc   = 128   ! the number of points in real space in each dimension
  INTEGER,        PARAMETER :: sc2  = sc/2  ! In this example sc is the same in each dimension
  REAL,           PARAMETER :: lper = 6.0   ! The periodic repeat distance in each dimension

  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(sc,sc,sc) :: fin, fout ! Contains the data to be transformed and the output matrix

  TYPE(C_PTR)      :: plan
  INTEGER          :: i, j, k, nup
  REAL             :: delta, delta3, x, y, z, g
  REAL             :: dk, kx, ky, kz, ksq, kmag, ghat
  REAL, PARAMETER  :: pi = 4.0 * ATAN( 1.0 ) ! decay parameter of Gaussian is set to pi

  delta      = lper / REAL (sc)     ! periodic repeat distance in x, y and z
  delta3     = delta ** 3           ! factor for forward transform


  WRITE(*,*) 'Input Gaussian'
  DO i = 1, sc

     IF ( i  <= sc2 ) THEN
        x  = REAL( i - 1 ) * delta
     ELSE
        x  = REAL( i - sc - 1 ) * delta
     END IF

     DO j = 1, sc

        IF ( j <= sc2 ) THEN
           y  = REAL( j - 1 ) * delta
        ELSE
           y  = REAL( j - sc - 1 ) * delta
        END IF

        DO k = 1, sc

           IF ( k <= sc2 ) THEN
              z = REAL( k - 1 ) * delta
           ELSE
              z = REAL( k - sc - 1 ) * delta
           END IF

           g          = EXP ( - pi * (x * x + y * y + z * z )) ! setup 3d Gaussian
           fin(i,j,k) = CMPLX ( g, 0.0 )                       ! feed into complex array for FFT

           !write some elements of data
           nup = FLOOR( SQRT (REAL( i*i +j*j + k*k ) ) )
           IF( nup <= 3 ) THEN
              WRITE(*,'(3I4, 2F10.4)') i, j, k, fin(i,j,k)
           END IF

        END DO
     END DO
  END DO
  WRITE(*,'(/)')

  ! Forward FFT

  plan = fftw_plan_dft_3d ( sc, sc, sc, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE) ! set up plan for the Fourier transform
  CALL fftw_execute_dft ( plan, fin, fout )                                     ! execute FFT
  CALL fftw_destroy_plan ( plan )                                               ! release plan

  dk = (2.0 * pi) / delta / REAL(sc) ! interval in reciprocal space

  WRITE (*, *) 'Reciprocal-space transform'

  DO i = 1, sc

     IF ( i <= sc2 ) THEN
        kx  = REAL( i - 1 ) * dk
     ELSE
        kx  = REAL( i - sc - 1 ) * dk
     END IF

     DO j = 1, sc

        IF ( j <= sc2 ) THEN
           ky = REAL( j - 1 ) * dk
        ELSE
           ky = REAL( j - sc - 1 ) * dk
        END IF

        DO k = 1, sc

           IF( k <= sc2 ) THEN
              kz = REAL( k - 1 ) * dk
           ELSE
              kz = REAL( k - sc - 1 )  * dk
           END IF

           ksq  = kx * kx + ky * ky + kz * kz  ! square of wave vector
           kmag = SQRT( ksq )                  ! modulus of wave vector
           ghat = EXP(- ksq / 4 / pi )         ! analytical transform of the Gaussian

           ! Write some elements of data in reciprocal space including factor of delta**3
           ! Compare with the analytical expression for the transform of the Gaussian test function

           nup  = FLOOR( SQRT (REAL( i*i +j*j + k*k ) ) )
           IF ( nup <= 3 ) THEN
              WRITE(*,'(3I4, 4F10.4)') i, j, k, kmag, delta3 * fout(i,j,k), ghat
           END IF
        END DO
     END DO
  END DO
  WRITE(*,'(/)')

  ! backward Fourier transform

  plan = fftw_plan_dft_3d ( sc, sc, sc, fout, fin, FFTW_BACKWARD, FFTW_ESTIMATE) ! set up plan for the Fourier transform
  CALL fftw_execute_dft ( plan, fout, fin )                                      ! execute FFT
  CALL fftw_destroy_plan ( plan )                                                ! release plan

  ! Write some elements of data in real space after the back transform including the normalising
  ! factor 1/n**3. These should be compared with the input data

  WRITE(*,*) 'Back Transform to real space'
  DO i = 1, sc
     DO j = 1, sc
        DO k = 1, sc
           nup = FLOOR( SQRT (REAL( i*i +j*j + k*k ) ) )
           IF( nup <= 3) THEN
              WRITE(*,'(3I4, 2F10.4)') i, j, k, fin(i,j,k)/REAL( sc*sc*sc  )
           END IF
        END DO
     END DO
  END DO

END PROGRAM fft3dwrap

