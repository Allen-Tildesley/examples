! grint_module.f90
! Routines and data for grint program
MODULE grint_module

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

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routine
  PUBLIC :: fit

CONTAINS

  SUBROUTINE fit ( x, y, c, func, func_derivs, fail )
    IMPLICIT NONE
    REAL,    DIMENSION(:), INTENT(in)    :: x    ! Abscissae (npts)
    REAL,    DIMENSION(:), INTENT(in)    :: y    ! Ordinates (npts)
    REAL,    DIMENSION(:), INTENT(inout) :: c    ! Coefficients in fit (nterms)
    LOGICAL,               INTENT(out)   :: fail ! Indicates success or failure of fit

    INTERFACE
       FUNCTION func ( x, c ) RESULT ( f )
         IMPLICIT NONE
         REAL                           :: f ! Returns fitting function
         REAL,               INTENT(in) :: x ! Abscissa
         REAL, DIMENSION(:), INTENT(in) :: c ! Coefficients
       END FUNCTION func
       FUNCTION func_derivs ( x, c ) RESULT ( d )
         IMPLICIT NONE
         REAL,                    INTENT(in) :: x ! Abscissa
         REAL, DIMENSION(:),      INTENT(in) :: c ! Coefficients
         REAL, DIMENSION(SIZE(c))            :: d ! Returns fitting function derivatives
       END FUNCTION func_derivs
    END INTERFACE

    ! This fitting routine traces its origins back to an early edition of
    ! "Data reduction and error analysis for the physical sciences" by PR Bevington,
    ! since when both the program and the book have evolved significantly
    ! Initial values for the coefficient array c must be provided.

    REAL, DIMENSION(SIZE(c))         :: c_new, beta, dy, sigma
    REAL, DIMENSION(SIZE(c),SIZE(c)) :: alpha, array
    REAL, DIMENSION(SIZE(x))         :: yfit

    INTEGER :: npts, nterms, nfree, i, j, k, t, iter
    REAL    :: lambda, chisq, chisq_red, chisq_old, change

    REAL,    PARAMETER :: tol = 1.e-6
    LOGICAL, PARAMETER :: verbose = .FALSE. ! Controls output

    npts = SIZE(x) ! Number of points
    IF ( SIZE(y) /= npts ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i5)') 'Array dimensioning error', npts, SIZE(y)
       STOP 'Error in fit'
    END IF

    nterms = SIZE(c) ! Number of terms in fitting function

    nfree = npts - nterms
    IF (nfree <= 0) THEN
       WRITE ( unit=error_unit, fmt='(a,3i5)') 'Too few degrees of freedom', npts, nterms, nfree
       STOP 'Error in fit'
    END IF

    iter = 0
    DO i = 1,npts
       yfit(i) = func(x(i),c)
    END DO
    chisq     = SUM ( (y-yfit)**2 )
    chisq_red = chisq / REAL(nfree)

    ! Write headings and initial values
    IF ( verbose ) THEN
       WRITE ( unit=output_unit, fmt='(a5)', advance='no' ) 'Iter'
       DO t = 1, nterms
          WRITE ( unit=output_unit, fmt='(a11,i1,a10,i1,a1)', advance='no' ) 'c', t, 'sigma(c', t, ')'
       END DO
       WRITE ( unit=output_unit, fmt='(3a12)' ) 'Red chisq', 'change', 'lambda'
       WRITE ( unit=output_unit, fmt='(i5)', advance='no' ) iter
       DO t = 1, nterms
          WRITE ( unit=output_unit, fmt='(f12.6,12x)', advance='no' ) c(t)
       END DO
       WRITE ( unit=output_unit, fmt='(es12.2)' ) chisq_red
    END IF

    ! Carry out fit
    lambda = 0.001

    DO ! Loop until change is below tol or failure to converge 

       iter = iter + 1
       chisq_old = chisq_red

       beta  = 0.0
       alpha = 0.0

       DO i = 1, npts
          dy = func_derivs(x(i),c)
          DO j = 1, nterms
             beta(j) = beta(j) + dy(j) * ( y(i) - yfit(i) )
             DO k = 1, j
                alpha(j,k) = alpha(j,k) + dy(j) * dy(k)
                alpha(k,j) = alpha(j,k)
             END DO
          END DO
       END DO

       DO j = 1, nterms
          DO k = 1, nterms
             array(j,k) = alpha(j,k) / SQRT( alpha(j,j) * alpha(k,k) )
          END DO
          array(j,j) = 1.0 + lambda
       END DO

       CALL matinv ( array )

       DO j = 1, nterms
          DO k = 1, nterms
             c(j) = c(j) + beta(k) * array(j,k) / SQRT ( alpha(j,j) * alpha(k,k) )
          END DO
       END DO
       c_new = c

       DO i = 1,npts
          yfit(i) = func(x(i),c_new)
       END DO
       chisq     = SUM ( (y-yfit)**2 )
       chisq_red = chisq / REAL(nfree)
       change    = ( chisq_old - chisq_red ) / chisq_old

       IF ( change > 0.0 ) THEN ! Better fit

          c = c_new
          FORALL ( t = 1:nterms ) sigma(t) = SQRT ( array(t,t) / alpha(t,t) )

          IF ( verbose ) THEN
             WRITE ( unit=output_unit, fmt='(i5)', advance='no' ) iter
             DO t = 1, nterms
                WRITE ( unit=output_unit, fmt='(2f12.6)', advance='no' ) c(t), sigma(t)
             END DO
             WRITE ( unit=output_unit, fmt='(3es12.2)' ) chisq_red, change, lambda
          END IF

          IF ( change < tol ) THEN ! Successful exit
             fail = .FALSE.
             EXIT
          END IF

          lambda = lambda/10.0 ! Improving: try again

       ELSE ! Worse fit

          IF ( lambda > 0.9 ) THEN ! Unsuccessful exit
             IF ( verbose ) WRITE ( unit=output_unit, fmt='(a)') '*** NOT CONVERGED ***'
             fail = .TRUE.
             EXIT
          ENDIF

          lambda = lambda*10.0 ! Worsening: try again

       END IF

    END DO ! End loop until change below tol or failure to converge

  END SUBROUTINE fit

  SUBROUTINE matinv ( arr )
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(inout) :: arr

    ! Invert matrix by Gauss method

    REAL,    DIMENSION(SIZE(arr,1),SIZE(arr,1)) :: a
    REAL,    DIMENSION(SIZE(arr,1))             :: temp
    INTEGER, DIMENSION(SIZE(arr,1))             :: pivot

    INTEGER :: i, j, k, m, n
    REAL    :: c, d

    n = SIZE(arr,1)
    IF ( n /= SIZE(arr,2) ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i5)') 'Array not square', n, SIZE(arr,2)
       STOP 'Error in matinv'
    END IF

    a = arr ! Working copy

    pivot = [ (i, i = 1, n) ]

    DO k = 1, n
       m = k - 1 + MAXLOC(ABS(a(k:n,k)),dim=1)

       IF (m /= k) THEN
          pivot([m,k]) = pivot([k,m]) ! Swap ( pivot(m), pivot(k) )
          a([m,k],:)   = a([k,m],:)   ! Swap ( a(m,:),   a(k,:) )
       END IF
       d = 1.0 / a(k,k)

       temp = a(:,k)
       DO j = 1, n
          c      = a(k,j)*d
          a(:,j) = a(:,j)-temp*c
          a(k,j) = c
       END DO
       a(:,k) = temp*(-d)
       a(k,k) = d
    END DO

    arr(:,pivot) = a

  END SUBROUTINE matinv

END MODULE grint_module
