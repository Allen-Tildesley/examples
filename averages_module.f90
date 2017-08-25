! averages_module.f90
! Calculation of run averages with output to output_unit
MODULE averages_module

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

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit, iostat_end, iostat_eor

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: run_begin, run_end, blk_begin, blk_end, blk_add, time_stamp

  ! Public data
  INTEGER, PARAMETER, PUBLIC :: avg = 0, msd = 1, cke = 2 ! Options for averaging methods

  ! Private data
  INTEGER,          PARAMETER :: col_width = 15              ! Must be large enough to allow sensible format
  INTEGER,          PARAMETER :: nam_width = 2*col_width+1   ! At most two column widths plus spacer
  CHARACTER(len=*), PARAMETER :: colf_fmt  = 'f15.6'         ! Format for floats; we assume that 6 dp will be sufficient
  CHARACTER(len=*), PARAMETER :: cole_fmt  = 'es15.4'        ! Alternative format for floats, suitable for small numbers
  CHARACTER(len=*), PARAMETER :: head_fmt  = '(*(1x,a15))'   ! Format for heading strings
  CHARACTER(len=*), PARAMETER :: col1a_fmt = '(a15)'         ! Format for column 1 strings
  CHARACTER(len=*), PARAMETER :: col1i_fmt = '(i15)'         ! Format for column 1 integers
  CHARACTER(len=*), PARAMETER :: sngl_fmt  = '(a,t40,f15.6)' ! Format for single line output

  INTEGER,                                             SAVE :: n_avg, line_width
  CHARACTER(len=col_width), DIMENSION(:), ALLOCATABLE, SAVE :: headings, subheads
  REAL,                     DIMENSION(:), ALLOCATABLE, SAVE :: blk_avg, blk_msd, run_avg, run_err, add
  INTEGER,                  DIMENSION(:), ALLOCATABLE, SAVE :: method
  REAL,                                                SAVE :: run_nrm, blk_nrm
  CHARACTER(len=:),                       ALLOCATABLE, SAVE :: line_fmt ! Format for single line of output

  ! Public derived type for variables to average
  TYPE, PUBLIC :: variable_type
     CHARACTER(len=nam_width) :: nam                ! Name to be used in headings
     REAL                     :: val                ! Instantaneous value to be averaged
     INTEGER                  :: method = avg       ! Selects method: avg (default), msd, or cke
     REAL                     :: add = 0.0          ! Constant to use in msd method (default zero)
     LOGICAL                  :: e_format = .FALSE. ! Flag to indicate scientific format (small numbers)
     LOGICAL                  :: instant = .TRUE.   ! Flag to indicate whether to print instantaneous value
  END TYPE variable_type

CONTAINS

  SUBROUTINE time_stamp
    IMPLICIT NONE

    CHARACTER(len=8)  :: date
    CHARACTER(len=10) :: time
    REAL              :: cpu

    CALL DATE_AND_TIME ( date, time )
    CALL CPU_TIME ( cpu )
    WRITE ( unit=output_unit, fmt='(a,t45,a4,a1,a2,a1,a2)' ) 'Date: ', date(1:4), '/', date(5:6), '/', date(7:8)
    WRITE ( unit=output_unit, fmt='(a,t47,a2,a1,a2,a1,a2)' ) 'Time: ', time(1:2), ':', time(3:4), ':', time(5:6)
    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'          ) 'CPU time: ', cpu

  END SUBROUTINE time_stamp

  SUBROUTINE run_begin ( variables )
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(:), INTENT(in) :: variables   ! Variables to be averaged

    ! Set up averaging variables based on supplied arrays of names & write headings

    INTEGER :: i, space

    n_avg = SIZE ( variables ) ! Set the number of variables to average

    ALLOCATE ( headings(n_avg), subheads(n_avg) )
    ALLOCATE ( blk_avg(n_avg), blk_msd(n_avg) )
    ALLOCATE ( run_avg(n_avg), run_err(n_avg) )
    ALLOCATE ( method(n_avg), add(n_avg) )

    ! First column plus a column for each variable; allow one space between columns
    line_width = col_width + n_avg * ( col_width + 1 )

    line_fmt = '(' ! Fortran 2003 will automatically allocate and reallocate as we build this up

    ! Store variable names locally in tidied-up format
    ! Attempt to split name string at first space
    ! Build up format string for line of averages
    DO i = 1, n_avg
       space = SCAN ( TRIM(variables(i)%nam), ' ' )
       IF ( space > 0 ) THEN
          headings(i) = variables(i)%nam(1:space-1)
          subheads(i) = variables(i)%nam(space+1:)
       ELSE
          headings(i) = variables(i)%nam(1:col_width)
          subheads(i) = ''
       END IF
       headings(i) = ADJUSTR ( headings(i) )
       subheads(i) = ADJUSTR ( subheads(i) )

       IF ( variables(i)%e_format ) THEN
          line_fmt = line_fmt // '1x,' // cole_fmt // ','
       ELSE
          line_fmt = line_fmt // '1x,' // colf_fmt // ','
       END IF

    END DO

    ! Replace last comma in format string by closing parenthesis
    i = LEN(line_fmt)
    line_fmt(i:i) = ')'

    ! Store method options and add-constants locally
    method = variables%method
    add    = variables%add

    ! Zero averages and error accumulators at start of run
    run_nrm = 0.0
    run_avg = 0.0
    run_err = 0.0

    ! Write initial instantaneous values
    IF ( ANY(variables(:)%instant) ) WRITE ( unit=output_unit, fmt='(a)' ) 'Initial values'
    DO i = 1, n_avg
       IF ( variables(i)%instant ) THEN
          WRITE ( unit=output_unit, fmt=sngl_fmt ) variables(i)%nam, variables(i)%val
       END IF
    END DO

    ! Write headings
    WRITE ( unit=output_unit, fmt=* )
    WRITE ( unit=output_unit, fmt='(/,a)' ) 'Run begins'
    CALL time_stamp
    WRITE ( unit=output_unit, fmt=* )
    WRITE ( unit=output_unit, fmt='(a)'                   ) REPEAT ( '=', line_width ) 
    WRITE ( unit=output_unit, fmt=col1a_fmt, advance='no' ) 'Block'
    WRITE ( unit=output_unit, fmt=head_fmt                ) headings
    WRITE ( unit=output_unit, fmt=col1a_fmt, advance='no' ) '     '
    WRITE ( unit=output_unit, fmt=head_fmt                ) subheads
    WRITE ( unit=output_unit, fmt='(a)'                   ) REPEAT ( '-', line_width )

  END SUBROUTINE run_begin

  SUBROUTINE blk_begin
    IMPLICIT NONE

    ! Zero averaging variables at start of each block
    blk_nrm = 0.0
    blk_avg = 0.0
    blk_msd = 0.0

  END SUBROUTINE blk_begin

  SUBROUTINE blk_add ( variables )
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(:), INTENT(in) :: variables ! Instantaneous values of variables

    ! Increment block-average variables

    IF ( SIZE(variables) /= n_avg ) THEN ! Check for inconsistency in calling program
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Mismatched variable arrays', n_avg, SIZE(variables)
       STOP 'Error in blk_add'
    END IF ! End check for inconsistency in calling program

    blk_avg = blk_avg + variables%val    ! Increment block averages of variables
    blk_msd = blk_msd + variables%val**2 ! Increment block averages of squared variables
    blk_nrm = blk_nrm + 1.0              ! Increment block normalizer

  END SUBROUTINE blk_add

  SUBROUTINE blk_end ( blk )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: blk ! Block number

    ! Write out averages at end of every block

    IF ( blk_nrm < 0.5 ) THEN ! Check for no accumulation; should never happen
       WRITE ( unit=error_unit, fmt='(a,f15.6)' ) 'Block accumulation error', blk_nrm
       STOP 'Error in blk_end'
    END IF ! End check for no accumulation

    blk_avg = blk_avg / blk_nrm ! Normalize block averages
    blk_msd = blk_msd / blk_nrm ! Normalize block averages of squared variables

    ! Replace blk_avg by mean-squared deviations plus optional constant where required
    WHERE ( method == msd .OR. method == cke ) blk_avg = add + blk_msd - blk_avg**2
    IF ( ANY ( method == cke ) ) CALL cke_calc ! Call special routine for Cv from KE fluctuations

    run_avg = run_avg + blk_avg    ! Increment run averages
    run_err = run_err + blk_avg**2 ! Increment error accumulators
    run_nrm = run_nrm + 1.0        ! Increment run normalizer

    ! Write out block averages
    WRITE ( unit=output_unit, fmt=col1i_fmt, advance='no' ) blk
    WRITE ( unit=output_unit, fmt=line_fmt                ) blk_avg

  END SUBROUTINE blk_end

  SUBROUTINE run_end ( variables )
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(:), INTENT(in) :: variables ! Instantaneous variables

    ! Write out averages and error estimates at end of run
    ! NB, these are the crudest possible error estimates, based on the wholly unjustified
    ! assumption that the blocks are statistically independent
    ! For a discussion of errors, see Chapter 8 and the error_calc.f90 example

    INTEGER :: i

    IF ( run_nrm < 0.5 ) THEN ! Check for no accumulation; should never happen
       WRITE ( unit=error_unit, fmt='(a,f15.6)' ) 'Run accumulation error', run_nrm
       STOP 'Error in run_end'
    END IF ! End check for no accumulation

    run_avg = run_avg / run_nrm    ! Normalize run averages
    run_err = run_err / run_nrm    ! Normalize error estimates
    run_err = run_err - run_avg**2 ! Compute fluctuations of block averages

    WHERE ( run_err > 0.0 ) ! Guard against roundoff
       run_err = SQRT ( run_err / run_nrm ) ! Normalize and get estimated errors
    END WHERE ! End guard against roundoff

    WRITE ( unit=output_unit, fmt='(a)'                   ) REPEAT('-', line_width )
    WRITE ( unit=output_unit, fmt=col1a_fmt, advance='no' ) 'Run averages'
    WRITE ( unit=output_unit, fmt=line_fmt                ) run_avg
    WRITE ( unit=output_unit, fmt=col1a_fmt, advance='no' ) 'Run errors'
    WRITE ( unit=output_unit, fmt=line_fmt                ) run_err
    WRITE ( unit=output_unit, fmt='(a)'                   ) REPEAT('=', line_width )
    WRITE ( unit=output_unit, fmt=* )
    WRITE ( unit=output_unit, fmt='(a)' ) 'Run ends'
    CALL time_stamp
    WRITE ( unit=output_unit, fmt=* )

    ! Write final instantaneous values
    IF ( ANY(variables(:)%instant) ) WRITE ( unit=output_unit, fmt='(/,a)' ) 'Final values'
    DO i = 1, SIZE ( variables )
       IF ( variables(i)%instant ) THEN
          WRITE ( unit=output_unit, fmt=sngl_fmt ) variables(i)%nam, variables(i)%val
       END IF
    END DO

    DEALLOCATE ( headings, blk_avg, blk_msd, run_avg, run_err, method, add )

  END SUBROUTINE run_end

  SUBROUTINE cke_calc
    IMPLICIT NONE

    INTEGER :: i
    LOGICAL :: found
    REAL    :: temperature

    ! Locate variable corresponding to kinetic temperature

    found = .FALSE.
    DO i = 1, n_avg
       IF ( INDEX ( headings(i)//subheads(i),'T' ) /= 0 .AND. INDEX ( headings(i)//subheads(i), 'kin' ) /= 0 ) THEN
          temperature = blk_avg(i)
          found       = .TRUE.
          EXIT
       END IF
    END DO

    IF ( .NOT. found ) THEN
       WRITE ( unit=error_unit, fmt='(a)' ) 'Could not find T kin variable'
       STOP 'Error in cke_calc'
    END IF

    ! Apply special fluctuation formula for microcanonical ensemble heat capacity
    ! blk_avg(i) should contain mean-squared total KE, divided by N
    
    DO i = 1, n_avg
       IF ( method(i) == cke ) THEN
          blk_avg(i) = 9.0 / ( 6.0 - 4.0 * blk_avg(i) / temperature**2 )
       END IF
    END DO

  END SUBROUTINE cke_calc
  
END MODULE averages_module
