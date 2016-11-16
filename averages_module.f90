! averages_module.f90
! Calculation of run averages with output to output_unit
MODULE averages_module

  ! We use the standard error_unit for error messages
  ! but allow the output_unit to be passed in as an argument, when needed
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit, iostat_end, iostat_eor

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: run_begin, run_end, blk_begin, blk_end, blk_add, time_stamp, write_variables

  ! Private data
  INTEGER,           PARAMETER :: col_width = 15              ! Must be large enough to allow sensible format
  CHARACTER(len=13), PARAMETER :: colf_fmt  = '(*(1x,f15.5))' ! Format for floats; we assume that 5 dp will be sufficient
  CHARACTER(len=11), PARAMETER :: cola_fmt  = '(*(1x,a15))'   ! Format for strings
  CHARACTER(len=5),  PARAMETER :: col1a_fmt = '(a15)'         ! Format for column 1 strings
  CHARACTER(len=5),  parameter :: col1i_fmt = '(i15)'         ! Format for column 1 integers

  INTEGER,                                             SAVE :: n_avg, line_width
  CHARACTER(len=col_width), DIMENSION(:), ALLOCATABLE, SAVE :: headings
  REAL,                     DIMENSION(:), ALLOCATABLE, SAVE :: blk_avg, blk_msd, run_avg, run_err
  LOGICAL,                  DIMENSION(:), ALLOCATABLE, SAVE :: msd
  REAL,                                                SAVE :: run_nrm, blk_nrm

  ! Public derived type for variables to average
  TYPE, PUBLIC :: variable_type
     CHARACTER(len=col_width) :: nam ! Name to be used in headings
     REAL                     :: val ! Instantaneous value to be averaged
     LOGICAL                  :: msd = .false. ! Flag indicating mean square difference required
  END TYPE variable_type

CONTAINS

  SUBROUTINE time_stamp ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit

    CHARACTER(len=8)  :: date
    CHARACTER(len=10) :: time
    REAL              :: cpu

    CALL DATE_AND_TIME ( date, time )
    CALL CPU_TIME ( cpu )
    WRITE ( unit=output_unit, fmt='(a,t45,a4,a1,a2,a1,a2)' ) 'Date: ', date(1:4), '/', date(5:6), '/', date(7:8)
    WRITE ( unit=output_unit, fmt='(a,t47,a2,a1,a2,a1,a2)' ) 'Time: ', time(1:2), ':', time(3:4), ':', time(5:6)
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)'          ) 'CPU time: ', cpu

  END SUBROUTINE time_stamp

  SUBROUTINE run_begin ( output_unit, variables )
    IMPLICIT NONE
    INTEGER,                           INTENT(in) :: output_unit
    TYPE(variable_type), DIMENSION(:), INTENT(in) :: variables   ! Variables to be averaged

    ! Set up averaging variables based on supplied arrays of names & write headings

    INTEGER :: i

    n_avg = SIZE ( variables ) ! Set the number of variables to average

    ALLOCATE ( headings(n_avg) )
    ALLOCATE ( blk_avg(n_avg), blk_msd(n_avg) )
    ALLOCATE ( run_avg(n_avg), run_err(n_avg) )
    ALLOCATE ( msd(n_avg) )

    ! First column plus a column for each variable; allow one space between columns 
    line_width = col_width + n_avg * ( col_width + 1 )

    ! Store variable names locally in tidied-up format
    headings = variables%nam
    DO i = 1, n_avg
       headings(i) = ADJUSTR ( headings(i) )
    END DO

    ! Store msd flags locally
    msd = variables%msd

    ! Zero averages and error accumulators at start of run
    run_nrm = 0.0
    run_avg = 0.0
    run_err = 0.0

    ! Write headings
    WRITE ( unit=output_unit, fmt='(a)'                   ) REPEAT ( '=', line_width ) 
    WRITE ( unit=output_unit, fmt=col1a_fmt, advance='no' ) 'Block'
    WRITE ( unit=output_unit, fmt=cola_fmt                ) headings
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

  SUBROUTINE blk_end ( blk, output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: blk, output_unit

    ! Write out averages at end of every block

    IF ( blk_nrm < 0.5 ) THEN ! Check for no accumulation; should never happen
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'Block accumulation error', blk_nrm
       STOP 'Error in blk_end'
    END IF ! End check for no accumulation

    blk_avg = blk_avg / blk_nrm ! Normalize block averages
    blk_msd = blk_msd / blk_nrm ! Normalize block averages of squared variables

    ! Replace blk_avg by mean-squared deviations where required
    WHERE ( msd ) blk_avg = blk_msd - blk_avg**2

    run_avg = run_avg + blk_avg    ! Increment run averages
    run_err = run_err + blk_avg**2 ! Increment error accumulators
    run_nrm = run_nrm + 1.0        ! Increment run normalizer

    ! Write out block averages
    WRITE ( unit=output_unit, fmt=col1i_fmt, advance='no' ) blk
    WRITE ( unit=output_unit, fmt=colf_fmt                ) blk_avg

  END SUBROUTINE blk_end

  SUBROUTINE run_end ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit

    ! Write out averages and error estimates at end of run
    ! NB, these are the crudest possible error estimates, based on the wholly unjustified
    ! assumption that the blocks are statistically independent
    ! For a discussion of errors, see Chapter 8 and the error_calc.f90 example

    IF ( run_nrm < 0.5 ) THEN ! Check for no accumulation; should never happen
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'Run accumulation error', run_nrm
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
    WRITE ( unit=output_unit, fmt=colf_fmt                ) run_avg
    WRITE ( unit=output_unit, fmt=col1a_fmt, advance='no' ) 'Run errors'
    WRITE ( unit=output_unit, fmt=colf_fmt                ) run_err
    WRITE ( unit=output_unit, fmt='(a)'                   ) REPEAT('=', line_width )

    DEALLOCATE ( headings, blk_avg, blk_msd, run_avg, run_err, msd )

  END SUBROUTINE run_end

  SUBROUTINE write_variables ( output_unit, variables )
    IMPLICIT NONE
    INTEGER,                           INTENT(in) :: output_unit
    TYPE(variable_type), DIMENSION(:), INTENT(in) :: variables   ! Variables to be written

    ! Writes out instantaneous values

    INTEGER :: i

    DO i = 1, SIZE(variables)
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) variables(i)%nam, variables(i)%val
    END DO

  END SUBROUTINE write_variables

END MODULE averages_module
