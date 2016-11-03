! averages_module.f90
! Calculation of run averages with output to output_unit
MODULE averages_module

  ! We use the standard error_unit for error messages
  ! but allow the output_unit to be passed in as an argument, when needed
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit, iostat_end, iostat_eor

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: run_begin, run_end, blk_begin, blk_end, blk_add, time_stamp

  ! Private data
  INTEGER,                                     SAVE :: nvariables, col_width, line_width
  CHARACTER(len=:), DIMENSION(:), ALLOCATABLE, SAVE :: variable_names
  REAL,             DIMENSION(:), ALLOCATABLE, SAVE :: blk_averages, run_averages, run_errors
  REAL,                                        SAVE :: run_norm, blk_norm
  CHARACTER(len=:),               ALLOCATABLE, SAVE :: col_fmt

CONTAINS

  SUBROUTINE time_stamp ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit

    CHARACTER(len=8)    :: date
    CHARACTER(len=10)   :: time
    REAL                :: cpu

    CALL DATE_AND_TIME ( date, time )
    CALL CPU_TIME ( cpu )
    WRITE ( unit=output_unit, fmt='(a,t45,a4,a1,a2,a1,a2)' ) 'Date: ', date(1:4), '/', date(5:6), '/', date(7:8)
    WRITE ( unit=output_unit, fmt='(a,t47,a2,a1,a2,a1,a2)' ) 'Time: ', time(1:2), ':', time(3:4), ':', time(5:6)
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)'          ) 'CPU time: ', cpu

  END SUBROUTINE time_stamp

  SUBROUTINE run_begin ( names )
    IMPLICIT NONE
    CHARACTER(len=*), DIMENSION(:), INTENT(in) :: names

    ! Set up averaging variables based on supplied array of names

    INTEGER            :: i
    CHARACTER(len=2)   :: ch_col_width       ! Column-width part of format string for columns
    INTEGER, PARAMETER :: max_col_width = 30 ! Must be < 100 to fit in ch_col_width
    INTEGER, PARAMETER :: min_col_width = 15 ! Must be large enough to allow sensible f##.5 format
    INTEGER, PARAMETER :: col_one_width = 15 ! First column width for block number etc

    nvariables = SIZE ( names ) ! This is how we set the number of variables to average

    col_width = MAX ( min_col_width, LEN(names) ) ! We don't allow columns to be too narrow

    IF ( col_width > max_col_width ) THEN ! We don't allow columns to be too large
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'col_width too large', col_width
       STOP 'Error in run_begin'
    END IF

    ! Set up format string for columns; in Fortran 2003 col_fmt is allocated automatically
    WRITE ( ch_col_width, fmt='(i2)' ) col_width  ! Character string for ## part of f##.5
    col_fmt = '(*(1x,f' // ch_col_width // '.5))' ! We assume that 5 dp will be sufficient

    ! First column plus a column for each variable; allow one space between columns 
    line_width = col_one_width + nvariables * ( col_width + 1 )

    ALLOCATE ( CHARACTER(col_width) :: variable_names(nvariables) )
    ALLOCATE ( blk_averages(nvariables) )
    ALLOCATE ( run_averages(nvariables) )
    ALLOCATE ( run_errors  (nvariables) )

    ! Store variable names in tidied-up format
    variable_names = names
    DO i = 1, nvariables
       variable_names(i) = ADJUSTR ( variable_names(i) )
    END DO

    ! Zero averages and error accumulators at start of run
    run_norm     = 0.0
    run_averages = 0.0
    run_errors   = 0.0

  END SUBROUTINE run_begin

  SUBROUTINE blk_begin
    IMPLICIT NONE

    ! Zero averaging variables at start of each block
    blk_norm     = 0.0
    blk_averages = 0.0

  END SUBROUTINE blk_begin

  SUBROUTINE blk_add ( variables )
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(in) :: variables

    ! Increment block-average variables

    IF ( SIZE(variables) /= nvariables ) THEN ! Check for inconsistency in calling program
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Mismatched variable arrays', nvariables, SIZE(variables)
       STOP 'Error in blk_add'
    END IF ! End check for inconsistency in calling program

    blk_averages = blk_averages + variables ! Increment block averages
    blk_norm     = blk_norm + 1.0           ! Increment block normalizer

  END SUBROUTINE blk_add

  SUBROUTINE blk_end ( blk, output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: blk, output_unit

    ! Write out averages at end of every block

    LOGICAL, SAVE :: first_call = .TRUE.

    IF ( blk_norm < 0.5 ) THEN ! Check for no accumulation; should never happen
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'Block accumulation error', blk_norm
       STOP 'Error in blk_end'
    END IF ! End check for no accumulation

    blk_averages = blk_averages / blk_norm ! Normalize block averages

    run_averages = run_averages + blk_averages    ! Increment run averages
    run_errors   = run_errors   + blk_averages**2 ! Increment error accumulators
    run_norm     = run_norm     + 1.0             ! Increment run normalizer

    IF ( first_call ) THEN ! Write headings on first call
       WRITE ( unit=output_unit, fmt='(a)'                 ) REPEAT ( '=', line_width ) 
       WRITE ( unit=output_unit, fmt='(a15)', advance='no' ) 'Block'
       WRITE ( unit=output_unit, fmt='(*(1x,a))'           ) variable_names
       WRITE ( unit=output_unit, fmt='(a)'                 ) REPEAT ( '=', line_width )
       first_call = .FALSE.
    END IF ! End write headings on first call

    ! Write out block averages
    WRITE ( unit=output_unit, fmt='(i15)', advance='no' ) blk
    WRITE ( unit=output_unit, fmt=col_fmt               ) blk_averages

  END SUBROUTINE blk_end

  SUBROUTINE run_end ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit

    ! Write out averages and error estimates at end of run
    ! NB, these are the crudest possible error estimates, based on the wholly unjustified
    ! assumption that the blocks are statistically independent
    ! For a discussion of errors, see Chapter 8 and the error_calc.f90 example

    IF ( run_norm < 0.5 ) THEN ! Check for no accumulation; should never happen
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'Run accumulation error', run_norm
       STOP 'Error in run_end'
    END IF ! End check for no accumulation

    run_averages = run_averages / run_norm      ! Normalize run averages
    run_errors   = run_errors / run_norm        ! Normalize error estimates
    run_errors   = run_errors - run_averages**2 ! Compute fluctuations of block averages

    WHERE ( run_errors > 0.0 ) ! Guard against roundoff
       run_errors = SQRT ( run_errors / run_norm ) ! Normalize and get estimated errors
    END WHERE ! End guard against roundoff

    WRITE ( unit=output_unit, fmt='(a)'                  ) REPEAT('-', line_width )
    WRITE ( unit=output_unit, fmt='(a15)', advance='no'  ) 'Run averages'
    WRITE ( unit=output_unit, fmt=col_fmt                ) run_averages
    WRITE ( unit=output_unit, fmt='(a15)', advance='no'  ) 'Run errors'
    WRITE ( unit=output_unit, fmt=col_fmt                ) run_errors
    WRITE ( unit=output_unit, fmt='(a)'                  ) REPEAT('=', line_width )

    DEALLOCATE ( variable_names, blk_averages, run_averages, run_errors, col_fmt )

  END SUBROUTINE run_end

END MODULE averages_module
