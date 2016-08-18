! averages_module.f90
! Calculation of run averages with output to output_unit
MODULE averages_module

  ! We use the standard error_unit for error messages
  ! but allow the output_unit to be passed in as an argument, when needed
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit, iostat_end, iostat_eor

  IMPLICIT NONE
  PRIVATE

  ! Calculation of averages with output to output_unit
  PUBLIC :: run_begin, run_end, blk_begin, blk_end, blk_add, time_stamp

  ! These variables (run averages etc) are only accessed within this module
  INTEGER,                                      SAVE :: nvariables
  CHARACTER(len=15), DIMENSION(:), ALLOCATABLE, SAVE :: variable_names
  REAL,              DIMENSION(:), ALLOCATABLE, SAVE :: blk_averages, run_averages, errors
  REAL,                                         SAVE :: run_norm, blk_norm

CONTAINS
  
  ! Routines associated with quantities to be averaged and output

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

  SUBROUTINE run_begin ( names ) ! Set up averaging variables based on supplied names
    CHARACTER(len=15), DIMENSION(:), INTENT(in) :: names

    nvariables = SIZE ( names )
    ALLOCATE ( variable_names(nvariables) )
    ALLOCATE ( blk_averages(nvariables) )
    ALLOCATE ( run_averages(nvariables) )
    ALLOCATE ( errors(nvariables) )

    variable_names = names
    run_norm       = 0.0
    run_averages   = 0.0
    errors         = 0.0

  END SUBROUTINE run_begin

  SUBROUTINE blk_begin ! Zero averaging variables at start of each block
    blk_norm     = 0.0
    blk_averages = 0.0
  END SUBROUTINE blk_begin

  SUBROUTINE blk_add ( variables ) ! Increment block-average variables
    REAL, DIMENSION(:), INTENT(in) :: variables

    IF ( SIZE(variables) /= nvariables ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Mismatched variable arrays', nvariables, SIZE(variables)
       STOP 'Error in blk_add'
    END IF

    blk_averages = blk_averages + variables ! Increment block averages
    blk_norm     = blk_norm + 1.0           ! Increment block normalizer
  END SUBROUTINE blk_add

  SUBROUTINE blk_end ( blk, output_unit ) ! Write out block averages
    INTEGER, INTENT(in) :: blk, output_unit

    LOGICAL, SAVE :: first_call = .TRUE.

    blk_averages = blk_averages / blk_norm     ! Normalize block averages
    run_averages = run_averages + blk_averages ! Increment run averages
    errors       = errors + blk_averages**2    ! Increment error accumulators
    run_norm     = run_norm + 1.0              ! Increment run normalizer

    IF ( first_call ) THEN  ! Write headings
       WRITE ( unit=output_unit, fmt='(a)'        ) REPEAT ( '=', 16*(nvariables+1) ) 
       WRITE ( unit=output_unit, fmt='(*(1x,a15))') 'Block', ADJUSTR ( variable_names )
       WRITE ( unit=output_unit, fmt='(a)'        ) REPEAT ( '=', 16*(nvariables+1) )
       first_call = .FALSE.
    END IF

    ! Write out block averages
    WRITE ( unit=output_unit, fmt='(1x,i15,*(1x,f15.5))') blk, blk_averages

  END SUBROUTINE blk_end

  SUBROUTINE run_end ( output_unit ) ! Write out run averages
    INTEGER, INTENT(in) :: output_unit

    run_averages = run_averages / run_norm  ! Normalize run averages
    errors       = errors / run_norm        ! Normalize error estimates
    errors       = errors - run_averages**2 ! Compute fluctuations
    WHERE ( errors > 0.0 )
       errors = SQRT ( errors / run_norm ) ! Normalize and get estimated errors
    END WHERE

    WRITE ( unit=output_unit, fmt='(a)'                  ) REPEAT('-',16*(nvariables+1))
    WRITE ( unit=output_unit, fmt='(1x,a15,*(1x,f15.5))' ) 'Run averages', run_averages
    WRITE ( unit=output_unit, fmt='(1x,a15,*(1x,f15.5))' ) 'Run errors', errors
    WRITE ( unit=output_unit, fmt='(a)'                  ) REPEAT('=',16*(nvariables+1))

    DEALLOCATE ( variable_names, blk_averages, run_averages, errors )

  END SUBROUTINE run_end

END MODULE averages_module
