MODULE utility_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_cnf_atoms, write_cnf_atoms
  PUBLIC :: run_begin, run_end, blk_begin, blk_end, stp_end

  INTEGER,                                      SAVE :: nvariables
  CHARACTER(len=10), DIMENSION(:), ALLOCATABLE, SAVE :: variable_names
  REAL,              DIMENSION(:), ALLOCATABLE, SAVE :: blk_averages, run_averages, errors
  REAL,                                         SAVE :: run_norm, blk_norm

CONTAINS

  SUBROUTINE read_cnf_atoms ( filename, n, box, r, v ) ! Read in atomic configuration
    CHARACTER(len=*),               INTENT(in)    :: filename
    INTEGER,                        INTENT(inout) :: n
    REAL,                           INTENT(out)   :: box
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r, v

    INTEGER :: cnf_unit, ioerr, i

    OPEN(newunit=cnf_unit,file=filename,status='old',action='read',iostat=ioerr)
    IF ( ioerr /= 0 ) STOP 'Error opening file in read_cnf_atoms'
    READ(cnf_unit,*) n
    READ(cnf_unit,*) box

    IF ( PRESENT (r ) ) THEN
       IF ( PRESENT ( v ) ) THEN
          DO i = 1, n
             READ(cnf_unit,*) r(:,i), v(:,i) ! positions and velocities
          END DO
       ELSE
          DO i = 1, n
             READ(cnf_unit,*) r(:,i) ! positions
          END DO
       END IF
    END IF
    CLOSE(unit=cnf_unit)

  END SUBROUTINE read_cnf_atoms

  SUBROUTINE write_cnf_atoms ( filename, n, box, r, v )
    CHARACTER(len=*),               INTENT(in) :: filename
    INTEGER,                        INTENT(in) :: n
    REAL,                           INTENT(in) :: box
    REAL, DIMENSION(:,:),           INTENT(in) :: r
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: v

    INTEGER :: cnf_unit, ioerr, i

    OPEN(newunit=cnf_unit,file=filename,status='replace',iostat=ioerr)
    IF ( ioerr /= 0 ) STOP 'Error opening file in write_cnf_atoms'
    WRITE(cnf_unit,'(i15)'  ) n
    WRITE(cnf_unit,'(f15.8)') box
    IF ( PRESENT ( v ) ) THEN
       DO i = 1, n
          WRITE(cnf_unit,'(6f15.10)') r(:,i), v(:,i) ! positions and velocities
       END DO
    ELSE
       DO i = 1, n
          WRITE(cnf_unit,'(3f15.10)') r(:,i) ! positions
       END DO
    END IF
    CLOSE(unit=cnf_unit)

  END SUBROUTINE write_cnf_atoms

  SUBROUTINE run_begin ( names )
    CHARACTER(len=10), DIMENSION(:), INTENT(in) :: names

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

  SUBROUTINE blk_begin
    blk_norm     = 0.0
    blk_averages = 0.0
  END SUBROUTINE blk_begin

  SUBROUTINE stp_end ( variables )
    REAL, DIMENSION(:), INTENT(in) :: variables

    IF ( SIZE(variables) /= nvariables ) STOP 'mismatched variable arrays in stp_end'
    blk_averages = blk_averages + variables ! Increment block averages
    blk_norm     = blk_norm + 1.0           ! Increment block normalizer
  END SUBROUTINE stp_end

  SUBROUTINE blk_end ( blk )
    INTEGER, INTENT(in) :: blk

    LOGICAL, SAVE :: first_call = .TRUE.

    blk_averages = blk_averages / blk_norm     ! Normalize block averages
    run_averages = run_averages + blk_averages ! Increment run averages
    errors       = errors + blk_averages**2    ! Increment error accumulators
    run_norm     = run_norm + 1.0              ! Increment run normalizer

    IF ( first_call ) THEN  ! Write headings
       WRITE(*,'(*(a15))') REPEAT ( '=', 15*(nvariables+1) ) 
       WRITE(*,'(*(5x,a10))') 'Block     ', variable_names
       WRITE(*,'(*(a15))') REPEAT ( '=', 15*(nvariables+1) )
       first_call = .FALSE.
    END IF

    ! Write out block averages
    WRITE(*,'(5x,i10,*(5x,f10.4))') blk, blk_averages

  END SUBROUTINE blk_end

  SUBROUTINE run_end

    run_averages = run_averages / run_norm  ! Normalize run averages
    errors       = errors / run_norm        ! Normalize error estimates
    errors       = errors - run_averages**2 ! Compute fluctuations
    WHERE ( errors > 0.0 )
       errors = SQRT ( errors / run_norm ) ! Normalize and get estimated errors
    END WHERE

    WRITE(*,'(*(a15))') REPEAT('-',15*(nvariables+1))
    WRITE(*,'(a15,*(5x,f10.4))') 'Run averages', run_averages
    WRITE(*,'(a15,*(5x,f10.4))') 'Run errors', errors
    WRITE(*,'(*(a15))') REPEAT('=',15*(nvariables+1))

  END SUBROUTINE run_end

END MODULE utility_module
