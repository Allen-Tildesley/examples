MODULE io_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_cnf_atoms, write_cnf_atoms

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

END MODULE io_module
