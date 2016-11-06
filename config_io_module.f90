! config_io_module.f90
! routines for atomic/molecular configuration I/O
MODULE config_io_module

  ! We use the standard error_unit for error messages
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit, iostat_end, iostat_eor

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: read_cnf_atoms, write_cnf_atoms, read_cnf_mols, write_cnf_mols

CONTAINS

  SUBROUTINE read_cnf_atoms ( filename, n, box, r, v ) ! Read in atomic configuration
    IMPLICIT NONE
    CHARACTER(len=*),               INTENT(in)    :: filename ! Supplied filename
    INTEGER,                        INTENT(inout) :: n        ! Number of atoms
    REAL,                           INTENT(out)   :: box      ! Simulation box length (assumed cubic)
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r        ! Atomic positions (3,n)
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: v        ! Atomic velocities (3,n)

    ! The first call of this routine is just to get n and box, used to allocate arrays
    ! The second call attempts to read in the atomic positions and optionally velocities

    INTEGER :: cnf_unit, ioerr, i

    ! Open given filename, will terminate on any errors

    OPEN ( newunit=cnf_unit, file=filename, status='old', action='read', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
       STOP 'Error in read_cnf_atoms'
    END IF

    ! Read number of atoms from first record
    
    READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) n
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading n from ', filename, ioerr
       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
       STOP 'Error in read_cnf_atoms'
    END IF

    ! Read box length from second record
    
    READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) box
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading box from ', filename, ioerr
       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
       STOP 'Error in read_cnf_atoms'
    END IF

    ! Expected format is one record per atom containing either r(:,i), v(:,i) or just r(:,i)

    IF ( PRESENT ( r ) ) THEN

       IF ( n /= SIZE ( r, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
          STOP 'Error in read_cnf_atoms'
       END IF

       IF ( PRESENT ( v ) ) THEN

          IF ( n /= SIZE ( v, dim=2 ) ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
             STOP 'Error in read_cnf_atoms'
          END IF

          ! Read positions, velocities
          DO i = 1, n
             READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i), v(:,i)
             IF ( ioerr /= 0 ) THEN
                WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r, v from ', filename, ioerr
                IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                STOP 'Error in read_cnf_atoms'
             END IF
          END DO

       ELSE

          ! Read positions
          DO i = 1, n
             READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i)
             IF ( ioerr /= 0 ) THEN
                WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r from ', filename, ioerr
                IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                STOP 'Error in read_cnf_atoms'
             END IF
          END DO

       END IF

    END IF

    CLOSE ( unit=cnf_unit )

  END SUBROUTINE read_cnf_atoms

  SUBROUTINE write_cnf_atoms ( filename, n, box, r, v ) ! Write out atomic configuration
    IMPLICIT NONE
    CHARACTER(len=*),               INTENT(in) :: filename ! Supplied filename
    INTEGER,                        INTENT(in) :: n        ! Number of atoms
    REAL,                           INTENT(in) :: box      ! Simulation box length (assumed cubic)
    REAL, DIMENSION(:,:),           INTENT(in) :: r        ! Atomic positions (3,n)
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: v        ! Atomic velocities

    INTEGER :: cnf_unit, ioerr, i

    ! Open given filename, replacing it if it already exists, will terminate on any errors
    OPEN ( newunit=cnf_unit, file=filename, status='replace', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
       STOP 'Error in write_cnf_atoms'
    END IF

    ! Write number of atoms to first record, box length to second record
    WRITE ( unit=cnf_unit, fmt='(i15)'  ) n
    WRITE ( unit=cnf_unit, fmt='(f15.8)') box

    IF ( n /= SIZE ( r, dim=2 ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
       STOP 'Error in write_cnf_atoms'
    END IF

    IF ( PRESENT ( v ) ) THEN

       IF ( n /= SIZE ( v, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
          STOP 'Error in write_cnf_atoms'
       END IF

       ! Write positions, velocities
       DO i = 1, n
          WRITE ( unit=cnf_unit, fmt='(*(f15.10))' ) r(:,i), v(:,i)
       END DO

    ELSE

       ! Write positions
       DO i = 1, n
          WRITE ( unit=cnf_unit, fmt='(*(f15.10))' ) r(:,i)
       END DO

    END IF

    CLOSE ( unit=cnf_unit )

  END SUBROUTINE write_cnf_atoms

  SUBROUTINE read_cnf_mols ( filename, n, box, r, e, v, w ) ! Read in molecular configuration
    IMPLICIT NONE
    CHARACTER(len=*),               INTENT(in)    :: filename ! Supplied filename
    INTEGER,                        INTENT(inout) :: n        ! Number of molecules
    REAL,                           INTENT(out)   :: box      ! Simulation box length (assumed cubic)
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r        ! Molecular positions (3,n)
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: e        ! Orientation vectors (3,n) or quaternions (4,n)
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: v        ! Molecular velocities (3,n)
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: w        ! Angular velocities (3,n)

    ! The first call of this routine is just to get n and box, used to allocate arrays
    ! The second call attempts to read in the molecular positions, orientations
    ! and optionally velocities, angular velocities

    INTEGER :: cnf_unit, ioerr, i

    ! Open given filename, will terminate on any errors

    OPEN ( newunit=cnf_unit, file=filename, status='old', action='read', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
       STOP 'Error opening file in read_cnf_mols'
    END IF

    ! Read number of molecules from first record
    
    READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) n
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading n from ', filename, ioerr
       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
       STOP 'Error in read_cnf_mols'
    END IF

    ! Read box length from second record
    
    READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) box
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading box from ', filename, ioerr
       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
       STOP 'Error in read_cnf_mols'
    END IF

    ! Expected format is one line per atom containing either r(:,i), e(:,i), v(:,i), w(:,i)  or just r(:,i), e(:,i)
    ! The first dimension of the e array can be 3 elements (vector) or 4 elements (quaternion)

    IF ( PRESENT ( r ) ) THEN

       IF ( .NOT. PRESENT ( e )    ) THEN
          WRITE ( unit=error_unit, fmt='(a,a,i15)') 'r and e arguments must be present together'
          STOP 'Error in read_cnf_mols'
       END IF
       IF ( n /= SIZE ( r, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
          STOP 'Error in read_cnf_mols'
       END IF
       IF ( n /= SIZE ( e, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of e ', n, SIZE ( e, dim=2 )
          STOP 'Error in read_cnf_mols'
       END IF

       IF ( PRESENT ( v ) ) THEN

          IF ( .NOT. PRESENT ( w )    ) THEN
             WRITE ( unit=error_unit, fmt='(a,a,i15)') 'v and w arguments must be present together'
             STOP 'Error in read_cnf_mols'
          END IF
          IF ( n /= SIZE ( v, dim=2 ) ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
             STOP 'Error in read_cnf_mols'
          END IF
          IF ( n /= SIZE ( w, dim=2 ) ) THEN
             WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of w ', n, SIZE ( w, dim=2 )
             STOP 'Error in read_cnf_mols'
          END IF

          ! Read positions, orientation vectors or quaternions, velocities, angular velocities
          DO i = 1, n
             READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i), e(:,i), v(:,i), w(:,i)
             IF ( ioerr /= 0 ) THEN
                WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r, e, v, w from ', filename, ioerr
                IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                STOP 'Error in read_cnf_mols'
             END IF
          END DO

       ELSE

          ! Read positions, orientation vectors or quaternions
          DO i = 1, n
             READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i), e(:,i)
             IF ( ioerr /= 0 ) THEN
                WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r, e from ', filename, ioerr
                IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                STOP 'Error in read_cnf_mols'
             END IF
          END DO

       END IF

    END IF

    CLOSE ( unit=cnf_unit )

  END SUBROUTINE read_cnf_mols

  SUBROUTINE write_cnf_mols ( filename, n, box, r, e, v, w ) ! Write out molecular configuration
    IMPLICIT NONE
    CHARACTER(len=*),               INTENT(in) :: filename ! Supplied filename
    INTEGER,                        INTENT(in) :: n        ! Number of molecules
    REAL,                           INTENT(in) :: box      ! Simulation box length (assumed cubic)
    REAL, DIMENSION(:,:),           INTENT(in) :: r        ! Molecular positions (3,n)
    REAL, DIMENSION(:,:),           INTENT(in) :: e        ! Orientation vectors (3,n) or quaternions (0:3,n)
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: v        ! Molecular velocities (3,n)
    REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: w        ! Angular velocities (3,n)

    INTEGER :: cnf_unit, ioerr, i

    ! Open given filename, replacing it if it already exists, will terminate on any errors
    OPEN ( newunit=cnf_unit, file=filename, status='replace', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
       STOP 'Error in write_cnf_mols'
    END IF

    ! Write number of molecules to first record, box length to second record
    WRITE(cnf_unit,'(i15)'  ) n
    WRITE(cnf_unit,'(f15.8)') box

    IF ( n /= SIZE ( r, dim=2 ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
       STOP 'Error in write_cnf_mols'
    END IF
    IF ( n /= SIZE ( e, dim=2 ) ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of e ', n, SIZE ( e, dim=2 )
       STOP 'Error in write_cnf_mols'
    END IF

    IF ( PRESENT ( v ) ) THEN
       IF ( .NOT. PRESENT ( w )    ) THEN
          WRITE ( unit=error_unit, fmt='(a,a,i15)') 'v and w arguments must be present together'
          STOP 'Error in write_cnf_mols'
       END IF
       IF ( n /= SIZE ( v, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
          STOP 'Error in write_cnf_mols'
       END IF
       IF ( n /= SIZE ( w, dim=2 ) ) THEN
          WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of w ', n, SIZE ( w, dim=2 )
          STOP 'Error in write_cnf_mols'
       END IF

       ! Write positions, orientation vectors or quaternions, velocities, angular velocities
       DO i = 1, n
          WRITE ( unit=cnf_unit, fmt='(*(f15.10))') r(:,i), e(:,i), v(:,i), w(:,i)
       END DO

    ELSE

       ! Write positions, orientation vectors or quaternions
       DO i = 1, n
          WRITE ( unit=cnf_unit, fmt='(*(f15.10))') r(:,i), e(:,i)
       END DO

    END IF

    CLOSE ( unit=cnf_unit )

  END SUBROUTINE write_cnf_mols

END MODULE config_io_module
