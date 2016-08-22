! mc_lj_ll_module.f90
! Energy and move routines for MC, LJ potential, link-lists
MODULE mc_lj_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, lt, ne, gt
  PUBLIC :: allocate_arrays, deallocate_arrays, resize, energy_1, energy, energy_lrc
  PUBLIC :: move, create, destroy

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,:)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: j_list ! list of j-neighbours

  INTEGER, PARAMETER :: lt = -1, ne = 0, gt = 1 ! j-range options

CONTAINS

  SUBROUTINE allocate_arrays ( box, r_cut )
    USE link_list_module, ONLY : initialize_list
    REAL, INTENT(in) :: box   ! simulation box length
    REAL, INTENT(in) :: r_cut ! potential cutoff distance

    REAL :: r_cut_box
    
    ALLOCATE ( r(3,n), j_list(n) )

    r_cut_box = r_cut / box
    IF ( r_cut_box > 0.5 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_cut/box too large ', r_cut_box
       STOP 'Error in allocate_arrays'
    END IF

    CALL initialize_list ( n, r_cut_box )

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    USE link_list_module, ONLY : finalize_list
    DEALLOCATE ( r, j_list )
    CALL finalize_list
  END SUBROUTINE deallocate_arrays

  SUBROUTINE resize ! reallocates r array, twice as large

    ! This is employed by mc_zvt_lj, grand canonical ensemble

    REAL, DIMENSION(:,:), ALLOCATABLE :: tmp
    INTEGER                           :: n_old, n_new

    n_old = SIZE(r,dim=2)
    n_new = 2*n_old
    WRITE( unit=output_unit, fmt='(a,i10,a,i10)' ) 'Reallocating r from old ', n_old, ' to ', n_new

    ALLOCATE ( tmp(3,n_new) ) ! new size for r
    tmp(:,1:n_old) = r(:,:)   ! copy elements across

    CALL MOVE_ALLOC ( tmp, r )

  END SUBROUTINE resize

  SUBROUTINE energy ( box, r_cut, overlap, pot, vir )
    USE link_list_module, ONLY : make_list

    REAL,    INTENT(in)  :: box        ! simulation box length
    REAL,    INTENT(in)  :: r_cut      ! potential cutoff
    LOGICAL, INTENT(out) :: overlap    ! shows if an overlap was detected
    REAL,    INTENT(out) :: pot, vir   ! potential and virial 

    ! Calculates potential and virial for whole system
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the values of pot and vir should not be used
    ! Actual calculation is performed by subroutine energy_1

    REAL               :: pot_i, vir_i
    INTEGER            :: i
    LOGICAL, SAVE      :: first_call = .TRUE.

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy'
    END IF
    IF ( first_call ) THEN
       r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
       CALL make_list ( n, r )
       first_call = .FALSE.
    END IF

    overlap  = .FALSE.
    pot      = 0.0
    vir      = 0.0

    DO i = 1, n
       CALL energy_1 ( r(:,i), i, gt, box, r_cut, overlap, pot_i, vir_i )
       IF ( overlap ) EXIT ! jump out of loop
       pot  = pot  + pot_i
       vir  = vir  + vir_i
    END DO

  END SUBROUTINE energy

  SUBROUTINE energy_1 ( ri, i, j_range, box, r_cut, overlap, pot, vir )

    REAL, DIMENSION(3), INTENT(in)  :: ri         ! coordinates of atom of interest
    INTEGER,            INTENT(in)  :: i, j_range ! index, and partner index range
    REAL,               INTENT(in)  :: box        ! simulation box length
    REAL,               INTENT(in)  :: r_cut      ! potential cutoff distance
    LOGICAL,            INTENT(out) :: overlap    ! shows if an overlap was detected
    REAL,               INTENT(out) :: pot, vir   ! potential and virial

    ! Calculates potential energy and virial of atom in ri
    ! with j/=i, j>i, or j<i depending on j_range
    ! Includes a check for overlap (potential too high) to avoid overflow
    ! If overlap==.true., the values of pot and vir should not be used
    ! It is assumed that r has been divided by box
    ! Results are in LJ units where sigma = 1, epsilon = 1

    INTEGER            :: j, jj, nj
    REAL               :: r_cut_box, r_cut_box_sq, box_sq
    REAL               :: sr2, sr6, rij_sq
    REAL, DIMENSION(3) :: rij
    REAL, PARAMETER    :: sr2_overlap = 1.8 ! overlap threshold

    IF ( n > SIZE(r,dim=2) ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Array bounds error for r', n, SIZE(r,dim=2)
       STOP 'Error in energy_1'
    END IF

    r_cut_box    = r_cut / box
    r_cut_box_sq = r_cut_box**2
    box_sq       = box**2

    pot     = 0.0
    vir     = 0.0
    overlap = .FALSE.

    CALL get_neighbours ( i, j_range, nj )

    DO jj = 1, nj
       j = j_list(jj)

       IF ( i == j ) CYCLE

       rij(:) = ri(:) - r(:,j)
       rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundaries in box=1 units
       rij_sq = SUM ( rij**2 )

       IF ( rij_sq < r_cut_box_sq ) THEN

          rij_sq = rij_sq * box_sq ! now in sigma=1 units
          sr2 = 1.0 / rij_sq       ! (sigma/rij)**2

          IF ( sr2 > sr2_overlap ) THEN
             overlap = .TRUE.
             EXIT ! jump out of loop
          END IF

          sr6 = sr2**3
          pot = pot + sr6**2 - sr6
          vir = vir + 2.0 * sr6**2 - sr6

       END IF

    END DO
    pot = 4.0 * pot 
    vir = 24.0 * vir
    vir = vir / 3.0

  END SUBROUTINE energy_1

  SUBROUTINE get_neighbours ( i, op, nj )
    USE link_list_module, ONLY : sc, list, head, c

    ! Arguments
    INTEGER, INTENT(in)  :: i  ! particle whose neighbours are required
    INTEGER, INTENT(in)  :: op ! ne or gt, determining the range of neighbours
    INTEGER, INTENT(out) :: nj ! number of j-partners

    ! Set up vectors to each cell in neighbourhood of 3x3x3 cells in cubic lattice
    INTEGER, PARAMETER :: nk = 13 
    INTEGER, DIMENSION(3,-nk:nk), PARAMETER :: d = RESHAPE( [ &
         &   -1,-1,-1,    0,-1,-1,    1,-1,-1, &
         &   -1, 1,-1,    0, 1,-1,    1, 1,-1, &
         &   -1, 0,-1,    1, 0,-1,    0, 0,-1, &
         &    0,-1, 0,    1,-1, 0,   -1,-1, 0, &
         &   -1, 0, 0,    0, 0, 0,    1, 0, 0, &
         &    1, 1, 0,   -1, 1, 0,    0, 1, 0, &
         &    0, 0, 1,   -1, 0, 1,    1, 0, 1, &
         &   -1,-1, 1,    0,-1, 1,    1,-1, 1, &
         &   -1, 1, 1,    0, 1, 1,    1, 1, 1    ], [ 3, 2*nk+1 ] )

    ! Local variables
    INTEGER :: k1, k2, k, j
    INTEGER, DIMENSION(3) :: cj

    SELECT CASE ( op )
    CASE ( ne ) ! check every atom other than i in all cells
       k1 = -nk
       k2 =  nk
    CASE ( gt ) ! check half neighbour cells and j downlist from i in current cell
       k1 = 0
       k2 = nk
    CASE default ! should never happen
       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Argument op incorrect', op
       STOP 'Error in get_neighbours'
    END SELECT

    nj = 0 ! Will store number of neighbours found

    DO k = k1, k2 ! Begin loop over neighbouring cells

       cj(:) = c(:,i) + d(:,k)      ! Neighbour cell index
       cj(:) = MODULO ( cj(:), sc ) ! Periodic boundary correction

       IF ( k == 0 .AND. op == gt ) THEN
          j = list(i) ! check down-list from i in i-cell
       ELSE
          j = head(cj(1),cj(2),cj(3)) ! check entire j-cell
       END IF

       DO ! Begin loop over j atoms in list

          IF ( j == 0 ) EXIT

          IF ( j /= i ) THEN
             nj         = nj + 1 ! increment count of j atoms
             j_list(nj) = j      ! store new j atom
          END IF
          j = list(j) ! Next atom in j cell

       ENDDO ! End loop over j atoms in list

    ENDDO ! End loop over neighbouring cells 

  END SUBROUTINE get_neighbours

  SUBROUTINE energy_lrc ( n, box, r_cut, pot, vir )
    INTEGER, INTENT(in)  :: n        ! number of atoms
    REAL,    INTENT(in)  :: box      ! simulation box length
    REAL,    INTENT(in)  :: r_cut    ! cutoff distance
    REAL,    INTENT(out) :: pot, vir ! potential and virial

    ! Calculates long-range corrections for Lennard-Jones potential and virial
    ! These are the corrections to the total values
    ! r_cut, box, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL               :: sr3, density
    REAL, PARAMETER    :: pi = 4.0 * ATAN(1.0)

    sr3 = ( 1.0 / r_cut ) ** 3
    pot = (8.0/9.0)  * sr3**3  -(8.0/3.0)  * sr3
    vir = (32.0/9.0) * sr3**3  -(32.0/6.0) * sr3

    density =  REAL(n) / box**3
    pot     = pot * pi * density * real(n)
    vir     = vir * pi * density * real(n)

  END SUBROUTINE energy_lrc

  SUBROUTINE move ( i, ri )
    USE link_list_module, ONLY : c_index, move_in_list
    INTEGER,               INTENT(in) :: i
    REAL,    DIMENSION(3), INTENT(in) :: ri

    INTEGER, DIMENSION(3) :: ci

    r(:,i) = ri                ! New position
    ci(:)  = c_index ( ri(:) ) ! New cell index
    CALL move_in_list ( i, ci(:) )

  END SUBROUTINE move

  SUBROUTINE create ( ri )
    USE link_list_module, ONLY : c_index, create_in_list
    REAL, DIMENSION(3), INTENT(in) :: ri

    INTEGER, DIMENSION(3) :: ci

    n      = n+1               ! increase number of atoms
    r(:,n) = ri(:)             ! add new atom at the end
    ci(:)  = c_index ( ri(:) ) ! New cell index
    CALL create_in_list ( n, ci )

  END SUBROUTINE create

  SUBROUTINE destroy ( i )
    USE link_list_module, ONLY : destroy_in_list, move_in_list, c
    INTEGER, INTENT(in) :: i

    r(:,i) = r(:,n) ! replace atom i coordinates with atom n
    CALL destroy_in_list ( n, c(:,n) )
    CALL move_in_list ( i, c(:,n) )
    n = n - 1  ! reduce number of atoms

  END SUBROUTINE destroy

END MODULE mc_lj_module
