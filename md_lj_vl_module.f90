! md_lj_vl_module.f90 (used by md_nve_lj.f90)
! Molecular dynamics simulation, Lennard-Jones atoms, Verlet neighbour list
! This is a drop-in replacement for md_lj_module in md_lj_module.f90
MODULE md_lj_module

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n, r, v, f
  PUBLIC :: initialize, finalize, force, energy_lrc

  INTEGER                              :: n ! number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r ! positions (3,:)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: v ! velocities (3,:)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f ! forces (3,:)

  REAL,    DIMENSION(:,:), ALLOCATABLE :: r_save ! saved positions for list (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: dr     ! displacements (3,n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: point  ! index to neighbour list (n)
  INTEGER, DIMENSION(:),   ALLOCATABLE :: list   ! Verlet neighbour list

CONTAINS

  SUBROUTINE initialize
    ALLOCATE ( r(3,n), v(3,n), f(3,n) )

    ! Size allocated to list is just a guess
    ALLOCATE ( point(n), list(25*n) )
  END SUBROUTINE initialize

  SUBROUTINE finalize
    DEALLOCATE ( r, v, f )
    DEALLOCATE ( point, list )
  END SUBROUTINE finalize

  SUBROUTINE resize ! reallocates list array, somewhat larger
    
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp
    INTEGER                            :: n_old, n_new

    n_old = SIZE(list)
    n_new = CEILING( 1.25 * REAL(n_old) )
    WRITE(*,'(a,i5,a,i5)') 'Warning: reallocating list array, from old size = ', n_old, ' to ', n_new

    ALLOCATE ( tmp(n_new) ) ! new size for list
    tmp(1:n_old) = list(:)  ! copy elements across

    CALL move_ALLOC ( tmp, list )

  END SUBROUTINE resize
  
  SUBROUTINE force ( sigma, r_cut, pot, pot_sh, vir )
    REAL, INTENT(in)  :: sigma, r_cut ! potential parameters
    REAL, INTENT(out) :: pot          ! total potential energy
    REAL, INTENT(out) :: pot_sh       ! potential shifted to zero at cutoff
    REAL, INTENT(out) :: vir          ! virial

    ! Calculates potential (unshifted and shifted), virial and forces
    ! It is assumed that potential parameters and positions are in units where box = 1
    ! The Lennard-Jones energy parameter is taken to be epsilon = 1
    ! Forces are calculated in units where box = 1 and epsilon = 1

    INTEGER            :: i, j, k, n_cut
    REAL               :: r_skin, r_list, dr_sq_max
    REAL               :: r_cut_sq, sigma_sq, rij_sq, sr2, sr6, sr12, potij, virij
    REAL, DIMENSION(3) :: rij, fij

    LOGICAL, SAVE :: first_call = .TRUE.

    r_skin = sigma * 1.5 ! somewhat arbitrary
    r_list = r_cut + r_skin

    r_cut_sq  = r_cut ** 2
    sigma_sq  = sigma ** 2

    IF ( first_call ) THEN
       CALL make_list ( r_list )
       r_save     = r
       first_call = .FALSE.
    ELSE
       dr = r - r_save                         ! displacement since last list update
       dr = dr - ANINT ( dr )                  ! periodic boundaries in box=1
       dr_sq_max = MAXVAL ( SUM(dr**2,dim=1) ) ! squared maximum displacement
       IF ( 4.0*dr_sq_max > r_skin ** 2 ) THEN
          CALL make_list ( r_list )
          r_save = r
       END IF
    END IF

    f     = 0.0
    pot   = 0.0
    vir   = 0.0
    n_cut = 0

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO k =  point(i), point(i+1) - 1 ! Begin inner loop over neighbour atoms (if any)

          j = list(k)

          rij(:) = r(:,i) - r(:,j)           ! separation vector
          rij(:) = rij(:) - ANINT ( rij(:) ) ! periodic boundary conditions in box=1 units
          rij_sq = SUM ( rij**2 )            ! squared separation

          IF ( rij_sq < r_cut_sq ) THEN

             sr2    = sigma_sq / rij_sq
             sr6    = sr2 ** 3
             sr12   = sr6 ** 2
             potij  = sr12 - sr6
             virij  = potij + sr12
             pot    = pot + potij
             vir    = vir + virij
             fij    = rij * virij / rij_sq
             f(:,i) = f(:,i) + fij
             f(:,j) = f(:,j) - fij
             n_cut  = n_cut + 1

          END IF

       END DO ! End inner loop over neighbour atoms (if any)

    END DO ! End outer loop over atoms

    ! calculate shifted potential
    sr2    = sigma_sq / r_cut_sq
    sr6    = sr2 ** 3
    sr12   = sr6 **2
    potij  = sr12 - sr6
    pot_sh = pot - REAL ( n_cut ) * potij

    ! multiply results by numerical factors

    f      = f      * 24.0
    pot    = pot    * 4.0
    pot_sh = pot_sh * 4.0
    vir    = vir    * 24.0 / 3.0

  END SUBROUTINE force

  SUBROUTINE make_list ( r_list )
    REAL, INTENT(in) :: r_list

    INTEGER            :: i, j, k
    REAL               :: r_list_sq, rij_sq
    REAL, DIMENSION(3) :: rij

    k = 0
    r_list_sq = r_list ** 2

    DO i = 1, n - 1 ! Begin outer loop over atoms

       point(i) = k + 1

       DO j = i + 1, n ! Begin inner loop over partner atoms

          rij(:)  = r(:,i) - r(:,j)
          rij(:)  = rij(:) - ANINT ( rij(:) )
          rij_sq = SUM ( rij**2 )

          IF ( rij_sq < r_list_sq ) THEN

             k = k + 1
             IF ( k > SIZE(list) ) CALL resize
             list(k) = j

          END IF

       END DO ! End inner loop over partner atoms

    END DO ! End outer loop over atoms

    point(n) = k + 1

  END SUBROUTINE make_list

  SUBROUTINE energy_lrc ( n, sigma, r_cut, pot, vir )
    INTEGER, INTENT(in)  :: n            ! number of atoms
    REAL,    INTENT(in)  :: sigma, r_cut ! LJ potential parameters
    REAL,    INTENT(out) :: pot, vir     ! potential and virial

    ! Calculates long-range corrections for Lennard-Jones potential and virial
    ! These are the corrections to the total values
    ! It is assumed that sigma and r_cut are in units where box = 1
    ! Results are in LJ units where sigma = 1, epsilon = 1
    
    REAL               :: sr3, density
    REAL, DIMENSION(2) :: pot2_lrc, vir2_lrc
    REAL, PARAMETER    :: pi = 4.0 * ATAN(1.0)

    sr3         = ( sigma / r_cut ) ** 3
    density     =  REAL(n)*sigma**3
    pot2_lrc(1) =  REAL(n)*(8.0/9.0)  * pi * density * sr3**3 ! LJ12 term
    pot2_lrc(2) = -REAL(n)*(8.0/3.0)  * pi * density * sr3    ! LJ6  term
    vir2_lrc(1) =  REAL(n)*(32.0/9.0) * pi * density * sr3**3 ! LJ12 term
    vir2_lrc(2) = -REAL(n)*(32.0/6.0) * pi * density * sr3    ! LJ6  term
    pot = SUM ( pot2_lrc )
    vir = SUM ( vir2_lrc )

  END SUBROUTINE energy_lrc

END MODULE md_lj_module
