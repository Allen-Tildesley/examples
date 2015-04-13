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

CONTAINS

  SUBROUTINE initialize ( r_cut )
    USE verlet_list_module, only : initialize_list
    REAL, INTENT(in) :: r_cut

    ALLOCATE ( r(3,n), v(3,n), f(3,n) )

    CALL initialize_list ( n, r_cut )
  END SUBROUTINE initialize

  SUBROUTINE finalize
    USE verlet_list_module, only : finalize_list
    DEALLOCATE ( r, v, f )
    call finalize_list
  END SUBROUTINE finalize
  
  SUBROUTINE force ( sigma, r_cut, pot, pot_sh, vir )
    USE verlet_list_module, ONLY : point, list, make_list

    REAL, INTENT(in)  :: sigma, r_cut ! potential parameters
    REAL, INTENT(out) :: pot          ! total potential energy
    REAL, INTENT(out) :: pot_sh       ! potential shifted to zero at cutoff
    REAL, INTENT(out) :: vir          ! virial

    ! Calculates potential (unshifted and shifted), virial and forces
    ! It is assumed that potential parameters and positions are in units where box = 1
    ! The Lennard-Jones energy parameter is taken to be epsilon = 1
    ! Forces are calculated in units where box = 1 and epsilon = 1

    INTEGER            :: i, j, k, n_cut
    REAL               :: r_cut_sq, sigma_sq, rij_sq, sr2, sr6, sr12, potij, virij
    REAL, DIMENSION(3) :: rij, fij

    r_cut_sq  = r_cut ** 2
    sigma_sq  = sigma ** 2

    CALL make_list ( n, r )
    
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
