! md_lj_llle_module.f90 (uses link_list_module.f90)
! Link-list algorithm for Lees-Edwards boundaries
MODULE md_lj_le_module

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
    USE link_list_module, ONLY : initialize_list
    REAL, INTENT(in) :: r_cut

    ALLOCATE ( r(3,n), v(3,n), f(3,n) )
    CALL initialize_list ( n, r_cut )

  END SUBROUTINE initialize

  SUBROUTINE finalize
    USE link_list_module, ONLY : finalize_list
    DEALLOCATE ( r, v, f )
    CALL finalize_list
  END SUBROUTINE finalize

  SUBROUTINE force ( sigma, r_cut, strain, pot, pot_sh, vir )
    USE link_list_module, ONLY : make_list, nc, head, list

    REAL, INTENT(in)  :: sigma, r_cut ! potential parameters
    REAL, intent(in)  :: strain       ! shear strain
    REAL, INTENT(out) :: pot          ! total potential energy
    REAL, INTENT(out) :: pot_sh       ! potential shifted to zero at cutoff
    REAL, INTENT(out) :: vir          ! virial

    ! Calculates potential (unshifted and shifted), virial and forces
    ! It is assumed that potential parameters and positions are in units where box = 1
    ! The Lennard-Jones energy parameter is taken to be epsilon = 1
    ! Forces are calculated in units where box = 1 and epsilon = 1
    ! Uses link lists
    ! Lees-Edwards boundaries, in sliding brick arrangement
    ! Flow direction is x (=1)
    ! Gradient direction is y (=2)
    ! Vorticity direction is z (=3)

    INTEGER               :: i, j, n_cut, ci1, ci2, ci3, k, k_max, shift
    INTEGER, DIMENSION(3) :: ci, cj
    REAL                  :: r_cut_sq, sigma_sq, rij_sq, sr2, sr6, sr12, potij, virij
    REAL,    DIMENSION(3) :: rij, fij

    ! Set up vectors to half the cells in neighbourhood of 3x3x3 cells in cubic lattice
    ! The cells are chosen so that if (d1,d2,d3) appears, then (-d1,-d2,-d3) does not.
    ! For sheared boundaries, with the chosen axes, it is most convenient to make this set
    ! include all the cells with d2=+1 and none of the ones with d2=-1
    INTEGER, PARAMETER :: nk = 13, nk_extra = nk+3
    INTEGER, DIMENSION(3,0:nk_extra), PARAMETER :: d = RESHAPE( [ &
         &    0, 0, 0,    1, 0, 0,    1, 0, 1,   -1, 0, 1,  0, 0, 1, & ! 5 cells with d2=0
         &    1, 1,-1,    1, 1, 0,    1, 1, 1,   & ! 3 cells with d1= 1, d2=1
         &    0, 1,-1,    0, 1, 0,    0, 1, 1,   & ! 3 cells with d1= 0, d2=1
         &   -1, 1,-1,   -1, 1, 0,   -1, 1, 1,   & ! 3 cells with d1=-1, d2=1
         &   -2, 1,-1,   -2, 1, 0,   -2, 1, 1 ], & ! 3 cells with d1=-2, d2=1 (extra cells, for shear)
         &  [ 3, nk_extra+1 ] )
    INTEGER, DIMENSION(3,0:nk_extra) :: dd ! will hold sheared cell indices

    r_cut_sq = r_cut ** 2
    sigma_sq = sigma ** 2

    f     = 0.0
    pot   = 0.0
    vir   = 0.0
    n_cut = 0

    ! Lees-Edwards periodic boundaries
    r(1,:) = r(1,:) - ANINT ( r(2,:) ) * strain ! Extra correction
    r(:,:) = r(:,:) - ANINT ( r(:,:) )          ! Standard correction
    CALL make_list ( n, r )

    shift  = floor ( strain * REAL ( nc ) ) ! strain measured in cell lengths

    ! Triple loop over cells
    DO ci1 = 0, nc-1
       DO ci2 = 0, nc-1
          DO ci3 = 0, nc-1
             ci(:) = [ ci1, ci2, ci3 ]
             i = head(ci1,ci2,ci3)

             ! Set up correct neighbour cell indices
             IF ( ci2 == nc-1 ) THEN                       ! top layer
                dd(:,0:4)        = d(:,0:4)                ! five cells do not need adjustment
                dd(:,5:nk_extra) = d(:,5:nk_extra) - shift ! remaining cells need adjustment
                k_max = nk_extra                           ! extra cells need to be checked
             ELSE                                          ! not top layer
                dd(:,0:nk) = d(:,0:nk)                     ! standard list copied
                k_max = nk                                 ! no extra cells need checking
             END IF
             
             DO ! Begin loop over i atoms in list
                IF ( i == 0 ) EXIT ! end of link list

                ! Loop over neighbouring cells
                DO k = 0, k_max

                   IF ( k == 0 ) THEN
                      j = list(i) ! Look downlist from i in current cell
                   ELSE
                      cj(:) = ci(:) + dd(:,k)         ! Neighbour cell index
                      cj(:) = MODULO ( cj(:), nc )    ! Periodic boundary correction
                      j     = head(cj(1),cj(2),cj(3)) ! Look at all atoms in neighbour cell
                   END IF

                   DO ! Begin loop over j atoms in list
                      IF ( j == 0 ) EXIT ! end of link list
                      IF ( j == i ) STOP 'Index error in force' ! This should never happen

                      rij(:) = r(:,i) - r(:,j)                    ! separation vector
                      rij(1) = rij(1) - ANINT ( rij(2) ) * strain ! Extra correction
                      rij(:) = rij(:) - ANINT ( rij(:) )          ! periodic boundary conditions
                      rij_sq = SUM ( rij**2 )                     ! squared separation

                      IF ( rij_sq < r_cut_sq ) THEN ! check within cutoff

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

                      END IF ! end check within cutoff

                      j = list(j) ! Next atom in j cell
                   END DO         ! End loop over j atoms in list

                END DO
                ! End loop over neighbouring cells

                i = list(i) ! Next atom in i cell
             END DO         ! End loop over i atoms in list

          END DO
       END DO
    END DO
    ! End triple loop over cells

    ! Calculate shifted potential
    sr2    = sigma_sq / r_cut_sq
    sr6    = sr2 ** 3
    sr12   = sr6 **2
    potij  = sr12 - sr6
    pot_sh = pot - REAL ( n_cut ) * potij

    ! Multiply results by numerical factors
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

END MODULE md_lj_le_module

