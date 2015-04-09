! force_verlet_list.f90
! Lennard-Jones force routine using Verlet neighbour list
SUBROUTINE force ( sigma, r_cut, pot, pot_sh, vir )
  REAL, INTENT(in)  :: sigma, r_cut ! potential parameters
  REAL, INTENT(out) :: pot          ! total potential energy
  REAL, INTENT(out) :: pot_sh       ! potential shifted to zero at cutoff
  REAL, INTENT(out) :: vir          ! virial

  ! This is a drop-in replacement for the Lennard-Jones force routine
  ! provided in md_lj_module.f90 and tt is assumed that the routine
  ! has access to the module variables r, f, n 

  ! Calculates potential (unshifted and shifted), virial and forces
  ! It is assumed that potential parameters and positions are in units where box = 1
  ! The Lennard-Jones energy parameter is taken to be epsilon = 1
  ! Forces are calculated in units where box = 1 and epsilon = 1

  REAL, save :: r_list
  REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: r0, dr
  INTEGER, DIMENSION(:), ALLOCATABLE, save :: point   ! index to neighbour list (n)
  INTEGER, DIMENSION(:), ALLOCATABLE, save :: list    ! Verlet neighbour list (?)

  REAL        rijsq, sr2, sr6, vij, virij, fij
  REAL        sigma_sq, r_cut_sq, r_list_sq
  REAL        rxij, ryij, rzij, fxij, fyij, fzij
  INTEGER     i, j, nlist
  INTEGER     jbeg, jend, jnab

  r_list = r_cut * 1.5 ! somewhat arbitrary
  maxnab = 25*n
  ALLOCATE ( point(n), list(maxnab) )

  r_list_sq = r_list **2 
  sigma_sq  = sigma ** 2
  r_cut_sq  = r_cut ** 2
  skin_sq = ( r_list - r_cut ) ** 2

  dr = r - r0 ! displacement since last update
  dr = dr - ANINT ( dr ) ! periodic boundaries in box=1
  dr_sq_max = MAXVAL ( SUM(dr**2,dim=1) )

  IF ( 4.0*dr_sq_max > skin_sq ) THEN ! List needs updating
     r0 = r ! save current configuration for later checking

     nlist = 0

     DO i = 1, n - 1 ! Begin outer loop over atoms

        point(i) = nlist + 1

        DO j = i + 1, n ! Begin inner loop over partner atoms

           rij(:)  = r(:,i) - r(:,j)
           rij(:)  = rij(:) - ANINT ( rij(:) )
           rij_sq = SUM ( rij**2 )

           IF ( rijsq < r_list_sq ) THEN

              nlist = nlist + 1
              if ( nlist > maxnab )  STOP 'list too small'
              list(nlist) = j

           END IF
        end do ! End inner loop over partner atoms
     end do
     point(n) = nlist + 1
  end if

  ! Use the list to calculate interactions

  f = 0.0
  pot = 0.0
  vir = 0.0
  n_cut = 0
  
  DO i = 1, n - 1 ! Begin outer loop over atoms

     DO jnab = point(i), point(i+1) - 1 ! Begin loop over neighbours if any

           j     = list(jnab)

           rij(:) = r(:,i) - r(:,j)
           rij(:) = rij(:) - ANINT( rij(:) )
           rij_sq = SUM ( rij**2 )

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

end do ! End loop over neighbours

end do ! End loop over atoms

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

           RETURN
        END DO



        SUBROUTINE check ( r_cut, r_list, update )

          COMMON / block1 / rx, ry, rz, fx, fy, fz
          COMMON / block3 / rx0, ry0, rz0

          c    *******************************************************************
          c    ** decides whether the list needs to be reconstructed.           **
          c    **                                                               **
          c    ** principal variables:                                          **
          c    **                                                               **
          c    ** REAL     rx(n),ry(n),rz(n)     atom positions                 **
          c    ** REAL     rx0(n),ry0(n),rz0(n)  coordinates at last update     **
          c    ** REAL     r_list                 radius of verlet list          **
          c    ** REAL     r_cut                  cutoff distance for forces     **
          c    ** REAL     dispmx                largest displacement           **
          c    ** INTEGER  n                     number of atoms                **
          c    ** LOGICAL  update                IF true the list is updated    **
          c    **                                                               **
          c    ** usage:                                                        **
          c    **                                                               **
          c    ** check is called to set update before every CALL to force.     **
          c    *******************************************************************

          INTEGER     n
          PARAMETER ( n = 108 )

          REAL        rx(n), ry(n), rz(n), fx(n), fy(n), fz(n)
          REAL        rx0(n), ry0(n), rz0(n)
          REAL        r_cut, r_list
          LOGICAL     update

          REAL        dispmx
          INTEGER     i

          c    *******************************************************************

          c    ** calculate maximum displacement since last update **

          dispmx = 0.0

          DO 30 i = 1, n

             dispmx = MAX ( ABS ( rx(i) - rx0(i) ), dispmx )
             dispmx = MAX ( ABS ( ry(i) - ry0(i) ), dispmx )
             dispmx = MAX ( ABS ( rz(i) - rz0(i) ), dispmx )

30           CONTINUE

             c    ** a conservative test of the list skin crossing **

             dispmx = 2.0 * SQRT  ( 3.0 * dispmx ** 2 )

             update = ( dispmx .GT. ( r_list - r_cut ) )

             RETURN
          END DO

