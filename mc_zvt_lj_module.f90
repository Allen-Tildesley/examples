! mc_zvt_lj_module.f90  (used by mc_zvt_lj.f90)
! Monte Carlo simulation, constant-zVT (grand) ensemble, Lennard-Jones atoms
module mc_zvt_lj_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nmax, n, r

  INTEGER                              :: nmax, n ! max and actual number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r       ! positions (3,nmax)

!*****************************************************************
! attempted creations and destructions in grand canonical mc.   **
!                                                               **
! these routines allow for a trial destruction or creation in a **
! grand canonical monte carlo program.                          **
!                                                               **
! principal variables:                                          **
!                                                               **
! integer n                   number of atoms before the trial  **
! integer ntrial              number of atoms during the trial  **
! real    rxnew,rynew,rznew   position for addition of atom     **
! real    v                   potential energy + lrc            **
! real    w                   virial + lrc                      **
! real    deltv               change in energy                  **
! real    deltw               change in virial                  **
! real    temp                reduced temperature               **
! real    z                   absolute activity coefficient     **
! real    sigma               lennard jones diameter            **
! real    rcut                reduced cutoff distance           **
! real    rmin                reduced minimum separation        **
! logical ovrlap              true for substantial atom overlap **
! logical create              true for an accepted creation     **
! logical ghost               true for an accepted destruction  **
!                                                               **
! usage:                                                        **
!                                                               **
! routines in and out should be called with equal probability   **
! in a grand canonical monte carlo simulation. if a trial       **
! creation is accepted then create is set to true. if a trial   **
! destruction is accepted then ghost is set to true. the        **
! routines are written for lennard-jones atoms. the box is of   **
! unit length, all distances are scaled to the box length.      **
! trial inputs which result in a separation of less than        **
! 0.5*sigma are rejected. the long-range corrections are        **
! included in v and w. all accumulators are updated in the main **
! part of the program which is not given here.                  **
!*****************************************************************

  ! create
  
  IF ( n+1 > nmax ) CALL resize

        CALL random_NUMBER ( ri ) ! three uniform random numbers in range (0,1)
        ri = ri - 0.5             ! now in range (-0.5,+0.5)
        CALL energy_1 ( ri, n+1, ne, sigma, r_cut, del_pot, del_vir, overlap )
        ! add change in long range corrections n**2 -> (n+1)**2
        del_pot = del_pot + REAL ( 2*n+1 ) * pot_lrc_on ( sigma, r_cut )
        del_vir = del_vir + REAL ( 2*n+1 ) * vir_lrc_on ( sigma, r_cut )
           IF ( .NOT. overlap ) THEN ! consider non-overlapping configuration
              delta = del_pot / temperature - LOG ( z / REAL ( n+1 ) )
              IF ( metropolis ( delta ) ) THEN    ! accept Metropolis test
                 n = n+1                ! increase number of atoms
                 r(:,n) = ri            ! add new atom coordinates
                 pot = pot + del_pot    ! update total potential energy
                 vir = vir + vir_pot    ! update total virial
                 ncreate = ncreate + 1  ! increment creation move counter
              END IF ! reject Metropolis test
           END IF ! reject overlapping configuration

! destroy

        i = random_integer ( 1, n )

        CALL  energy_1 ( r(:,i), i, ne, sigma, r_cut, del_pot, del_vir, overlap )
! add change in long-range corrections  n**2 -> (n-1)**2
        del_pot = del_pot + REAL ( 2*n-1 ) * pot_lrc_on ( sigma, r_cut )
        del_vir = del_vir + REAL ( 2*n-1 ) * vir_lrc_on ( sigma, r_cut )
! change sign for a removal
        del_pot = -del_pot
        del_vir = -del_vir
        IF ( overlap ) STOP 'Overlap found on particle removal'
        delta = del_pot/temperature - LOG ( REAL ( n ) / z )
        IF ( metropolis ( delta ) ) THEN    ! accept Metropolis test
           r(:,i) = r(:,n) ! replace atom i with atom n
           n = n - 1       ! reduce number of atoms
                 pot = pot + del_pot    ! update total potential energy
                 vir = vir + vir_pot    ! update total virial
                 ndestroy = ndestroy + 1  ! increment destruction move counter
           
              END IF ! reject Metropolis test


SUBROUTINE resize ! reallocates r array, twice as large
REAL,DIMENSION(:,:),ALLOCATABLE :: tmp_arr

WRITE(*,'(a,advance='no')') 'Warning: doubling size of r array, '
ALLOCATE ( tmp_arr ( SIZE(r,dim=1), 2*SIZE(r,dim=2) ) ) ! temporary array, new size
tmp_arr(:,1:nmax)=r ! copy all elements of r
DEALLOCATE(r)
ALLOCATE(r,source=tmp_arr) ! reallocate and copy tmp_arr to r
DEALLOCATE(tmp_arr) ! cleanup
nmax = SIZE(r,dim=2)
WRITE(*,'(a,i5)') 'new nmax = ', nmax

END SUBROUTINE resize


  FUNCTION pot_lrc_on ( sigma, r_cut ) ! Long-range correction to potential per atom / N
    REAL             :: pot_lrc_on     ! Function result in LJ (sigma=1) units
    REAL, INTENT(in) :: sigma, r_cut

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3        = ( sigma / r_cut ) ** 3
    pot_lrc_on = (8.0/9.0) * pi * (sigma**3) * ( sr3**3 - 3.0*sr3 )

  END FUNCTION pot_lrc_on

  FUNCTION vir_lrc_on ( sigma, r_cut ) ! Long-range correction to virial per atom / N
    REAL             :: vir_lrc_on     ! Function result in LJ (sigma=1) units
    REAL, INTENT(in) :: sigma, r_cut

    REAL            :: sr3
    REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

    sr3        = ( sigma / r_cut ) ** 3
    vir_lrc_on = (32.0/9.0) * pi * (sigma**3) * ( sr3**3 - 1.5*sr3 )

  END FUNCTION vir_lrc_on
