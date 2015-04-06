! mc_zvt_lj.f90 (uses mc_lj_module.f90, utility_module.f90)
! Monte Carlo, zVT (grand) ensemble, Lennard-Jones atoms
PROGRAM mc_zvt_lj
  USE utility_module, ONLY : metropolis, read_cnf_atoms, write_cnf_atoms, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add, random_integer
  USE mc_lj_module,   ONLY : energy_1, energy, energy_lrc, n, r, ne
  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1

  ! Most important variables
  REAL                  :: sigma       ! atomic diameter (in units where box=1)
  REAL                  :: box         ! box length (in units where sigma=1)
  REAL                  :: density     ! reduced density n*sigma**3/box**3
  REAL                  :: dr_max      ! maximum MC displacement
  REAL                  :: temperature ! specified temperature
  REAL                  :: activity    ! specified activity z
  REAL                  :: r_cut       ! potential cutoff distance
  REAL                  :: pot         ! total potential energy
  REAL                  :: vir         ! total virial
  REAL, DIMENSION(-1:1) :: move_ratio  ! acceptance ratio of moves (to be averaged)
  REAL                  :: pressure    ! pressure (LJ sigma=1 units, to be averaged)
  REAL                  :: potential   ! potential energy per atom (LJ sigma=1 units, to be averaged)

  LOGICAL                  :: overlap
  INTEGER                  :: blk, stp, i, nstep, nblock
  INTEGER                  :: try, ntry
  INTEGER, DIMENSION(-1:1) :: tries, moves ! count tries and moves of different kinds
  REAL                     :: pot_old, pot_new, pot_lrc, del_pot, vir_old, vir_new, vir_lrc, del_vir, delta
  REAL                     :: prob_move, prob_create
  REAL, DIMENSION(3)       :: ri   ! position of atom i
  REAL, DIMENSION(3)       :: zeta ! random numbers

  CHARACTER(len=13), PARAMETER :: cnf_prefix = 'mc_zvt_lj.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /run_parameters/ nblock, nstep, temperature, activity, prob_move, r_cut, dr_max

  WRITE(*,'(''mc_zvt_lj'')')
  WRITE(*,'(''Monte Carlo, constant-zVT, Lennard-Jones'')')
  WRITE(*,'(''Results in units epsilon = sigma = 1'')')

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  activity    = 1.0
  prob_move   = 0.34
  r_cut       = 2.5
  dr_max      = 0.15
  READ(*,nml=run_parameters)
  prob_create = (1.0-prob_move)/2.0
  WRITE(*,'(''Number of blocks'',             t40,i15)'  ) nblock
  WRITE(*,'(''Number of steps per block'',    t40,i15)'  ) nstep
  WRITE(*,'(''Temperature'',                  t40,f15.5)') temperature
  WRITE(*,'(''Activity'',                     t40,f15.5)') activity
  WRITE(*,'(''Probability of move'',          t40,f15.5)') prob_move
  WRITE(*,'(''Probability of create/destroy'',t40,f15.5)') prob_create
  WRITE(*,'(''Potential cutoff distance'',    t40,f15.5)') r_cut
  WRITE(*,'(''Maximum displacement'',         t40,f15.5)') dr_max

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box )
  WRITE(*,'(''Number of particles'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma   = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density

  ALLOCATE ( r(3,n*2) ) ! allocate plenty of spare space

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r )

  ! Convert to box units
  r(:,1:n) = r(:,1:n) / box
  r(:,1:n) = r(:,1:n) - ANINT ( r(:,1:n) ) ! Periodic boundaries
  sigma  = sigma / box
  r_cut  = r_cut / box
  dr_max = dr_max / box
  WRITE(*,'(''sigma (in box units)'',t40,f15.5)') sigma
  WRITE(*,'(''r_cut (in box units)'',t40,f15.5)') r_cut
  IF ( r_cut > 0.5 ) STOP 'r_cut too large '

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) STOP 'Overlap in initial configuration'
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(*,'(''Initial potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Initial pressure (sigma units)'',        t40,f15.5)') pressure

  CALL run_begin ( ['Move ratio','+1 ratio  ','-1 ratio  ', &
       &            'Potential ','Pressure  ','Density   '] ) ! must all be character*10 constants

  ntry = n ! each step consists of ntry tries (during which n might vary)

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        tries = 0
        moves = 0

        DO try = 1, ntry ! Begin loop over tries of different kinds

           CALL random_NUMBER ( zeta(1) ) ! uniform random numbers in RANGE (0,1)

           IF ( zeta(1) < prob_move ) THEN ! try particle move

              tries(0) = tries(0) + 1
              CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
              zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

              i = random_integer ( 1, n )
              ri(:) = r(:,i)
              CALL  energy_1 ( ri, i, ne, sigma, r_cut, overlap, pot_old, vir_old )
              IF ( overlap ) STOP 'Overlap in current configuration'
              ri(:) = ri(:) + zeta * dr_max   ! trial move to new position
              ri(:) = ri(:) - ANINT ( ri(:) ) ! periodic boundary correction
              CALL  energy_1 ( ri, i, ne, sigma, r_cut, overlap, pot_new, vir_new )

              IF ( .NOT. overlap ) THEN ! consider non-overlapping configuration
                 delta = ( pot_new - pot_old ) / temperature
                 IF ( metropolis ( delta ) ) THEN    ! accept Metropolis test
                    pot    = pot + pot_new - pot_old ! update potential energy
                    vir    = vir + vir_new - vir_old ! update virial
                    r(:,i) = ri(:)                   ! update position
                    moves(0)  = moves(0) + 1         ! increment move counter
                 END IF ! reject Metropolis test
              END IF ! reject overlapping configuration

           ELSE IF ( zeta(1) < prob_move + prob_create ) THEN ! try create
              tries(1) = tries(1) + 1

              IF ( n+1 > SIZE(r,dim=2) ) CALL resize

              CALL random_NUMBER ( ri ) ! three uniform random numbers in range (0,1)
              ri = ri - 0.5             ! now in range (-0.5,+0.5)
              CALL energy_1 ( ri, n+1, ne, sigma, r_cut, overlap, del_pot, del_vir )
              CALL energy_lrc ( n+1, sigma, r_cut, pot_lrc, vir_lrc ) ! LRC for n+1 atoms
              del_pot = del_pot + pot_lrc
              del_vir = del_vir + vir_lrc
              CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc ) ! LRC for n atoms
              del_pot = del_pot - pot_lrc
              del_vir = del_vir - vir_lrc

              IF ( .NOT. overlap ) THEN ! consider non-overlapping configuration
                 delta = del_pot / temperature - LOG ( activity / REAL ( n+1 ) )
                 IF ( metropolis ( delta ) ) THEN ! accept Metropolis test
                    n        = n+1           ! increase number of atoms
                    r(:,n)   = ri            ! add new atom coordinates
                    pot      = pot + del_pot ! update total potential energy
                    vir      = vir + del_vir ! update total virial
                    moves(1) = moves(1) + 1  ! increment creation move counter
                 END IF ! reject Metropolis test
              END IF ! reject overlapping configuration

           ELSE ! try destroy
              tries(-1) = tries(-1) + 1
              i = random_integer ( 1, n )

              CALL  energy_1 ( r(:,i), i, ne, sigma, r_cut, overlap, del_pot, del_vir )
              CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc ) ! LRC for n atoms
              del_pot = del_pot + pot_lrc
              del_vir = del_vir + vir_lrc
              CALL energy_lrc ( n-1, sigma, r_cut, pot_lrc, vir_lrc ) ! LRC for n-1 atoms
              del_pot = del_pot - pot_lrc
              del_vir = del_vir - vir_lrc
              del_pot = -del_pot ! change sign for a removal
              del_vir = -del_vir ! change sign for a removal

              IF ( overlap ) STOP 'Overlap found on particle removal'
              delta = del_pot/temperature - LOG ( REAL ( n ) / activity )
              IF ( metropolis ( delta ) ) THEN ! accept Metropolis test
                 r(:,i)    = r(:,n)            ! replace atom i with atom n
                 n         = n - 1             ! reduce number of atoms
                 pot       = pot + del_pot     ! update total potential energy
                 vir       = vir + del_vir     ! update total virial
                 moves(-1) = moves(-1) + 1     ! increment destruction move counter
              END IF ! reject Metropolis test

           END IF ! end choice of move type

        END DO ! End loop over tries of different kinds

        ! Calculate all variables for this step
        move_ratio = REAL(moves) / REAL(tries)
        potential  = pot / REAL(n)
        density    = REAL(n)*sigma**3
        pressure   = density * temperature + vir / box**3
        CALL blk_add ( [move_ratio(0),move_ratio(1),move_ratio(-1), &
             &          potential,pressure,density] )

     END DO ! End loop over steps

     CALL blk_end ( blk )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end

  potential = pot / REAL ( n )
  density   = REAL(n)*sigma**3
  pressure  = density * temperature + vir / box**3
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final reduced density'',               t40,f15.5)') density
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) STOP 'Overlap in final configuration'
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(*,'(''Final check'')')
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final reduced density'',               t40,f15.5)') density
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )

  DEALLOCATE ( r )

CONTAINS

  SUBROUTINE resize ! reallocates r array, twice as large
    USE mc_lj_module, ONLY : r
    IMPLICIT NONE
    REAL, DIMENSION(:,:), ALLOCATABLE :: tmp
    INTEGER                           :: n_old, n_new

    n_old = SIZE(r,dim=2)
    n_new = 2*n_old
    WRITE(*,'(a,i5)',advance='no') 'Warning: reallocating r array, from old size = ', n_old
    ALLOCATE ( tmp ( 3, n_new ) ) ! temporary array, new size
    tmp(:,1:n_old) = r            ! copy all elements of r
    DEALLOCATE ( r )              ! deallocate r
    ALLOCATE ( r(3,n_new) )       ! reallocate r
    r(:,1:n_old) = tmp(:,1:n_old) ! copy elements back
    DEALLOCATE ( tmp )            ! cleanup
    WRITE(*,'(a,i5)') ' to new size = ', n_new

  END SUBROUTINE resize

END PROGRAM mc_zvt_lj

