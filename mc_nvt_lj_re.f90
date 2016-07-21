! mc_nvt_lj_re.f90
! Monte Carlo, NVT ensemble, Lennard-Jones, replica exchange
PROGRAM mc_nvt_lj_re
  USE utility_module, ONLY : metropolis, read_cnf_atoms, write_cnf_atoms, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_lj_module,   ONLY : initialize, finalize, energy_1, energy, energy_lrc, move, &
       &                     n, r, ne
  USE mpi
  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists
  ! Exchanges configurations with neighbouring temperatures

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1

  ! Most important variables
  REAL :: sigma       ! atomic diameter (in units where box=1)
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: dr_max      ! maximum MC displacement
  REAL :: temperature ! specified temperature
  REAL :: r_cut       ! potential cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: vir         ! total virial
  REAL :: move_ratio  ! acceptance ratio of moves (to be averaged)
  REAL :: swap_ratio  ! acceptance ratio of swaps with higher neighbour (to be averaged)
  REAL :: pressure    ! pressure (LJ sigma=1 units, to be averaged)
  REAL :: potential   ! potential energy per atom (LJ sigma=1 units, to be averaged)

  REAL, DIMENSION(:), ALLOCATABLE :: every_temperature, every_dr_max

  LOGICAL            :: overlap, swap
  INTEGER            :: blk, stp, i, nstep, nblock, moves, swap_interval, swapped, updown
  REAL               :: pot_old, pot_new, pot_lrc, vir_old, vir_new, vir_lrc, delta
  REAL               :: beta, other_beta, other_pot, ran
  REAL, DIMENSION(3) :: ri   ! position of atom i
  REAL, DIMENSION(3) :: zeta ! random numbers

  INTEGER :: log_unit ! log file unit number 
  CHARACTER(len=3),  PARAMETER :: run_prefix = 'mc_'
  CHARACTER(len=7),  PARAMETER :: cnfinp_tag = '.cnfinp', cnfout_tag = '.cnfout'
  CHARACTER(len=7)             :: cnfsav_tag = '.cnfsav' ! may be overwritten with block number
  CHARACTER(len=3)             :: rank_tag = 'xxx' ! will be overwritten with processor rank

  INTEGER                             :: m      ! MPI processor rank (id of this process)
  INTEGER                             :: p      ! MPI world size (number of processes)
  INTEGER                             :: error  ! MPI error return
  INTEGER                             :: msg_error  ! MPI error return
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: msg_status ! MPI status return
  INTEGER, PARAMETER                  :: msg1_id = 999, msg2_id = 888, msg3_id = 777

  NAMELIST /params/ nblock, nstep, swap_interval, every_temperature, r_cut, every_dr_max

  CALL MPI_Init(error)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,m,error)
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,error)
  WRITE(rank_tag,'(i3.3)') m

  OPEN(newunit=log_unit, file = run_prefix // rank_tag // '.log' )
  WRITE(log_unit,'(''mc_nvt_lj_re'')')
  WRITE(log_unit,'(''Monte Carlo, constant-NVT, Lennard-Jones, replica exchange'')')
  WRITE(log_unit,'(''Results in units epsilon = sigma = 1'')')
  WRITE(log_unit,'(a,t40,i15)') 'This is processor rank', m
  WRITE(log_unit,'(a,t40,i15)') 'Number of processors is', p

  CALL RANDOM_SEED () ! Initialize random number generator
  CALL random_NUMBER ( zeta )
  WRITE(log_unit,'(a,t40,3f15.5)') 'First three random numbers', zeta

  ALLOCATE ( every_temperature(0:p-1), every_dr_max(0:p-1) )

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  swap_interval = 20
  r_cut       = 2.5
  every_temperature = [ ( 0.7*(1.05)**i, i = 0, p-1 ) ]
  every_dr_max      = [ ( 0.15*(1.05)**i, i = 0, p-1 ) ]
  READ(*,nml=params)
  temperature = every_temperature(m) ! temperature for this process
  dr_max = every_dr_max(m) ! max displacement for this process
  beta = 1.0 / temperature

  WRITE(log_unit,'(''Number of blocks'',         t40,i15)'  ) nblock
  WRITE(log_unit,'(''Number of steps per block'',t40,i15)'  ) nstep
  WRITE(log_unit,'(''Replica exchange swap interval'',t40,i15)'  ) swap_interval
  WRITE(log_unit,'(''Temperature'',              t40,f15.5)') temperature
  WRITE(log_unit,'(''Potential cutoff distance'',t40,f15.5)') r_cut
  WRITE(log_unit,'(''Maximum displacement'',     t40,f15.5)') dr_max

  CALL read_cnf_atoms ( run_prefix//rank_tag//cnfinp_tag, n, box )
  WRITE(log_unit,'(''Number of particles'', t40,i15)'  ) n
  WRITE(log_unit,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE(log_unit,'(''Reduced density'',t40,f15.5)') density

  ! Convert run and potential parameters to box units
  sigma  = sigma / box
  r_cut  = r_cut / box
  dr_max = dr_max / box
  WRITE(log_unit,'(''sigma (in box units)'',t40,f15.5)') sigma
  WRITE(log_unit,'(''r_cut (in box units)'',t40,f15.5)') r_cut
  IF ( r_cut > 0.5 ) STOP 'r_cut too large '

  CALL initialize ( r_cut ) ! Allocate r

  CALL read_cnf_atoms ( run_prefix//rank_tag//cnfinp_tag, n, box, r )

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) STOP 'Overlap in initial configuration'
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(log_unit,'(''Initial potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(log_unit,'(''Initial pressure (sigma units)'',        t40,f15.5)') pressure

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Swap ratio', 'Potential', 'Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

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
                 CALL move ( i, ri )              ! update position
                 moves  = moves + 1               ! increment move counter
              END IF ! reject Metropolis test
           END IF ! reject overlapping configuration

        END DO ! End loop over atoms

        IF ( MOD(stp, swap_interval) == 0 ) THEN ! test for swap interval
           swapped = 0
           DO updown = 1, 2 ! loop to look one way then the other
              IF ( MOD(m,2) == MOD(updown,2) ) THEN ! look up, partner is m+1
                 IF ( m+1 < p ) THEN ! ensure partner exists
                    other_beta = 1.0 / every_temperature(m+1)
                    CALL MPI_Sendrecv ( pot,       1, MPI_REAL, m+1, msg1_id, &
                         &              other_pot, 1, MPI_REAL, m+1, msg2_id, &
                         &              MPI_COMM_WORLD, msg_status, msg_error )
                    delta = -(beta - other_beta) * (pot - other_pot)
                    CALL random_NUMBER(ran)
                    swap = EXP(-delta) > ran
                    CALL MPI_Send ( swap, 1, MPI_LOGICAL, m+1, msg3_id, &
                         &          MPI_COMM_WORLD, msg_error )
                    IF ( swap ) THEN ! exchange configurations
                       CALL MPI_Sendrecv_replace ( r, 3*n, MPI_REAL, m+1, msg1_id, m+1, msg2_id, &
                            &                      MPI_COMM_WORLD, msg_status, msg_error )
                       pot = other_pot
                       swapped = 1
                    END IF ! end exchange configurations
                 END IF ! end ensure partner exists
              ELSE ! look down, partner is m-1
                 IF ( m-1 >= 0 ) THEN ! ensure partner exists
                    other_beta = 1.0 / every_temperature(m-1)
                    CALL MPI_Sendrecv (  pot,       1, MPI_REAL, m-1, msg2_id, &
                         &               other_pot, 1, MPI_REAL, m-1, msg1_id, &
                         &               MPI_COMM_WORLD, msg_status, msg_error )
                    CALL MPI_Recv ( swap, 1, MPI_LOGICAL, m-1, msg3_id, &
                         &  MPI_COMM_WORLD, msg_status, msg_error )
                    IF ( swap ) THEN ! exchange configurations
                       CALL MPI_Sendrecv_replace ( r, 3*n, MPI_REAL, m-1, msg2_id, m-1, msg1_id, &
                            &                      MPI_COMM_WORLD,msg_status,msg_error)
                       pot = other_pot
                    END IF ! end exchange configurations
                 END IF ! end ensure partner exists
              END IF ! end choice of which way to look
           END DO ! end loop to look one way then the other
        END IF ! end test for swap interval

        ! Calculate all variables for this step
        move_ratio = REAL(moves) / REAL(n)
        swap_ratio = REAL(swapped*swap_interval) ! correct for only attempting at intervals
        potential  = pot / REAL(n)
        pressure   = density * temperature + vir / box**3
        CALL blk_add ( [move_ratio,swap_ratio,potential,pressure] )

     END DO ! End loop over steps

     CALL blk_end ( blk )
     IF ( nblock < 1000 ) WRITE(cnfsav_tag(5:7),'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( run_prefix//rank_tag//cnfsav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end

  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(log_unit,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(log_unit,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) STOP 'Overlap in final configuration'
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(log_unit,'(''Final check'')')
  WRITE(log_unit,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(log_unit,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL write_cnf_atoms ( run_prefix//rank_tag//cnfout_tag, n, box, r*box )

  CALL finalize

  CLOSE ( log_unit )
  CALL MPI_Finalize(error)

END PROGRAM mc_nvt_lj_re

