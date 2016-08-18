! mc_nvt_lj_re.f90
! Monte Carlo, NVT ensemble, Lennard-Jones, replica exchange
PROGRAM mc_nvt_lj_re

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE utility_module,   ONLY : metropolis
  USE mc_lj_module,     ONLY : initialize, finalize, energy_1, energy, energy_lrc, move, &
       &                       n, r, ne
  USE mpi

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists
  ! Uses MPI to run replica exchange of configurations with neighbouring temperatures
  ! Assume at most 100 processes numbered 0 to 99

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults
  ! Note that the various processes running this program will share standard input and error units
  ! but will write "standard" output to a file std##.out where ## is the process rank
  ! Similarly, configurations are read, saved, and written to files named cnf##.inp etc

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

  REAL, DIMENSION(:), ALLOCATABLE :: every_temperature, every_beta, every_dr_max

  LOGICAL            :: overlap, swap
  INTEGER            :: blk, stp, i, nstep, nblock, moves, swap_interval, swapped, updown, ioerr
  REAL               :: pot_old, pot_new, pot_lrc, vir_old, vir_new, vir_lrc, delta
  REAL               :: beta, other_beta, other_pot
  REAL, DIMENSION(3) :: ri   ! position of atom i
  REAL, DIMENSION(3) :: zeta ! random numbers

  INTEGER                      :: output_unit                                  ! output file unit number 
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=6)             :: cnf_prefix = 'cnf##.', std_prefix = 'std##.' ! will have rank inserted
  CHARACTER(len=3)             :: sav_tag = 'sav'                              ! may be overwritten with block number
  CHARACTER(len=2)             :: m_tag                                        ! will contain rank number

  INTEGER                             :: m             ! MPI process rank (id of this process)
  INTEGER                             :: p             ! MPI world size (number of processes)
  INTEGER                             :: error         ! MPI error return
  INTEGER                             :: msg_error     ! MPI error return
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: msg_status    ! MPI status return
  INTEGER, PARAMETER                  :: msg1_id = 999 ! MPI message identifier
  INTEGER, PARAMETER                  :: msg2_id = 888 ! MPI message identifier
  INTEGER, PARAMETER                  :: msg3_id = 777 ! MPI message identifier

  NAMELIST /nml/ nblock, nstep, swap_interval, every_temperature, r_cut, every_dr_max

  CALL MPI_Init(error)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,m,error)
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,error)
  IF ( p > 100 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Number of proccesses is too large', p
     STOP 'Error in mc_nvt_lj_re'
  END IF
  WRITE(m_tag,fmt='(i2.2)') m ! convert rank into character form
  cnf_prefix(4:5) = m_tag     ! insert process rank into configuration filename
  std_prefix(4:5) = m_tag     ! insert process rank into standard output filename

  OPEN ( newunit=output_unit, file = std_prefix // out_tag, iostat=ioerr ) ! open file for standard output
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,a)') 'Could not open file ', std_prefix // out_tag
     STOP 'Error in mc_nvt_lj_re'
  END IF

  WRITE( unit=output_unit, fmt='(a)'        ) 'mc_nvt_lj_re'
  WRITE( unit=output_unit, fmt='(a)'        ) 'Monte Carlo, constant-NVT, Lennard-Jones, replica exchange'
  WRITE( unit=output_unit, fmt='(a)'        ) 'Results in units epsilon = sigma = 1'
  WRITE( unit=output_unit, fmt='(a,t40,i15)') 'This is process rank',   m
  WRITE( unit=output_unit, fmt='(a,t40,i15)') 'Number of processes is', p
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator
  CALL RANDOM_NUMBER ( zeta )
  WRITE( unit=output_unit, fmt='(a,t40,3f15.5)') 'First three random numbers', zeta

  ALLOCATE ( every_temperature(0:p-1), every_beta(0:p-1), every_dr_max(0:p-1) )

  ! Set sensible defaults for testing
  nblock            = 10
  nstep             = 1000
  swap_interval     = 20
  r_cut             = 2.5
  every_temperature = [ ( 0.7*(1.05)**i, i = 0, p-1 ) ]
  every_dr_max      = [ ( 0.15*(1.05)**i, i = 0, p-1 ) ]
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_lj_re'
  END IF
  every_beta  = 1.0 / every_temperature ! all the inverse temperatures
  temperature = every_temperature(m)    ! temperature for this process
  dr_max      = every_dr_max(m)         ! max displacement for this process
  beta        = every_beta(m)           ! inverse temperature for this process

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',               nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',      nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Replica exchange swap interval', swap_interval
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',                    temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance',      r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum displacement',           dr_max

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! first call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box (in sigma units)', box
  sigma   = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Reduced density', density

  ! Convert run and potential parameters to box units
  sigma  = sigma / box
  r_cut  = r_cut / box
  dr_max = dr_max / box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'sigma (in box units)', sigma
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'r_cut (in box units)', r_cut
  IF ( r_cut > 0.5 ) THEN
     WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_cut too large ', r_cut
     STOP 'Error in mc_nvt_lj_re'
  END IF

  CALL initialize ( r_cut ) ! Allocate r

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! second call is to get r

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_lj_re'
  END IF
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial potential energy (sigma units)', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial pressure (sigma units)',         pressure

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
           IF ( overlap ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_nvt_lj_re'
           END IF
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

        IF ( MOD (stp, swap_interval) == 0 ) THEN ! test for swap interval
           swapped = 0
           DO updown = 1, 2 ! loop to look one way then the other
              IF ( MOD(m,2) == MOD(updown,2) ) THEN ! look up, partner is m+1
                 IF ( m+1 < p ) THEN ! ensure partner exists
                    other_beta = every_beta(m+1)
                    CALL MPI_Sendrecv ( pot,       1, MPI_REAL, m+1, msg1_id, &
                         &              other_pot, 1, MPI_REAL, m+1, msg2_id, &
                         &              MPI_COMM_WORLD, msg_status, msg_error )
                    delta = -(beta - other_beta) * (pot - other_pot)
                    swap  = metropolis ( delta )
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
                    other_beta = every_beta(m-1)
                    CALL MPI_Sendrecv (  pot,       1, MPI_REAL, m-1, msg2_id, &
                         &               other_pot, 1, MPI_REAL, m-1, msg1_id, &
                         &               MPI_COMM_WORLD, msg_status, msg_error )
                    CALL MPI_Recv ( swap, 1, MPI_LOGICAL, m-1, msg3_id, &
                         &          MPI_COMM_WORLD, msg_status, msg_error )
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
        swap_ratio = REAL(swapped*swap_interval) ! factor for only attempting at intervals
        potential  = pot / REAL(n)
        pressure   = density * temperature + vir / box**3
        CALL blk_add ( [move_ratio,swap_ratio,potential,pressure] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk        ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final potential energy (sigma units)', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final pressure (sigma units)',         pressure

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_lj_re'
  END IF
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a)'           ) 'Final check'
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final potential energy (sigma units)', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final pressure (sigma units)',         pressure

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )
  CALL time_stamp ( output_unit )

  CALL finalize
  DEALLOCATE ( every_temperature, every_beta, every_dr_max )

  CLOSE ( unit=output_unit )
  CALL MPI_Finalize(error)

END PROGRAM mc_nvt_lj_re

