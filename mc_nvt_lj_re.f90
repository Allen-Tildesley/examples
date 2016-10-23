! mc_nvt_lj_re.f90
! Monte Carlo, NVT ensemble, replica exchange
PROGRAM mc_nvt_lj_re

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis, random_translate_vector
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       potential_1, potential, move, n, r, &
       &                       potential_type, OPERATOR(+), OPERATOR(-)
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
  ! but will write "standard output" to a file std##.out where ## is the process rank
  ! Similarly, configurations are read, saved, and written to files named cnf##.inp etc

  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in mc_module

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: density     ! Density
  REAL :: dr_max      ! Maximum MC displacement
  REAL :: temperature ! Specified temperature
  REAL :: r_cut       ! Potential cutoff distance

  ! Quantities to be averaged
  REAL :: m_ratio ! Acceptance ratio of moves
  REAL :: x_ratio ! Acceptance ratio of exchange swaps with higher neighbour
  REAL :: en_c    ! Internal energy per atom for simulated, cut, potential
  REAL :: en      ! Internal energy per atom for full potential with LRC
  REAL :: p_c     ! Pressure for simulated, cut, potential
  REAL :: p       ! Pressure for full potential with LRC

  ! Composite interaction = pot & vir & overlap variables
  TYPE(potential_type) :: total, partial_old, partial_new

  REAL, DIMENSION(:), ALLOCATABLE :: every_temperature, every_beta, every_dr_max

  LOGICAL            :: swap
  INTEGER            :: blk, stp, i, nstep, nblock, moves, swap_interval, swapped, updown, ioerr
  REAL               :: beta, other_beta, other_pot_c, delta, zeta
  REAL, DIMENSION(3) :: ri

  INTEGER                      :: output_unit                                  ! output file unit number 
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=6)             :: cnf_prefix = 'cnf##.', std_prefix = 'std##.' ! will have rank inserted
  CHARACTER(len=3)             :: sav_tag = 'sav'                              ! may be overwritten with block number
  CHARACTER(len=2)             :: m_tag                                        ! will contain rank number

  INTEGER                             :: m             ! MPI process rank (id of this process)
  INTEGER                             :: nproc         ! MPI world size (number of processes)
  INTEGER                             :: error         ! MPI error return
  INTEGER                             :: msg_error     ! MPI error return
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: msg_status    ! MPI status return
  INTEGER, PARAMETER                  :: msg1_id = 999 ! MPI message identifier
  INTEGER, PARAMETER                  :: msg2_id = 888 ! MPI message identifier
  INTEGER, PARAMETER                  :: msg3_id = 777 ! MPI message identifier
  INTEGER, PARAMETER                  :: msg4_id = 666 ! MPI message identifier

  NAMELIST /nml/ nblock, nstep, swap_interval, every_temperature, r_cut, every_dr_max

  CALL MPI_Init(error)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,m,error)
  CALL MPI_Comm_size(MPI_COMM_WORLD,nproc,error)
  IF ( nproc > 100 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Number of processes is too large', nproc
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

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_nvt_lj_re'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT, replica exchange'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Simulation uses cut (but not shifted) potential'
  CALL introduction ( output_unit )
  WRITE( unit=output_unit, fmt='(a,t40,i15)') 'This is process rank',   m
  WRITE( unit=output_unit, fmt='(a,t40,i15)') 'Number of processes is', nproc
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator
  CALL RANDOM_NUMBER ( zeta )
  WRITE( unit=output_unit, fmt='(a,t40,f15.5)') 'First random number', zeta
  WRITE( unit=output_unit, fmt='(a)'          ) 'Should be different for different processes!'

  ALLOCATE ( every_temperature(0:nproc-1), every_beta(0:nproc-1), every_dr_max(0:nproc-1) )

  ! Set sensible defaults for testing
  nblock            = 10
  nstep             = 1000
  swap_interval     = 20
  r_cut             = 2.5
  every_temperature = [ ( 0.7*(1.05)**i,  i = 0, nproc-1 ) ] ! just empirical settings for LJ
  every_dr_max      = [ ( 0.15*(1.05)**i, i = 0, nproc-1 ) ] ! just empirical settings for LJ
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_lj_re'
  END IF
  every_beta  = 1.0 / every_temperature ! All the inverse temperatures
  temperature = every_temperature(m)    ! Temperature for this process
  dr_max      = every_dr_max(m)         ! Max displacement for this process
  beta        = every_beta(m)           ! Inverse temperature for this process

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',               nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',      nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Replica exchange swap interval', swap_interval
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',                    temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance',      r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum displacement',           dr_max

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! first call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density

  CALL allocate_arrays ( box, r_cut ) ! Allocate r

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! second call is to get r

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial energy and overlap check
  total = potential ( box, r_cut )
  IF ( total%overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_lj_re'
  END IF
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Swap ratio', &
       &            'E/N (cut)', 'P (cut)', 'E/N (full)', 'P (full)' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           partial_old = potential_1 ( r(:,i), i, box, r_cut ) ! Old atom potential, virial etc

           IF ( partial_old%overlap ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_nvt_lj_re'
           END IF

           ri(:) = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )                        ! Periodic boundary correction

           partial_new = potential_1 ( ri, i, box, r_cut ) ! New atom potential, virial etc

           IF ( .NOT. partial_new%overlap ) THEN ! Test for non-overlapping configuration

              delta = partial_new%pot_c - partial_old%pot_c ! Use cut (but not shifted) potential
              delta = delta / temperature

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 total = total + partial_new - partial_old ! Update total values
                 CALL move ( i, ri )                       ! Update position
                 moves = moves + 1                         ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End test for overlapping configuration

        END DO ! End loop over atoms

        m_ratio  = REAL(moves) / REAL(n)

        IF ( MOD (stp, swap_interval) == 0 ) THEN ! Test for swap interval

           swapped = 0

           DO updown = 1, 2 ! Loop to look one way then the other

              IF ( MOD(m,2) == MOD(updown,2) ) THEN ! Look up, partner is m+1

                 IF ( m+1 < nproc ) THEN ! Ensure partner exists
                    other_beta = every_beta(m+1) ! We already know the other beta
                    CALL MPI_Recv ( other_pot_c, 1, MPI_REAL, m+1, msg1_id, &
                         &          MPI_COMM_WORLD, msg_status, msg_error ) ! Receive pot_c from other process

                    delta = -(beta - other_beta) * ( total%pot_c - other_pot_c ) ! Delta for Metropolis decision
                    swap  = metropolis ( delta ) ! Decision taken on this process
                    CALL MPI_Send ( swap, 1, MPI_LOGICAL, m+1, msg2_id, &
                         &          MPI_COMM_WORLD, msg_error ) ! Send decision to other process

                    IF ( swap ) THEN ! Exchange configurations
                       CALL MPI_Sendrecv_replace ( r, 3*n, MPI_REAL, m+1, msg3_id, m+1, msg4_id, &
                            &                      MPI_COMM_WORLD, msg_status, msg_error )
                       total   = potential ( box, r_cut ) ! Alternatively, we could exchange this information
                       swapped = 1
                    END IF ! End exchange configurations

                 END IF ! End ensure partner exists

              ELSE ! Look down, partner is m-1

                 IF ( m-1 >= 0 ) THEN ! Ensure partner exists
                    CALL MPI_Send ( total%pot_c, 1, MPI_REAL, m-1, msg1_id, &
                         &          MPI_COMM_WORLD, msg_error ) ! Send pot_c to other process

                    CALL MPI_Recv ( swap, 1, MPI_LOGICAL, m-1, msg2_id, &
                         &          MPI_COMM_WORLD, msg_status, msg_error ) ! Receive decision from other process

                    IF ( swap ) THEN ! Exchange configurations
                       CALL MPI_Sendrecv_replace ( r, 3*n, MPI_REAL, m-1, msg4_id, m-1, msg3_id, &
                            &                      MPI_COMM_WORLD, msg_status, msg_error )
                       total = potential ( box, r_cut ) ! Alternatively, we could exchange this information
                    END IF ! End exchange configurations

                 END IF ! End ensure partner exists

              END IF ! End choice of which way to look

           END DO ! End loop to look one way then the other

        END IF ! End test for swap interval

        x_ratio = REAL(swapped*swap_interval) ! Include factor for only attempting at intervals

        ! Calculate all variables for this step
        CALL calculate ( )
        CALL blk_add ( [m_ratio,x_ratio,en_c,p_c,en,p] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk        ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  ! Double-check book-keeping for totals, and overlap
  total = potential ( box, r_cut )
  IF ( total%overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_lj_re'
  END IF
  CALL calculate ( 'Final check' )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  DEALLOCATE ( every_temperature, every_beta, every_dr_max )
  CALL conclusion ( output_unit )

  CLOSE ( unit=output_unit )
  CALL MPI_Finalize(error)

CONTAINS

  SUBROUTINE calculate ( string )
    USE mc_module, ONLY : potential_lrc, pressure_lrc, pressure_delta
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates the properties of interest from total values
    ! and optionally writes them out (e.g. at the start and end of the run)
    ! In this example we simulate using the cut (but not shifted) potential
    ! The values of < p_c >,  < en_c > and density should be consistent (for this potential)
    ! For comparison, long-range corrections are also applied to give
    ! estimates of < en > and < p > for the full (uncut) potential
    ! The value of the cut-and-shifted potential pot_s is not used, in this example

    en_c = total%pot_c / REAL ( n )                ! PE/N for cut (but not shifted) potential
    en_c = en_c + 1.5 * temperature                ! Add ideal gas contribution KE/N to give E_c/N
    en   = en_c + potential_lrc ( density, r_cut ) ! Add long-range contribution to give E/N estimate
    p_c  = total%vir / box**3                      ! Virial contribution to P_c
    p_c  = p_c + density * temperature             ! Add ideal gas contribution to P_c
    p    = p_c + pressure_lrc ( density, r_cut )   ! Add long-range contribution to give P
    p_c  = p_c + pressure_delta ( density, r_cut ) ! Add delta correction to P_c (not needed for P)

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (cut)',  en_c
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (cut)',    p_c
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (full)', en
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (full)',   p
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_nvt_lj_re

