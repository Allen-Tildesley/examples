! mc_zvt_lj.f90
! Monte Carlo, zVT (grand) ensemble
PROGRAM mc_zvt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis, random_integer
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       resize, energy_1, energy, energy_lrc, &
       &                       move, create, destroy, n, r, potovr

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts grand canonical Monte Carlo at the given temperature and activity
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Note that long-range corrections are not included in the acceptance/rejection
  ! of creation and destruction moves, but are added to the simulation averages

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in mc_module

  ! Most important variables
  REAL :: box         ! box length
  REAL :: density     ! density
  REAL :: dr_max      ! maximum MC displacement
  REAL :: temperature ! specified temperature
  REAL :: activity    ! specified activity z
  REAL :: r_cut       ! potential cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: vir         ! total virial
  REAL :: m_ratio     ! acceptance ratio of moves (to be averaged)
  REAL :: c_ratio     ! acceptance ratio of creation attempts (to be averaged)
  REAL :: d_ratio     ! acceptance ratio of destruction attempts (to be averaged)
  REAL :: pres_virial ! virial pressure (to be averaged)
  REAL :: potential   ! potential energy per atom (to be averaged)

  TYPE(potovr) :: eng_old, eng_new ! Composite energy = pot & vir & overlap variables

  INTEGER            :: blk, stp, i, nstep, nblock
  INTEGER            :: try, ntry, ioerr
  INTEGER            :: m_tries, m_moves ! count tries and moves
  INTEGER            :: c_tries, c_moves ! count tries and moves for creation
  INTEGER            :: d_tries, d_moves ! count tries and moves for destruction
  REAL               :: pot_lrc, vir_lrc, delta
  REAL               :: prob_move, prob_create
  REAL, DIMENSION(3) :: ri   ! position of atom i
  REAL, DIMENSION(3) :: zeta ! random numbers

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, activity, prob_move, r_cut, dr_max

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_zvt_lj'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-zVT ensemble'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  activity    = 1.0
  prob_move   = 0.34
  r_cut       = 2.5
  dr_max      = 0.15
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_zvt_lj'
  END IF
  prob_create = (1.0-prob_move)/2.0
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',              nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',     nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',                   temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Activity',                      activity
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Probability of move',           prob_move
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Probability of create/destroy', prob_create
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance',     r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum displacement',          dr_max

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! first call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box

  CALL allocate_arrays ( box, r_cut ) ! Allocate r

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! second call is to get r

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL resize ! Increase the size of the r array

  eng_old = energy ( box, r_cut )
  IF ( eng_old%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_zvt_lj'
  END IF
  pot = eng_old%pot
  vir = eng_old%vir
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: &
       &            'Move ratio', 'Create ratio', 'Destroy ratio', &
       &            'Potential', 'Virial Pressure', 'Density' ] )

  ntry = n ! each step consists of ntry tries (during which n might vary)

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        m_tries = 0
        m_moves = 0
        c_tries = 0
        c_moves = 0
        d_tries = 0
        d_moves = 0

        DO try = 1, ntry ! Begin loop over tries of different kinds

           CALL RANDOM_NUMBER ( zeta(1) ) ! uniform random numbers in RANGE (0,1)

           IF ( zeta(1) < prob_move ) THEN ! try particle move

              m_tries = m_tries + 1

              i       = random_integer ( 1, n )
              ri(:)   = r(:,i)
              eng_old = energy_1 ( ri, i, box, r_cut )

              IF ( eng_old%ovr ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in mc_zvt_lj'
              END IF

              CALL RANDOM_NUMBER ( zeta )           ! Three uniform random numbers in range (0,1)
              zeta    = 2.0*zeta - 1.0              ! now in range (-1,+1)
              ri(:)   = ri(:) + zeta * dr_max / box ! Trial move to new position (in box=1 units)
              ri(:)   = ri(:) - ANINT ( ri(:) )     ! Periodic boundary correction
              eng_new = energy_1 ( ri, i, box, r_cut )

              IF ( .NOT. eng_new%ovr ) THEN ! consider non-overlapping configuration
                 delta = ( eng_new%pot - eng_old%pot ) / temperature
                 IF ( metropolis ( delta ) ) THEN ! accept Metropolis test
                    pot = pot + eng_new%pot - eng_old%pot ! update potential energy
                    vir = vir + eng_new%vir - eng_old%vir ! update virial
                    CALL move ( i, ri )           ! update position
                    m_moves = m_moves + 1         ! increment move counter
                 END IF ! reject Metropolis test
              END IF ! reject overlapping configuration

           ELSE IF ( zeta(1) < prob_move + prob_create ) THEN ! try create
              c_tries = c_tries + 1

              IF ( n+1 > SIZE(r,dim=2) ) CALL resize

              CALL RANDOM_NUMBER ( ri ) ! Three uniform random numbers in range (0,1)
              ri = ri - 0.5             ! now in range (-0.5,+0.5)
              eng_new = energy_1 ( ri, n+1, box, r_cut )

              IF ( .NOT. eng_new%ovr ) THEN ! Consider non-overlapping configuration
                 delta = eng_new%pot / temperature - LOG ( activity / REAL ( n+1 ) )
                 IF ( metropolis ( delta ) ) THEN ! accept Metropolis test
                    CALL create ( ri )            ! create new particle
                    pot     = pot + eng_new%pot   ! update total potential energy
                    vir     = vir + eng_new%vir   ! update total virial
                    c_moves = c_moves + 1         ! increment creation move counter
                 END IF ! reject Metropolis test
              END IF ! reject overlapping configuration

           ELSE ! try destroy
              d_tries = d_tries + 1
              i       = random_integer ( 1, n ) ! Choose particle at random
              eng_old = energy_1 ( r(:,i), i, box, r_cut )

              IF ( eng_old%ovr ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap found on particle removal'
                 STOP 'Error in mc_zvt_lj'
              END IF

              delta = -eng_old%pot/temperature - LOG ( REAL ( n ) / activity )
              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 CALL destroy ( i )            ! Destroy chosen particle
                 pot     = pot - eng_old%pot   ! Update total potential energy
                 vir     = vir - eng_old%vir   ! Update total virial
                 d_moves = d_moves + 1         ! Increment destruction move counter
              END IF ! End accept Metropolis test

           END IF ! End choice of move type

        END DO ! End loop over tries of different kinds

        ! Calculate all variables for this step
        m_ratio     = REAL(m_moves) / REAL(m_tries)
        c_ratio     = REAL(c_moves) / REAL(c_tries)
        d_ratio     = REAL(d_moves) / REAL(d_tries)
        CALL calculate ( )
        CALL blk_add ( [m_ratio,c_ratio,d_ratio,potential,pres_virial,density] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk                   ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r(:,1:n)*box ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  eng_old = energy ( box, r_cut )
  IF ( eng_old%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_zvt_lj'
  END IF
  pot = eng_old%pot
  vir = eng_old%vir
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  CALL calculate ( 'Final check' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r(:,1:n)*box )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    density     = REAL(n) / box**3
    potential   = ( pot + pot_lrc ) / REAL ( n )
    pres_virial = density * temperature + ( vir + vir_lrc ) / box**3

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density',          density
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential energy', potential
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Virial pressure',  pres_virial
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_zvt_lj

