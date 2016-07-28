! mc_zvt_lj.f90
! Monte Carlo, zVT (grand) ensemble, Lennard-Jones atoms
PROGRAM mc_zvt_lj
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE utility_module, ONLY : metropolis, read_cnf_atoms, write_cnf_atoms, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add, random_integer
  USE mc_lj_module,   ONLY : initialize, finalize, resize, energy_1, energy, energy_lrc, &
       &                     move, create, destroy, &
       &                     n, r, ne
  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

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
  INTEGER                  :: try, ntry, ioerr
  INTEGER, DIMENSION(-1:1) :: tries, moves ! count tries and moves of different kinds
  REAL                     :: pot_old, pot_new, pot_lrc, del_pot, vir_old, vir_new, vir_lrc, del_vir, delta
  REAL                     :: prob_move, prob_create
  REAL, DIMENSION(3)       :: ri   ! position of atom i
  REAL, DIMENSION(3)       :: zeta ! random numbers

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, activity, prob_move, r_cut, dr_max

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_zvt_lj'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-zVT, Lennard-Jones'
  WRITE( unit=output_unit, fmt='(a)' ) 'Results in units epsilon = sigma = 1'

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
     STOP 'Error in mc_zvt_lj'
  END IF

  CALL initialize ( r_cut ) ! Allocate r

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! second call is to get r

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL resize ! Increase the size of the r array

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_zvt_lj'
  END IF
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial potential energy (sigma units)', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial pressure (sigma units)',         pressure

  CALL run_begin ( [ CHARACTER(len=15) :: &
       &            'Move ratio', 'Create ratio', 'Destroy ratio', &
       &            'Potential', 'Pressure', 'Density' ] )

  ntry = n ! each step consists of ntry tries (during which n might vary)

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        tries = 0
        moves = 0

        DO try = 1, ntry ! Begin loop over tries of different kinds

           CALL RANDOM_NUMBER ( zeta(1) ) ! uniform random numbers in RANGE (0,1)

           IF ( zeta(1) < prob_move ) THEN ! try particle move

              tries(0) = tries(0) + 1
              CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
              zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

              i = random_integer ( 1, n )
              ri(:) = r(:,i)
              CALL  energy_1 ( ri, i, ne, sigma, r_cut, overlap, pot_old, vir_old )
              IF ( overlap ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in mc_zvt_lj'
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
                    moves(0)  = moves(0) + 1         ! increment move counter
                 END IF ! reject Metropolis test
              END IF ! reject overlapping configuration

           ELSE IF ( zeta(1) < prob_move + prob_create ) THEN ! try create
              tries(1) = tries(1) + 1

              IF ( n+1 > SIZE(r,dim=2) ) CALL resize

              CALL RANDOM_NUMBER ( ri ) ! three uniform random numbers in range (0,1)
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
                    call create ( ri )            ! create new particle
                    pot      = pot + del_pot      ! update total potential energy
                    vir      = vir + del_vir      ! update total virial
                    moves(1) = moves(1) + 1       ! increment creation move counter
                 END IF ! reject Metropolis test
              END IF ! reject overlapping configuration

           ELSE ! try destroy
              tries(-1) = tries(-1) + 1
              i = random_integer ( 1, n )

              CALL energy_1 ( r(:,i), i, ne, sigma, r_cut, overlap, del_pot, del_vir )
              CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc ) ! LRC for n atoms
              del_pot = del_pot + pot_lrc
              del_vir = del_vir + vir_lrc
              CALL energy_lrc ( n-1, sigma, r_cut, pot_lrc, vir_lrc ) ! LRC for n-1 atoms
              del_pot = del_pot - pot_lrc
              del_vir = del_vir - vir_lrc
              del_pot = -del_pot ! change sign for a removal
              del_vir = -del_vir ! change sign for a removal

              IF ( overlap ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap found on particle removal'
                 STOP 'Error in mc_zvt_lj'
              END IF
              delta = del_pot/temperature - LOG ( REAL ( n ) / activity )
              IF ( metropolis ( delta ) ) THEN ! accept Metropolis test
                 call destroy ( i )            ! destroy chosen particle
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

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  potential = pot / REAL ( n )
  density   = REAL(n)*sigma**3
  pressure  = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final potential energy (sigma units)', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final reduced density',                density
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final pressure (sigma units)',         pressure

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_zvt_lj'
  END IF
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a)'           ) 'Final check'
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final potential energy (sigma units)', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final reduced density',                density
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final pressure (sigma units)',         pressure

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )

  CALL finalize

END PROGRAM mc_zvt_lj

