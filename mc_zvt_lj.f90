! mc_zvt_lj.f90
! Monte Carlo, zVT (grand) ensemble
PROGRAM mc_zvt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis, random_integer, random_translate_vector
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       potential_1, potential, &
       &                       move, create, destroy, n, r, &
       &                       potential_type, OPERATOR(+), OPERATOR(-)

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
  REAL :: box         ! Box length
  REAL :: dr_max      ! Maximum MC displacement
  REAL :: temperature ! Specified temperature
  REAL :: activity    ! Specified activity z
  REAL :: r_cut       ! Potential cutoff distance

  ! Quantities to be averaged
  REAL :: m_ratio ! Acceptance ratio of moves
  REAL :: c_ratio ! Acceptance ratio of creation attempts
  REAL :: d_ratio ! Acceptance ratio of destruction attempts
  REAL :: density ! Density
  REAL :: en_c    ! Internal energy per atom for simulated, cut, potential
  REAL :: en      ! Internal energy per atom for full potential with LRC
  REAL :: p_c     ! Pressure for simulated, cut, potential
  REAL :: p       ! Pressure for full potential with LRC

  ! Composite interaction = pot & vir & overlap variables
  TYPE(potential_type) :: total, partial_old, partial_new

  INTEGER            :: blk, stp, i, nstep, nblock
  INTEGER            :: try, ntry, ioerr
  INTEGER            :: m_tries, m_moves ! count tries and moves
  INTEGER            :: c_tries, c_moves ! count tries and moves for creation
  INTEGER            :: d_tries, d_moves ! count tries and moves for destruction
  REAL               :: prob_move, prob_create, delta, zeta
  REAL, DIMENSION(3) :: ri

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

  prob_create = (1.0-prob_move)/2.0 ! So that create and destroy have equal probabilities

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

  n = n * 2 ! Increase n for array allocation
  CALL allocate_arrays ( box, r_cut ) ! Allocate r
  n = n / 2 ! Restore value of n

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! second call is to get r

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial energy and overlap check
  total = potential ( box, r_cut ) 
  IF ( total%overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_zvt_lj'
  END IF
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Create ratio', 'Destroy ratio', &
       &            'Density', 'E/N (cut)', 'P (cut)', 'E/N (full)', 'P (full)' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        m_tries = 0
        m_moves = 0
        c_tries = 0
        c_moves = 0
        d_tries = 0
        d_moves = 0

        ntry = n ! Each step consists of ntry tries (during which n might vary)

        DO try = 1, ntry ! Begin loop over tries of different kinds

           CALL RANDOM_NUMBER ( zeta ) ! uniform random number in range (0,1)

           IF ( zeta < prob_move ) THEN ! Try particle move

              m_tries = m_tries + 1

              i = random_integer ( 1, n ) ! Choose moving particle at random

              partial_old = potential_1 ( r(:,i), i, box, r_cut ) ! Old atom potential, virial etc

              IF ( partial_old%overlap ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in mc_zvt_lj'
              END IF

              ri(:) = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
              ri(:) = ri(:) - ANINT ( ri(:) )                        ! Periodic boundary correction

              partial_new = potential_1 ( ri, i, box, r_cut ) ! New atom potential, virial etc

              IF ( .NOT. partial_new%overlap ) THEN ! Test for non-overlapping configuration

                 delta = partial_new%pot_c - partial_old%pot_c ! Use cut (but not shifted) potential
                 delta = delta / temperature                   ! Divide by temperature

                 IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                    total = total + partial_new - partial_old ! Update total values
                    CALL move ( i, ri )                       ! Update position
                    m_moves = m_moves + 1                     ! Increment move counter
                 END IF ! End accept Metropolis test

              END IF ! End test for overlapping configuration

           ELSE IF ( zeta < prob_move + prob_create ) THEN ! Try create

              c_tries = c_tries + 1

              IF ( n+1 > SIZE(r,dim=2) ) then
                 WRITE ( unit=error_unit, fmt='(a,2i5)') 'n has grown too large', n+1, SIZE(r,dim=2)
                 STOP 'Error in mc_zvt_lj'
              END IF

              CALL RANDOM_NUMBER ( ri ) ! Three uniform random numbers in range (0,1)
              ri = ri - 0.5             ! now in range (-0.5,+0.5) for box=1 units

              partial_new = potential_1 ( ri, n+1, box, r_cut ) ! New atom potential, virial, etc

              IF ( .NOT. partial_new%overlap ) THEN ! Test for non-overlapping configuration

                 delta = partial_new%pot_c / temperature                  ! Use cut (not shifted) potential
                 delta = delta - LOG ( activity * box**3 / REAL ( n+1 ) ) ! Activity term for creation

                 IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                    CALL create ( ri )            ! Create new particle
                    total   = total + partial_new ! Update total values
                    c_moves = c_moves + 1         ! Increment creation move counter
                 END IF ! End accept Metropolis test

              END IF ! End test for overlapping configuration

           ELSE ! try destroy

              d_tries = d_tries + 1

              i = random_integer ( 1, n ) ! Choose particle at random

              partial_old = potential_1 ( r(:,i), i, box, r_cut ) ! Old atom potential, virial, etc

              IF ( partial_old%overlap ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap found on particle removal'
                 STOP 'Error in mc_zvt_lj'
              END IF

              delta = -partial_old%pot_c / temperature                    ! Use cut (not shifted) potential
              delta = delta - LOG ( REAL ( n ) / ( activity * box**3 )  ) ! Activity term for destruction

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 CALL destroy ( i )            ! Destroy chosen particle
                 total   = total - partial_old ! Update total values
                 d_moves = d_moves + 1         ! Increment destruction move counter
              END IF ! End accept Metropolis test

           END IF ! End choice of move type

        END DO ! End loop over tries of different kinds

        IF ( m_tries > 0 ) THEN
           m_ratio = REAL(m_moves) / REAL(m_tries)
        ELSE
           m_ratio = 0.0
        END IF

        IF ( c_tries > 0 ) THEN
           c_ratio = REAL(c_moves) / REAL(c_tries)
        ELSE
           c_ratio = 0.0
        END IF

        IF ( d_tries > 0 ) THEN
           d_ratio = REAL(d_moves) / REAL(d_tries)
        ELSE
           d_ratio = 0.0
        END IF

        ! Calculate all variables for this step
        CALL calculate ( )
        CALL blk_add ( [m_ratio,c_ratio,d_ratio,density,en_c,p_c,en,p] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk                   ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r(:,1:n)*box ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  ! Double-check book-keeping of totals and overlap
  total = potential ( box, r_cut )
  IF ( total%overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_zvt_lj'
  END IF
  CALL calculate ( 'Final check' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r(:,1:n)*box )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    USE mc_module, ONLY : potential_lrc, pressure_lrc, pressure_delta
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates the properties of interest from total values
    ! and optionally writes them out (e.g. at the start and end of the run)
    ! In this example we simulate using the cut (but not shifted) potential
    ! The values of < p_c >,  < en_c > and < density > should be consistent (for this potential)
    ! For comparison, long-range corrections are also applied to give
    ! estimates of < en > and < p > for the full (uncut) potential
    ! The value of the cut-and-shifted potential pot_s is not used, in this example

    density = REAL(n) / box**3                        ! Number density N/V
    en_c    = total%pot_c / REAL ( n )                ! PE/N for cut (but not shifted) potential
    en_c    = en_c + 1.5 * temperature                ! Add ideal gas contribution KE/N to give E_c/N
    en      = en_c + potential_lrc ( density, r_cut ) ! Add long-range contribution to give E/N estimate
    p_c     = total%vir / box**3                      ! Virial contribution to P_c
    p_c     = p_c + density * temperature             ! Add ideal gas contribution to P_c
    p       = p_c + pressure_lrc ( density, r_cut )   ! Add long-range contribution to give P
    p_c     = p_c + pressure_delta ( density, r_cut ) ! Add delta correction to P_c (not needed for P)

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density',    density
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (cut)',  en_c
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (cut)',    p_c
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (full)', en
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (full)',   p
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_zvt_lj

