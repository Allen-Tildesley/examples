! mc_nvt_lj.f90
! Monte Carlo, NVT ensemble
PROGRAM mc_nvt_lj
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis, random_translate_vector
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       potential_1, potential, move, n, r, &
       &                       potential_type

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

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
  REAL :: en_c    ! Internal energy per atom for simulated, cut, potential
  REAL :: en_f    ! Internal energy per atom for full potential with LRC
  REAL :: p_c     ! Pressure for simulated, cut, potential
  REAL :: p_f     ! Pressure for full potential with LRC
  REAL :: tc      ! Configurational temperature

  ! Composite interaction = pot & vir & ovr variables
  TYPE(potential_type) :: total, partial_old, partial_new
  
  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL               :: delta
  REAL, DIMENSION(3) :: ri

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, r_cut, dr_max

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_nvt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Simulation uses cut (but not shifted) potential'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize! random number generator

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  r_cut       = 2.5
  dr_max      = 0.15
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_lj'
  END IF

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum displacement',      dr_max

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
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_lj'
  END IF
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'E/N (cut)', &
       &            'P (cut)', 'E/N (full)', 'P (full)', 'T (con)' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           partial_old = potential_1 ( r(:,i), i, box, r_cut ) ! Old atom potential, virial etc

           IF ( partial_old%ovr ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_nvt_lj'
           END IF

           ri(:) = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )                        ! Periodic boundary correction

           partial_new = potential_1 ( ri, i, box, r_cut ) ! New atom potential, virial etc

           IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

              delta = partial_new%pot - partial_old%pot ! Use cut (but not shifted) potential
              delta = delta / temperature

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 total = total + partial_new - partial_old ! Update total values
                 CALL move ( i, ri )                       ! Update position
                 moves = moves + 1                         ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End test for overlapping configuration

        END DO ! End loop over atoms

        m_ratio = REAL(moves) / REAL(n)

        CALL calculate ( )
        CALL blk_add ( [m_ratio,en_c,p_c,en_f,p_f,tc] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk        ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  ! Double-check book-keeping for totals, and overlap
  total = potential ( box, r_cut )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_lj'
  END IF
  CALL calculate ( 'Final check' )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    USE mc_module, ONLY : potential_lrc, pressure_lrc, pressure_delta, force_sq
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates the properties of interest from total values
    ! and optionally writes them out (e.g. at the start and end of the run)
    ! In this example we simulate using the cut (but not shifted) potential
    ! The values of < p_c >, < en_c > and density should be consistent (for this potential)
    ! For comparison, long-range corrections are also applied to give
    ! estimates of < en_f > and < p_f > for the full (uncut) potential
    ! The value of the cut-and-shifted potential is not used, in this example

    en_c = total%pot / REAL ( n )                  ! PE/N for cut (but not shifted) potential
    en_c = en_c + 1.5 * temperature                ! Add ideal gas contribution KE/N to give E_c/N
    en_f = en_c + potential_lrc ( density, r_cut ) ! Add long-range contribution to give E_f/N estimate

    p_c = total%vir / box**3                      ! Virial contribution to P_c
    p_c = p_c + density * temperature             ! Add ideal gas contribution to P_c
    p_f = p_c + pressure_lrc ( density, r_cut )   ! Add long-range contribution to give P_f
    p_c = p_c + pressure_delta ( density, r_cut ) ! Add delta correction to P_c (not needed for P_f)

    tc = force_sq ( box, r_cut ) / total%lap ! Configurational temperature

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (cut)',  en_c
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (cut)',    p_c
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (full)', en_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (full)',   p_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (con)',    tc
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_nvt_lj

