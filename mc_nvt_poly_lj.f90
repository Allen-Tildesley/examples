! mc_nvt_poly_lj.f90
! Monte Carlo, NVT ensemble, polyatomic molecule
PROGRAM mc_nvt_poly_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_mols, write_cnf_mols
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE maths_module,     ONLY : metropolis, random_rotate_quaternion, random_translate_vector
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       potential_1, potential, n, r, e, potential_type

  IMPLICIT NONE

  ! Takes in a configuration of polyatomic molecules (positions and quaternions)
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
  REAL :: dr_max      ! Maximum MC translational displacement
  REAL :: de_max      ! Maximum MC rotational displacement
  REAL :: temperature ! Specified temperature
  REAL :: r_cut       ! Potential cutoff distance

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  ! Composite interaction = pot & vir & ovr variables
  TYPE(potential_type) :: total, partial_old, partial_new

  INTEGER              :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL                 :: delta, m_ratio
  REAL, DIMENSION(3)   :: ri
  REAL, DIMENSION(0:3) :: ei

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, r_cut, dr_max, de_max

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_nvt_poly_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble, polyatomic molecule'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Simulation uses cut-and-shifted potential'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  r_cut       = 2.5
  dr_max      = 0.15
  de_max      = 0.1

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_poly_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum r displacement',    dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum e displacement',    de_max

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of molecules',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density of molecules',  REAL(n) / box**3
  CALL allocate_arrays ( box, r_cut )
  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! Second call is to get r and e
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial energy and overlap check
  total = potential ( box, r_cut )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_poly_lj'
  END IF
  CALL calculate ( 'Initial values' )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( output_unit, variables )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           partial_old = potential_1 ( r(:,i), e(:,i), i, box, r_cut ) ! Old molecule potential, virial, etc

           IF ( partial_old%ovr ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_nvt_poly_lj'
           END IF

           ri = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
           ri = ri - ANINT ( ri )                              ! Periodic boundary correction
           ei = random_rotate_quaternion ( de_max, e(:,i) )    ! Trial rotation

           partial_new = potential_1 ( ri, ei, i, box, r_cut ) ! New molecule potential, virial etc

           IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

              delta = partial_new%pot - partial_old%pot ! Use cut-and-shifted potential
              delta = delta / temperature               ! Divide by temperature

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 total  = total + partial_new - partial_old ! Update potential energy
                 r(:,i) = ri                                ! Update position
                 e(:,i) = ei                                ! Update quaternion
                 moves  = moves + 1                         ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End test for non-overlapping configuration

        END DO ! End loop over atoms

        m_ratio = REAL(moves) / REAL(n)

        ! Calculate and accumulate variables for this step
        CALL calculate ( )
        CALL blk_add ( variables )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )                             ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk              ! Number configuration by block
     CALL write_cnf_mols ( cnf_prefix//sav_tag, n, box, r*box, e ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit ) ! Output run averages

  CALL calculate ( 'Final values' )

  ! Double-check book-keeping for totals, and overlap
  total = potential ( box, r_cut )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_poly_lj'
  END IF
  CALL calculate ( 'Final check' )

  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e ) ! Write out final configuration
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    USE averages_module, ONLY : write_variables
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    ! In this example we simulate using the cut-and-shifted potential only
    ! The values of < p_s >, < e_s > and density should be consistent (for this potential)
    ! There are no long-range or delta corrections
    ! The value of the cut (but not shifted) potential is not used, in this example

    TYPE(variable_type) :: m_r, e_s, p_s
    REAL                :: vol, rho

    ! Preliminary calculations
    vol = box**3        ! Volume
    rho = REAL(n) / vol ! Density

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %msd: indicating if mean squared deviation required
    ! If not set below, %msd adopts its default value of .false.
    ! The %msd and %nam components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move acceptance ratio

    IF ( PRESENT ( string ) ) THEN ! The ratio is meaningless in this case
       m_r = variable_type ( nam = 'Move ratio', val = 0.0 )
    ELSE
       m_r = variable_type ( nam = 'Move ratio', val = m_ratio )
    END IF

    ! Internal energy per molecule (cut-and-shifted potential)
    ! Ideal gas contribution (assuming nonlinear molecules) plus total (cut-and-shifted) PE divided by N
    e_s = variable_type ( nam = 'E/N cut&shift', val = 3.0*temperature + total%pot/REAL(n) )

    ! Pressure (cut-and-shifted potential)
    ! Ideal gas contribution plus total virial divided by V
    p_s = variable_type ( nam = 'P cut&shift', val = rho*temperature + total%vir/vol )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ m_r, e_s, p_s ]

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       CALL write_variables ( output_unit, variables(2:) ) ! Don't write out move ratio
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_nvt_poly_lj

