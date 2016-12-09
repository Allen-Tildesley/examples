! mc_npt_hs.f90
! Monte Carlo simulation, NPT ensemble, hard spheres
PROGRAM mc_npt_hs

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE maths_module,     ONLY : metropolis, random_translate_vector
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       overlap_1, overlap, n_overlap, n, r

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at given NPT (the temperature is irrelevant)
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! We take kT=1 throughout defining the unit of energy
  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model,
  ! in this case, for hard spheres, sigma = 1

  ! The logarithm of the box length is sampled uniformly

  ! Most important variables
  REAL :: box      ! Box length
  REAL :: dr_max   ! Maximum MC displacement
  REAL :: db_max   ! Maximum MC box displacement
  REAL :: pressure ! Specified pressure

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL, DIMENSION(3) :: ri
  REAL               :: box_scale, box_new, den_scale, delta, zeta, m_ratio, v_ratio

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, dr_max, db_max, pressure

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_npt_hs'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NPT'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock   = 10
  nstep    = 10000
  dr_max   = 0.15
  db_max   = 0.005
  pressure = 4.0

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_npt_hs'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Pressure',                  pressure
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum displacement',      dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum box displacement',  db_max

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  CALL allocate_arrays ! Allocates r
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! Second call to get r
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial calculation and overlap check
  IF ( overlap ( box ) ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_npt_hs'
  END IF
  CALL calculate ( 'Initial values' )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( variables )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           ri(:) = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )                        ! Periodic boundary correction

           IF ( .NOT. overlap_1 ( ri, i, box ) ) THEN ! Accept
              r(:,i) = ri(:)     ! Update position
              moves  = moves + 1 ! Increment move counter
           END IF ! End accept

        END DO ! End loop over atoms

        m_ratio = REAL(moves) / REAL(n)

        v_ratio = 0.0                   ! Zero volume move counter
        CALL RANDOM_NUMBER ( zeta )     ! Uniform random number in range (0,1)
        zeta      = 2.0*zeta - 1.0      ! Now in range (-1,+1)
        box_scale = EXP ( zeta*db_max ) ! Sampling log(box) and log(vol) uniformly
        box_new   = box * box_scale     ! New box (in sigma units)
        den_scale = 1.0 / box_scale**3  ! Density scaling factor

        IF ( .NOT. overlap ( box_new ) ) THEN ! Test for non-overlapping configuration

           delta = pressure * ( box_new**3 - box**3 ) ! PV term (temperature = 1.0 )
           delta = delta + REAL(n+1) * LOG(den_scale) ! Factor (n+1) consistent with log(box) sampling

           IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
              box     = box_new  ! Update box
              v_ratio = 1.0      ! Set move counter
           END IF ! reject Metropolis test

        END IF ! End test for overlapping configuration

        ! Calculate and accumulate variables for this step
        CALL calculate ( )
        CALL blk_add ( variables )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                        ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ! Output run averages

  ! Final overlap check and pressure calculation
  IF ( overlap ( box ) ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_npt_hs'
  END IF
  CALL calculate ( 'Final values' )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE calculate ( string )
    USE averages_module, ONLY : write_variables
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: m_r, v_r, density
    REAL                :: vol, rho

    ! Preliminary calculations (m_ratio, v_ratio, box etc already known)
    vol = box**3        ! Volume
    rho = REAL(n) / vol ! Density

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move acceptance ratios

    IF ( PRESENT ( string ) ) THEN ! The ratios are meaningless in this case
       m_r = variable_type ( nam = 'Move ratio',   val = 0.0 )
       v_r = variable_type ( nam = 'Volume ratio', val = 0.0 )
    ELSE
       m_r = variable_type ( nam = 'Move ratio',   val = m_ratio )
       v_r = variable_type ( nam = 'Volume ratio', val = v_ratio )
    END IF

    ! Density
    density = variable_type ( nam = 'Density', val = rho )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ m_r, v_r, density ]

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       CALL write_variables ( variables(3:) ) ! Don't write out ratios
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_npt_hs
