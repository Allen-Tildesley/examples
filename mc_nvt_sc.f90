! mc_nvt_sc.f90
! Monte Carlo, NVT ensemble, linear hard molecules
PROGRAM mc_nvt_sc

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_mols, write_cnf_mols
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE maths_module,     ONLY : random_rotate_vector, random_translate_vector, orientational_order
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       overlap_1, overlap, n_overlap, n, r, e

  IMPLICIT NONE

  ! Takes in a configuration of linear molecules (positions and orientations)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo for hard particles (the temperature is irrelevant)
  ! Uses no special neighbour lists
  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in reduced units kT=1

  ! Despite the program name, there is nothing here specific to spherocylinders
  ! The model is defined in mc_module

  ! Most important variables
  REAL :: length      ! cylinder length (in units where sigma=1)
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: dr_max      ! maximum MC displacement
  REAL :: de_max      ! maximum MC rotation
  REAL :: eps_box     ! pressure scaling parameter

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL, DIMENSION(3) :: ri, ei
  real               :: m_ratio

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, dr_max, de_max, length, eps_box

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_nvt_sc'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT, hard linear molecules'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock  = 10
  nstep   = 1000
  length  = 3.0
  dr_max  = 0.15
  de_max  = 0.15
  eps_box = 0.005

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_sc'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',           nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',  nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Spherocylinder L/D ratio',   length
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum displacement',       dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Pressure scaling parameter', eps_box

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box (in sigma units)', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density',              REAL(n) / box**3
  CALL allocate_arrays
  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! Second call to get r and e
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial pressure and order calculation and overlap check
  IF ( overlap ( box, length ) ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_sc'
  END IF
  CALL calculate ( 'Initial values' )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( output_unit, variables )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           ri(:) = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )                        ! Periodic boundary correction
           ei(:) = random_rotate_vector ( de_max, e(:,i) )        ! Trial move to new orientation

           IF ( .NOT. overlap_1 ( ri, ei, i, box, length ) ) THEN ! Accept
              r(:,i) = ri(:)     ! Update position
              e(:,i) = ei(:)     ! Update orientation
              moves  = moves + 1 ! Increment move counter
           END IF ! End accept

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

  ! Final overlap check and pressure and order calculation
  IF ( overlap ( box, length ) ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_sc'
  END IF
  CALL calculate ( 'Final values' )

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

    TYPE(variable_type) :: m_r, p, order
    REAL                :: vol, rho, vir, ord

    ! Preliminary calculations (m_ratio, eps_box, box etc are already known)
    vol = box**3                                                           ! Volume
    rho = REAL(n) / vol                                                    ! Density
    vir = REAL ( n_overlap ( box/(1.0+eps_box), length ) ) / (3.0*eps_box) ! Virial
    ord = orientational_order ( e )                                        ! Order

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %msd: indicating if mean squared deviation required
    ! If not set below, %msd adopts its default value of .false.
    ! The %msd and %nam components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move acceptance ratio

    IF ( PRESENT ( string ) ) THEN ! The ratio is meaningless in this case
       m_r = variable_type ( nam = 'Move:ratio', val = 0.0 )
    ELSE
       m_r = variable_type ( nam = 'Move:ratio', val = m_ratio )
    END IF

    ! Pressure in units kT/sigma**3
    ! Ideal gas contribution plus total virial divided by V
    p = variable_type ( nam = 'P', val = rho + vir/vol )

    ! Orientational order parameter
    order = variable_type ( nam = 'Order:orientational', val = ord )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ m_r, p, order ]

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       CALL write_variables ( output_unit, variables(2:) ) ! Don't write out move ratio
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_nvt_sc
