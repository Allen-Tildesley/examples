! mc_nvt_sc.f90
! Monte Carlo, NVT ensemble, linear hard molecules
PROGRAM mc_nvt_sc

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_mols, write_cnf_mols
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
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
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: dr_max      ! maximum MC displacement
  REAL :: de_max      ! maximum MC rotation
  REAL :: epsilon     ! pressure scaling parameter

  ! Quantities to be averaged
  REAL :: move_ratio ! acceptance ratio of moves
  REAL :: pkt        ! pressure in units kT/sigma**3
  REAL :: order      ! orientational order parameter

  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL, DIMENSION(3) :: ri, ei

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, dr_max, de_max, length, epsilon

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_nvt_sc'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT, hard linear molecules'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! initialize random number generator

  ! Set sensible defaults for testing
  nblock  = 10
  nstep   = 1000
  length  = 3.0
  dr_max  = 0.15
  de_max  = 0.15
  epsilon = 0.005
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_sc'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',           nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',  nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Spherocylinder L/D ratio',   length
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum displacement',       dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Pressure scaling parameter', epsilon

  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box ) ! first call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box (in sigma units)', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Reduced density', density

  CALL allocate_arrays

  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! second call to get r and e

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  IF ( overlap ( box, length ) ) then
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_sc'
  END IF

  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'P/kT', 'Order' ] )

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

        move_ratio = REAL(moves) / REAL(n)

        ! Calculate all variables for this step
        call calculate ( )
        CALL blk_add ( [move_ratio,pkt,order] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk              ! number configuration by block
     CALL write_cnf_mols ( cnf_prefix//sav_tag, n, box, r*box, e ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  IF ( overlap ( box, length ) ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_sc'
  END IF

  CALL calculate ( 'Final values' )

  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates the properties of interest
    ! and optionally writes them out (e.g. at the start and end of the run)

    REAL :: box_scaled ! scaled box for pressure calculation

    box_scaled = box / (1.0+epsilon) 
    pkt        = REAL ( n_overlap ( box_scaled, length ) ) / (3.0*epsilon*box**3) ! virial part of P/kT
    pkt        = pkt + density                                                    ! add ideal gas part of P/kT
    order      = orientational_order ( e )

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P/kT',  pkt
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Order', order
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_nvt_sc
