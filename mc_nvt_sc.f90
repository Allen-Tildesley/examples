! mc_nvt_sc.f90
! Monte Carlo, NVT ensemble, hard spherocylinders
PROGRAM mc_nvt_sc
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE utility_module, ONLY : read_cnf_mols, write_cnf_mols, time_stamp, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add, &
       &                     random_rotate_vector, orientational_order
  USE mc_sc_module,   ONLY : allocate_arrays, deallocate_arrays, overlap_1, overlap, n_overlap, n, r, e, ne
  IMPLICIT NONE

  ! Takes in a configuration of linear molecules (positions and orientations)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo (the temperature is irrelevant)
  ! Uses no special neighbour lists
  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in reduced units sigma = 1 kT=1

  ! Most important variables
  REAL :: sigma       ! cylinder diameter (in units where box=1)
  REAL :: length      ! cylinder length (in units where box=1)
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: pressure    ! measured pressure in units kT/sigma**3
  REAL :: order       ! orientational order parameter
  REAL :: dr_max      ! maximum MC displacement
  REAL :: de_max      ! maximum MC rotation
  REAL :: epsilon     ! pressure scaling parameter
  REAL :: move_ratio  ! acceptance ratio of moves (to be averaged)

  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL, DIMENSION(3) :: ri, ei ! position and orientation of atom i
  REAL, DIMENSION(3) :: zeta   ! random numbers

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, dr_max, de_max, length, epsilon

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_nvt_sc'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT, hard spherocylinders'
  WRITE( unit=output_unit, fmt='(a)' ) 'Results in units sigma = 1'
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
  sigma = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Reduced density', density

  CALL allocate_arrays

  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! second call to get r and e

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
  sigma  = sigma / box
  length = length / box
  dr_max = dr_max / box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'sigma (in box units)',  sigma
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'length (in box units)', length
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'dr_max (in box units)', dr_max

  IF ( overlap ( sigma, length ) ) then
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_sc'
  END IF

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Pressure', 'P2 Order' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

           ri(:) = r(:,i) + zeta * dr_max                  ! trial move to new position
           ri(:) = ri(:) - ANINT ( ri(:) )                 ! periodic boundary correction
           ei(:) = random_rotate_vector ( de_max, e(:,i) ) ! trial move to new orientation

           IF ( .NOT. overlap_1 ( ri, ei, i, ne, sigma, length ) ) THEN ! accept
              r(:,i) = ri(:)     ! update position
              e(:,i) = ei(:)     ! update orientation
              moves  = moves + 1 ! increment move counter
           END IF ! reject overlapping configuration

        END DO ! End loop over atoms

        ! Calculate all variables for this step
        move_ratio = REAL(moves) / REAL(n)
        pressure = REAL ( n_overlap ( (1.0+epsilon)*sigma, (1.0+epsilon)*length ) ) / (3.0*epsilon) ! virial part
        pressure = density + pressure * sigma**3 ! convert to sigma units and add ideal gas part
        order    = orientational_order ( e )
        CALL blk_add ( [move_ratio,pressure,order] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk              ! number configuration by block
     CALL write_cnf_mols ( cnf_prefix//sav_tag, n, box, r*box, e ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  IF ( overlap ( sigma, length ) ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_sc'
  END IF
  
  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays

END PROGRAM mc_nvt_sc
