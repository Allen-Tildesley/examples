! mc_nvt_sc.f90 (uses mc_sc_module.f90, utility_module.f90)
! Monte Carlo, NVT ensemble, hard spherocylinders
PROGRAM mc_nvt_sc
  USE utility_module, ONLY : read_cnf_molecules, write_cnf_molecules, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add, &
       &                     random_rotate, orientational_order
  USE mc_sc_module,   ONLY : initialize, finalize, overlap_1, overlap, n_overlap, n, r, e, ne
  IMPLICIT NONE

  ! Takes in a configuration of linear molecules (positions and orientations)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo (the temperature is irrelevant)
  ! Uses no special neighbour lists

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

  INTEGER            :: blk, stp, i, nstep, nblock, moves
  REAL, DIMENSION(3) :: ri, ei ! position and orientation of atom i
  REAL, DIMENSION(3) :: zeta   ! random numbers

  CHARACTER(len=13), PARAMETER :: cnf_prefix = 'mc_nvt_sc.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /run_parameters/ nblock, nstep, dr_max, de_max, length, epsilon

  WRITE(*,'(''mc_nvt_sc'')')
  WRITE(*,'(''Monte Carlo, constant-NVT, hard spherocylinders'')')
  WRITE(*,'(''Results in units sigma = 1'')')

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible defaults for testing
  nblock  = 10
  nstep   = 1000
  length  = 3.0
  dr_max  = 0.15
  de_max  = 0.15
  epsilon = 0.005
  READ(*,nml=run_parameters)
  WRITE(*,'(''Number of blocks'',           t40,i15)'  ) nblock
  WRITE(*,'(''Number of steps per block'',  t40,i15)'  ) nstep
  WRITE(*,'(''Spherocylinder L/D ratio'',   t40,f15.5)') length
  WRITE(*,'(''Maximum displacement'',       t40,f15.5)') dr_max
  WRITE(*,'(''Pressure scaling parameter'', t40,f15.5)') epsilon

  CALL read_cnf_molecules ( cnf_prefix//inp_tag, n, box )
  WRITE(*,'(''Number of particles'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density

  CALL initialize

  CALL read_cnf_molecules ( cnf_prefix//inp_tag, n, box, r, e )

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
  sigma  = sigma / box
  length = length / box
  dr_max = dr_max / box
  WRITE(*,'(''sigma (in box units)'', t40,f15.5)') sigma
  WRITE(*,'(''length (in box units)'',t40,f15.5)') length
  WRITE(*,'(''dr_max (in box units)'',t40,f15.5)') dr_max

  IF ( overlap ( sigma, length ) ) STOP 'Overlap in initial configuration'

  CALL run_begin ( ['Move ratio','Pressure  ','P2 Order  '] ) ! must all be character*10 constants

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

           ri(:) = r(:,i) + zeta * dr_max           ! trial move to new position
           ri(:) = ri(:) - ANINT ( ri(:) )          ! periodic boundary correction
           ei(:) = random_rotate ( de_max, e(:,i) ) ! trial move to new orientation

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

     CALL blk_end ( blk )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk                   ! number configuration by block
     CALL write_cnf_molecules ( cnf_prefix//sav_tag, n, box, r*box, e ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end

  IF ( overlap ( sigma, length ) ) STOP 'Overlap in final configuration'

  CALL write_cnf_molecules ( cnf_prefix//out_tag, n, box, r*box, e )

  CALL finalize

END PROGRAM mc_nvt_sc
