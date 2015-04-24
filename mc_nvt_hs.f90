! mc_nvt_hs.f90
! Monte Carlo simulation, NVT ensemble, hard spheres
PROGRAM mc_nvt_hs
  USE utility_module, ONLY : read_cnf_atoms, write_cnf_atoms, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_hs_module,   ONLY : initialize, finalize, overlap_1, overlap, n_overlap, n, r, ne
  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo (the temperature is irrelevant)
  ! Uses no special neighbour lists

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in reduced units sigma = 1 kT=1

  ! Most important variables
  REAL :: sigma       ! atomic diameter (in units where box=1)
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: pressure    ! measured pressure in units kT/sigma**3
  REAL :: dr_max      ! maximum MC displacement
  REAL :: epsilon     ! pressure scaling parameter
  REAL :: move_ratio  ! acceptance ratio of moves (to be averaged)

  INTEGER            :: blk, stp, i, nstep, nblock, moves
  REAL, DIMENSION(3) :: ri   ! position of atom i
  REAL, DIMENSION(3) :: zeta ! random numbers

  CHARACTER(len=13), PARAMETER :: cnf_prefix = 'mc_nvt_hs.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /run_parameters/ nblock, nstep, dr_max, epsilon

  WRITE(*,'(''mc_nvt_hs'')')
  WRITE(*,'(''Monte Carlo, constant-NVT, hard spheres'')')
  WRITE(*,'(''Results in units sigma = 1'')')

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible defaults for testing
  nblock  = 10
  nstep   = 1000
  dr_max  = 0.15
  epsilon = 0.005
  READ(*,nml=run_parameters)
  WRITE(*,'(''Number of blocks'',           t40,i15)'  ) nblock
  WRITE(*,'(''Number of steps per block'',  t40,i15)'  ) nstep
  WRITE(*,'(''Maximum displacement'',       t40,f15.5)') dr_max
  WRITE(*,'(''Pressure scaling parameter'', t40,f15.5)') epsilon

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box )
  WRITE(*,'(''Number of particles'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density

  CALL initialize ! Allocates r
  
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r )

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
  sigma  = sigma / box
  dr_max = dr_max / box
  WRITE(*,'(''sigma (in box units)'', t40,f15.5)') sigma
  WRITE(*,'(''dr_max (in box units)'',t40,f15.5)') dr_max

  IF ( overlap ( sigma ) ) STOP 'Overlap in initial configuration'

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

           ri(:) = r(:,i) + zeta * dr_max  ! trial move to new position
           ri(:) = ri(:) - ANINT ( ri(:) ) ! periodic boundary correction

           IF ( .NOT. overlap_1 ( ri, i, ne, sigma ) ) THEN ! accept
              r(:,i) = ri(:)     ! update position
              moves  = moves + 1 ! increment move counter
           END IF ! reject overlapping configuration

        END DO ! End loop over atoms

        ! Calculate all variables for this step
        move_ratio = REAL(moves) / REAL(n)
        pressure = REAL ( n_overlap ( (1.0+epsilon)*sigma ) ) / (3.0*epsilon) ! virial part
        pressure = density + pressure * sigma**3 ! convert to sigma units and add ideal gas part
        CALL blk_add ( [move_ratio,pressure] )

     END DO ! End loop over steps

     CALL blk_end ( blk )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end

  IF ( overlap ( sigma ) ) STOP 'Overlap in final configuration'

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )

  CALL finalize ! Deallocates r

END PROGRAM mc_nvt_hs
