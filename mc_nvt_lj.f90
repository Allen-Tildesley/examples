! mc_nvt_lj.f90 (also uses mc_lj_module.f90 and utility_module.f90)
! Monte Carlo simulation, constant-NVT ensemble, Lennard-Jones atoms
PROGRAM mc_nvt_lj
  USE utility_module, ONLY : metropolis, read_cnf_atoms, write_cnf_atoms, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_lj_module,   ONLY : energy_1, energy, energy_lrc, n, r, ne
  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1

  ! Most important variables
  REAL :: sigma       ! atomic diameter (in units where box=1)
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: dr_max      ! maximum MC displacement
  REAL :: temperature ! specified temperature
  REAL :: r_cut       ! potential cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: vir         ! total virial
  REAL :: move_ratio  ! acceptance ratio of moves (to be averaged)
  REAL :: pressure    ! pressure (LJ sigma=1 units, to be averaged)
  REAL :: potential   ! potential energy per atom (LJ sigma=1 units, to be averaged)

  LOGICAL            :: overlap
  INTEGER            :: blk, stp, i, nstep, nblock, moves
  REAL               :: pot_old, pot_new, pot_lrc, vir_old, vir_new, vir_lrc, delta
  REAL, DIMENSION(3) :: ri   ! position of atom i
  REAL, DIMENSION(3) :: zeta ! random numbers

  CHARACTER(len=13), PARAMETER :: cnf_prefix = 'mc_nvt_lj.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /run_parameters/ nblock, nstep, temperature, r_cut, dr_max

  WRITE(*,'(''mc_nvt_lj'')')
  WRITE(*,'(''Monte Carlo, constant-NVT, Lennard-Jones'')')
  WRITE(*,'(''Results in units epsilon = sigma = 1'')')

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  r_cut       = 2.5
  dr_max      = 0.15
  READ(*,nml=run_parameters)
  WRITE(*,'(''Number of blocks'',         t40,i15)'  ) nblock
  WRITE(*,'(''Number of steps per block'',t40,i15)'  ) nstep
  WRITE(*,'(''Temperature'',              t40,f15.5)') temperature
  WRITE(*,'(''Potential cutoff distance'',t40,f15.5)') r_cut
  WRITE(*,'(''Maximum displacement'',     t40,f15.5)') dr_max

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box )
  WRITE(*,'(''Number of particles'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density

  ALLOCATE ( r(3,n) )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r )

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
  sigma  = sigma / box
  r_cut  = r_cut / box
  dr_max = dr_max / box
  WRITE(*,'(''sigma (in box units)'',t40,f15.5)') sigma
  WRITE(*,'(''r_cut (in box units)'',t40,f15.5)') r_cut
  IF ( r_cut > 0.5 ) STOP 'r_cut too large '

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) STOP 'Overlap in initial configuration'
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(*,'(''Initial potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Initial pressure (sigma units)'',        t40,f15.5)') pressure

  CALL run_begin ( ['Move ratio','Potential ','Pressure  '] ) ! must all be character*10 constants

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

           ri(:) = r(:,i)
           CALL  energy_1 ( ri, i, ne, sigma, r_cut, overlap, pot_old, vir_old )
           IF ( overlap ) STOP 'Overlap in current configuration'
           ri(:) = ri(:) + zeta * dr_max   ! trial move to new position
           ri(:) = ri(:) - ANINT ( ri(:) ) ! periodic boundary correction
           CALL  energy_1 ( ri, i, ne, sigma, r_cut, overlap, pot_new, vir_new )

           IF ( .NOT. overlap ) THEN ! consider non-overlapping configuration
              delta = ( pot_new - pot_old ) / temperature
              IF ( metropolis ( delta ) ) THEN    ! accept Metropolis test
                 pot    = pot + pot_new - pot_old ! update potential energy
                 vir    = vir + vir_new - vir_old ! update virial
                 r(:,i) = ri(:)                   ! update position
                 moves  = moves + 1               ! increment move counter
              END IF ! reject Metropolis test
           END IF ! reject overlapping configuration

        END DO ! End loop over atoms

        ! Calculate all variables for this step
        move_ratio = REAL(moves) / REAL(n)
        potential  = pot / REAL(n)
        pressure   = density * temperature + vir / box**3
        CALL blk_add ( [move_ratio,potential,pressure] )

     END DO ! End loop over steps

     CALL blk_end ( blk )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end

  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL energy ( sigma, r_cut, overlap, pot, vir )
  IF ( overlap ) STOP 'Overlap in final configuration'
  CALL energy_lrc ( n, sigma, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(*,'(''Final check'')')
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )

  DEALLOCATE ( r )

END PROGRAM mc_nvt_lj

