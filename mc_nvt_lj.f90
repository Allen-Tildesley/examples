! mc_nvt_lj.f90 (also uses mc_nvt_lj_module.f90 and utility_module.f90)
! Monte Carlo simulation, constant-NVT ensemble, Lennard-Jones atoms
PROGRAM mc_nvt_lj
  USE utility_module,   ONLY : read_cnf_atoms, write_cnf_atoms, run_begin, run_end, blk_begin, blk_end, stp_end
  USE mc_nvt_lj_module, ONLY : energy_i, energy, pot_lrc, vir_lrc
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
  REAL               :: pot_old, pot_new, vir_old, vir_new, exponent
  REAL, DIMENSION(3) :: ri   ! position of atom i
  REAL, DIMENSION(4) :: zeta ! random numbers

  REAL,              PARAMETER :: exponent_guard = 75.0
  CHARACTER(len=13), PARAMETER :: cnf_prefix = 'md_nvt_lj.cnf'
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
  WRITE(*,'(''Potential LRC (sigma units)'',t40,f15.5)') pot_lrc ( sigma, r_cut, density )
  WRITE(*,'(''Virial LRC (sigma units)'',   t40,f15.5)') vir_lrc ( sigma, r_cut, density )

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

  CALL energy ( sigma, r_cut, pot, vir, overlap )
  IF ( overlap ) STOP 'Overlap in initial configuration'
  potential = pot / REAL ( n ) + pot_lrc ( sigma, r_cut, density )
  pressure  = density * temperature + vir / box**3 + vir_lrc ( sigma, r_cut, density ) * density
  WRITE(*,'(''Initial potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Initial pressure (sigma units)'',        t40,f15.5)') pressure

  CALL run_begin ( ['Move ratio','Pressure  ','Potential '] ) ! must all be character*10 constants

  DO blk = 1, nblock

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           CALL RANDOM_NUMBER ( zeta )     ! four uniform random numbers in range (0,1)
           zeta(1:3) = 2.0*zeta(1:3) - 1.0 ! first three now in range (-1,+1)

           ri(:) = r(:,i)
           CALL  energy_i ( ri, i, j_ne_i, sigma, r_cut, pot_old, vir_old, overlap )
           IF ( overlap ) STOP 'Overlap in current configuration'
           ri(:) = ri(:) + zeta(1:3) * dr_max ! trial move to new position
           ri(:) = ri(:) - ANINT ( ri(:) )    ! periodic boundary correction
           CALL  energy_i ( ri, i, j_ne_i, sigma, r_cut, pot_new, vir_new, overlap )

           IF ( .NOT. overlap ) THEN ! consider non-overlapping configuration

              exponent = ( pot_new - pot_old ) / temperature ! goes into Boltzmann factor

              IF ( exponent < exponent_guard ) THEN ! consider not-too-high energy change

                 IF ( exponent <= 0.0 ) THEN ! accept because downhill
                    pot    = pot + pot_new - pot_old
                    vir    = vir + vir_new - vir_old
                    r(:,i) = ri(:)
                    moves  = moves + 1

                 ELSEIF ( EXP(-exponent) > zeta(4) ) THEN ! accept because Metropolis test
                    pot    = pot + pot_new - pot_old
                    vir    = vir + vir_new - vir_old
                    r(:,i) = ri(:)
                    moves  = moves + 1

                 END IF ! reject because Metropolis test

              END IF ! reject too-high energy change

           END IF ! reject overlapping configuration

        END DO ! End loop over atoms

        ! Calculate all variables for this step
        move_ratio = REAL(moves) / REAL(n)
        potential  = pot / REAL(n) + pot_lrc ( sigma, r_cut, density )
        pressure   = density * temperature + vir / box**3 + vir_lrc ( sigma, r_cut, density ) * density
        CALL stp_end ( [move_ratio,potential,pressure] )

     END DO ! End loop over steps

     CALL blk_end ( blk )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end

  potential  = pot / REAL ( n ) + pot_lrc ( sigma, r_cut, density )
  pressure   = density * temperature + vir / box**3 + vir_lrc ( sigma, r_cut, density ) * density
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL energy ( sigma, r_cut, pot, vir, overlap )
  IF ( overlap ) STOP 'Overlap in final configuration'
  potential  = pot / REAL ( n ) + pot_lrc ( sigma, r_cut, density )
  pressure   = density * temperature + vir / box**3 + vir_lrc ( sigma, r_cut, density ) * density
  WRITE(*,'(''Final check'')')
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )

  DEALLOCATE ( r )

END PROGRAM mc_nvt_lj

