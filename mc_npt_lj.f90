! mc_npt_lj.f90 (uses mc_lj_module.f90, utility_module.f90)
! Monte Carlo, NPT ensemble, Lennard-Jones atoms
PROGRAM mc_npt_lj
  USE utility_module, ONLY : metropolis, read_cnf_atoms, write_cnf_atoms, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_lj_module,   ONLY : initialize, finalize, energy_1, energy, energy_lrc, n, r, ne
  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature and pressure
  ! Uses no special neighbour lists

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1

  ! This program uses the known scaling of the separate parts of the LJ potential
  ! with distances (i.e. with box scaling) to handle the volume move.
  ! However, this requires us to scale the cutoff distance with the box
  ! The logarithm of the box length is sampled uniformly

  ! Most important variables
  REAL               :: sigma          ! atomic diameter (in units where box=1)
  REAL               :: box            ! box length (in units where sigma=1)
  REAL               :: dr_max         ! maximum MC particle displacement
  REAL               :: db_max         ! maximum MC box displacement
  REAL               :: temperature    ! specified temperature
  REAL               :: pressure_inp   ! specified pressure
  REAL               :: r_cut          ! potential cutoff distance
  REAL, DIMENSION(2) :: pot            ! total potential energy (LJ12 and LJ6 parts separate)
  REAL, DIMENSION(2) :: vir            ! total virial (LJ12 and LJ6 parts separate)
  REAL               :: move_ratio     ! acceptance ratio of moves (to be averaged)
  REAL               :: box_move_ratio ! acceptance ratio of box moves (to be averaged)
  REAL               :: density        ! reduced density n*sigma**3/box**3 (to be averaged)
  REAL               :: pressure       ! pressure (LJ sigma=1 units, to be averaged)
  REAL               :: potential      ! potential energy per atom (LJ sigma=1 units, to be averaged)

  LOGICAL            :: overlap
  INTEGER            :: blk, stp, i, nstep, nblock, moves
  REAL               :: sigma_scale, box_new, sigma_new, density_new, delta
  REAL, DIMENSION(2) :: pot_old, pot_new, pot_lrc, vir_old, vir_new, vir_lrc, lj_scale
  REAL, DIMENSION(3) :: ri   ! position of atom i
  REAL, DIMENSION(3) :: zeta ! random numbers

  CHARACTER(len=13), PARAMETER :: cnf_prefix = 'mc_npt_lj.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /run_parameters/ nblock, nstep, temperature, pressure_inp, r_cut, dr_max, db_max

  WRITE(*,'(''mc_npt_lj'')')
  WRITE(*,'(''Monte Carlo, constant-NPT, Lennard-Jones'')')
  WRITE(*,'(''Results in units epsilon = sigma = 1'')')

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible defaults for testing
  nblock       = 10
  nstep        = 1000
  temperature  = 0.7
  pressure_inp = 0.1
  r_cut        = 2.5
  dr_max       = 0.15
  db_max       = 0.025
  READ(*,nml=run_parameters)
  WRITE(*,'(''Number of blocks'',                      t40,i15)'  ) nblock
  WRITE(*,'(''Number of steps per block'',             t40,i15)'  ) nstep
  WRITE(*,'(''Temperature'',                           t40,f15.5)') temperature
  WRITE(*,'(''Pressure'',                              t40,f15.5)') pressure_inp
  WRITE(*,'(''Potential cutoff (sigma units)'',        t40,f15.5)') r_cut
  WRITE(*,'(''Maximum displacement (sigma units)'',    t40,f15.5)') dr_max
  WRITE(*,'(''Maximum box displacement (sigma units)'',t40,f15.5)') db_max

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box )
  WRITE(*,'(''Number of particles'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density

  CALL initialize ! Allocate r

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r )

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
  sigma  = 1.0 / box
  r_cut  = r_cut / box
  dr_max = dr_max / box
  WRITE(*,'(''sigma (in box units)'',t40,f15.5)') sigma
  IF ( r_cut > 0.5 ) STOP 'r_cut too large '

  CALL energy ( sigma, r_cut, overlap, pot2=pot, vir2=vir )
  IF ( overlap ) STOP 'Overlap in initial configuration'
  CALL energy_lrc ( n, sigma, r_cut, pot2=pot_lrc, vir2=vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = SUM ( pot ) / REAL ( n )
  pressure  = density * temperature + SUM ( vir ) / box**3
  WRITE(*,'(''Initial potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Initial pressure (sigma units)'',        t40,f15.5)') pressure

  CALL run_begin ( [ CHARACTER(len=15) :: &
       &            'Move ratio', 'Box ratio', 'Density', 'Potential', 'Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

           ri(:) = r(:,i)
           CALL  energy_1 ( ri, i, ne, sigma, r_cut, overlap, pot2=pot_old, vir2=vir_old )
           IF ( overlap ) STOP 'Overlap in current configuration'
           ri(:) = ri(:) + zeta * dr_max   ! trial move to new position
           ri(:) = ri(:) - ANINT ( ri(:) ) ! periodic boundary correction
           CALL  energy_1 ( ri, i, ne, sigma, r_cut, overlap, pot2=pot_new, vir2=vir_new )

           IF ( .NOT. overlap ) THEN ! consider non-overlapping configuration

              delta = SUM ( pot_new - pot_old ) / temperature

              IF (  metropolis ( delta )  ) THEN ! accept Metropolis test
                 pot = pot + pot_new - pot_old   ! update potential energy
                 vir = vir + vir_new - vir_old   ! update virial
                 r(:,i) = ri(:)                  ! update position
                 moves  = moves + 1              ! increment move counter
              END IF ! reject Metropolis test

           END IF ! reject overlapping configuration

        END DO ! End loop over atoms

        move_ratio = REAL(moves) / REAL(n)

        box_move_ratio = 0.0
        CALL RANDOM_NUMBER ( zeta(1) )                 ! uniform random number in range (0,1)
        zeta(1)     = 2.0*zeta(1) - 1.0                ! now in range (-1,+1)
        sigma_scale = EXP ( zeta(1)*db_max )           ! sampling log(sigma), log(box) and log(vol) uniformly
        sigma_new   = sigma * sigma_scale              ! new sigma (in box units)
        box_new     = box / sigma_scale                ! new box (in sigma units)
        density_new = density * sigma_scale**3         ! reduced density
        lj_scale    = [sigma_scale**12,sigma_scale**6] ! scaling for LJ potential
        pot_new     = pot * lj_scale                   ! scaled potential
        vir_new     = vir * lj_scale                   ! scaled virial

        delta =  ( SUM(pot_new-pot) + pressure_inp * ( box_new**3 - box**3 )  ) / temperature &
             &   + REAL(n+1) * LOG(density_new/density) ! factor (n+1) consistent with box scaling

        IF ( metropolis ( delta ) ) THEN ! accept because Metropolis test
           pot     = pot_new     ! update both parts of potential
           vir     = vir_new     ! update both parts of virial
           sigma   = sigma_new   ! update sigma
           box     = box_new     ! update box
           density = density_new ! update density
           box_move_ratio = 1.0  ! increment move counter
        END IF ! reject Metropolis test

        ! Calculate all variables for this step
        potential = SUM ( pot ) / REAL(n)
        pressure  = density * temperature + SUM ( vir ) / box**3
        CALL blk_add ( [move_ratio,box_move_ratio,density,potential,pressure] )

     END DO ! End loop over steps

     CALL blk_end ( blk )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end

  potential = SUM ( pot ) / REAL ( n )
  pressure  = density * temperature + SUM ( vir ) / box**3
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure
  WRITE(*,'(''Final density (sigma units)'',         t40,f15.5)') density

  CALL energy ( sigma, r_cut, overlap, pot2=pot, vir2=vir )
  IF ( overlap ) STOP 'Overlap in final configuration'
  CALL energy_lrc ( n, sigma, r_cut, pot2=pot_lrc, vir2=vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential = SUM ( pot ) / REAL ( n )
  pressure  = density * temperature + SUM ( vir ) / box**3
  WRITE(*,'(''Final check'')')
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure
  WRITE(*,'(''Final density (sigma units)'',         t40,f15.5)') REAL(n) * sigma **3

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )

  CALL finalize

END PROGRAM mc_npt_lj

