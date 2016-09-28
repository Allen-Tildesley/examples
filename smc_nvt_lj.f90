! smc_nvt_lj.f90
! Smart Monte Carlo, NVT ensemble
PROGRAM smc_nvt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : random_normals, metropolis
  USE smc_module,       ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, force1, energy_lrc, &
       &                       r, r_old, v, f, f1, move, n, ne

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Smart Monte Carlo using Hybrid Monte Carlo / Brownian Dynamics notation
  ! Uses no special neighbour lists
  ! Assume that a sweep consists of either
  ! (a) N successive single-particle moves
  ! (b) 1 multi-particle move involving a large fraction of atoms
  ! (large enough to justify calling the complete force routine)
  ! The ensemble corresponds to the shifted potential, not the simple cutoff potential

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in, and we assume mass=1 throughout
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  ! Most important variables
  REAL    :: box         ! box length
  REAL    :: density     ! density
  REAL    :: temperature ! temperature (specified)
  INTEGER :: move_mode   ! selects single- or multi-atom moves
  REAL    :: fraction    ! fraction of atoms to move in multi-atom move
  REAL    :: dt          ! time step
  REAL    :: r_cut       ! potential cutoff distance
  REAL    :: pot         ! total potential energy
  REAL    :: pot_sh      ! total shifted potential energy
  REAL    :: vir         ! total virial
  REAL    :: pres_virial ! virial pressure (to be averaged)
  REAL    :: energy      ! total energy per atom (to be averaged)
  REAL    :: energy_sh   ! total shifted energy per atom (to be averaged)

  INTEGER :: blk, stp, nstep, nblock, ioerr
  INTEGER :: i, n_move
  REAL    :: pot_lrc, vir_lrc, pot_old, pot_new, pot_sh_old, pot_sh_new, vir_old, vir_new
  REAL    :: kin_old, kin_new, delta, move_ratio

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number
  INTEGER,          PARAMETER :: single_atom = 1, multi_atom = 2

  NAMELIST /nml/ nblock, nstep, r_cut, dt, move_mode, temperature, fraction

  WRITE ( unit=output_unit, fmt='(a)' ) 'smc_nvt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Smart Monte Carlo, constant-NVT ensemble'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  r_cut       = 2.5
  dt          = 0.05
  temperature = 0.7
  move_mode   = 1
  fraction    = 0.5

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in smc_nvt_lj'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',               temperature
  IF ( move_mode == single_atom ) THEN
     WRITE ( unit=output_unit, fmt='(a,t40,a)' ) 'Move mode is ', 'single-atom'
  ELSE IF ( move_mode == multi_atom ) THEN
     WRITE ( unit=output_unit, fmt='(a,t40,a)'     ) 'Move mode is ', 'multi-atom'
     WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Fraction of atoms moving', fraction
     IF ( fraction < 0.0 .OR. fraction > 1.0 ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Error: fraction out of range'
        STOP 'Error in smc_nvt_lj'
     END IF
  ELSE
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error: move_mode out of range', move_mode
     STOP 'Error in smc_nvt_lj'
  END IF

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density

  CALL allocate_arrays ( box, r_cut )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! Second call gets r

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL force ( box, r_cut, pot, pot_sh, vir )
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  call calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move Ratio', 'Energy', 'Shifted Energy', 'Virial Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        SELECT CASE ( move_mode )

        CASE ( single_atom )

           n_move = 0
           DO i = 1, n ! loop over atoms
              r_old(:,i) = r(:,i)! Store old position of this atom
              CALL random_normals ( 0.0, SQRT(temperature), v(:,i) )          ! choose random momentum
              kin_old = 0.5*SUM(v(:,i)**2)                                    ! Compute old kinetic energy of this atom
              CALL force1 ( box, r_cut, i, ne, pot_old, pot_sh_old, vir_old ) ! old force f1 and potential etc
              v(:,i) = v(:,i) + 0.5 * dt * f1(:,i)                            ! Kick half-step for one atom
              r(:,i) = r(:,i) + dt * v(:,i) / box                             ! Drift step (positions in box=1 units)
              r(:,i) = r(:,i) - ANINT ( r(:,i) )                              ! Periodic boundaries (box=1 units)
              CALL force1 ( box, r_cut, i, ne, pot_new, pot_sh_new, vir_new ) ! new force f1 and potential etc
              v(:,i) = v(:,i) + 0.5 * dt * f1(:,i)                            ! Kick half-step for one atom
              kin_new = 0.5*SUM(v(:,i)**2)                                    ! Compute new kinetic energy of this atom

              delta = ( pot_sh_new - pot_sh_old + kin_new - kin_old ) / temperature
              IF ( metropolis ( delta ) ) THEN ! accept Metropolis test
                 pot    = pot + pot_new - pot_old          ! update potential energy
                 vir    = vir + vir_new - vir_old          ! update virial
                 pot_sh = pot_sh + pot_sh_new - pot_sh_old ! update shifted potential energy
                 n_move = n_move + 1                       ! update move counter
              ELSE
                 r(:,i) = r_old(:,i) ! restore position
              END IF ! reject Metropolis test

           END DO ! end loop over atoms
           move_ratio = REAL(n_move) / REAL(n)

        CASE ( multi_atom )

           CALL RANDOM_NUMBER ( v(1,:) )                         ! Select n uniform random numbers
           move = SPREAD ( v(1,:) < fraction, dim=2, ncopies=3 ) ! Construct mask for moving atoms

           r_old = r                                         ! Store old positions
           CALL random_normals ( 0.0, SQRT(temperature), v ) ! Choose random momenta
           kin_old    = 0.5*SUM(v**2,mask=move)              ! Compute old kinetic energy of moving atoms
           pot_old    = pot                                  ! Store old total potential energy
           vir_old    = vir                                  ! Store old total virial
           pot_sh_old = pot_sh                               ! Store old total shifted potential energy
           WHERE ( move )
              v = v + 0.5 * dt * f ! Kick half-step
              r = r + dt * v / box ! Drift step (positions in box=1 units)
              r = r - ANINT ( r )  ! Periodic boundaries (box=1 units)
           END WHERE
           CALL force ( box, r_cut, pot, pot_sh, vir ) ! Force evaluation, new total potential energy etc
           WHERE ( move ) v = v + 0.5 * dt * f         ! Kick half-step
           kin_new = 0.5*SUM(v**2,mask=move)           ! Compute new kinetic energy of moving atoms

           delta = ( pot_sh_new - pot_sh_old + kin_new - kin_old ) / temperature
           IF ( metropolis ( delta ) ) THEN  ! accept Metropolis test
              move_ratio = 1.0        ! set move counter
           ELSE
              r          = r_old      ! restore positions
              pot        = pot_old    ! restore potential energy
              vir        = vir_old    ! restore virial
              pot_sh     = pot_sh_old ! restore shifted potential energy
              move_ratio = 0.0        ! set move counter
           END IF ! reject Metropolis test

        END SELECT

        CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )

        ! Calculate all variables for this step
        call calculate ( )
        CALL blk_add ( [move_ratio,energy,energy_sh,pres_virial] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL force ( box, r_cut, pot, pot_sh, vir )
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  call calculate ( 'Final values' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS
  
  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER (len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates variables of interest and (optionally) writes them out
  energy      = 1.5*temperature + ( pot + pot_lrc ) / REAL ( n )
  energy_sh   = 1.5*temperature + pot_sh / REAL ( n )
  pres_virial = density * temperature + ( vir + vir_lrc ) / box**3

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Energy',          energy
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Shifted energy',  energy_sh
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Virial pressure', pres_virial
    END IF

  END SUBROUTINE calculate

END PROGRAM smc_nvt_lj

