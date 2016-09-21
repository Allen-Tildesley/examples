! mc_nvt_lj.f90
! Monte Carlo, NVT ensemble
PROGRAM mc_nvt_lj
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis
  USE mc_module,        ONLY : model_description, allocate_arrays, deallocate_arrays, &
       &                       energy_1, energy, energy_lrc, move, n, r, ne

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in mc_module

  ! Most important variables
  REAL :: box         ! box length
  REAL :: density     ! density
  REAL :: dr_max      ! maximum MC displacement
  REAL :: temperature ! specified temperature
  REAL :: r_cut       ! potential cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: vir         ! total virial
  REAL :: move_ratio  ! acceptance ratio of moves (to be averaged)
  REAL :: pres_virial ! virial pressure (to be averaged)
  REAL :: potential   ! potential energy per atom (to be averaged)

  LOGICAL            :: overlap
  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL               :: pot_old, pot_new, pot_lrc, vir_old, vir_new, vir_lrc, delta
  REAL, DIMENSION(3) :: ri   ! position of atom i
  REAL, DIMENSION(3) :: zeta ! random numbers

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, r_cut, dr_max

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_nvt_lj'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT'
  CALL model_description ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize! random number generator

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  r_cut       = 2.5
  dr_max      = 0.15
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_lj'
  END IF

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum displacement',      dr_max

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! first call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density

  CALL allocate_arrays ( box, r_cut ) ! Allocate r

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! second call is to get r

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  CALL energy ( box, r_cut, overlap, pot, vir )
  IF ( overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_lj'
  END IF
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  pot = pot + pot_lrc
  vir = vir + vir_lrc
  potential   = pot / REAL ( n )
  pres_virial = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial potential energy', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial virial pressure',  pres_virial

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Potential', 'Virial Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

           ri(:) = r(:,i)
           CALL  energy_1 ( ri, i, ne, box, r_cut, overlap, pot_old, vir_old )
           IF ( overlap ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_nvt_lj'
           END IF
           ri(:) = ri(:) + zeta * dr_max / box ! trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )     ! periodic boundary correction
           CALL  energy_1 ( ri, i, ne, box, r_cut, overlap, pot_new, vir_new )

           IF ( .NOT. overlap ) THEN ! consider non-overlapping configuration
              delta = ( pot_new - pot_old ) / temperature
              IF ( metropolis ( delta ) ) THEN    ! accept Metropolis test
                 pot    = pot + pot_new - pot_old ! update potential energy
                 vir    = vir + vir_new - vir_old ! update virial
                 CALL move ( i, ri )              ! update position
                 moves  = moves + 1               ! increment move counter
              END IF ! reject Metropolis test
           END IF ! reject overlapping configuration

        END DO ! End loop over atoms

        ! Calculate all variables for this step
        move_ratio  = REAL(moves) / REAL(n)
        potential   = pot / REAL(n)
        pres_virial = density * temperature + vir / box**3
        CALL blk_add ( [move_ratio,potential,pres_virial] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk        ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  potential   = pot / REAL ( n )
  pres_virial = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final potential energy', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final virial pressure',  pres_virial

  CALL energy ( box, r_cut, overlap, pot, vir )
  IF ( overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_lj'
  END IF
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  pot        = pot + pot_lrc
  vir        = vir + vir_lrc
  potential  = pot / REAL ( n )
  pres_virial = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a)'           ) 'Final check'
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final potential energy', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final virial pressure',  pres_virial

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays

END PROGRAM mc_nvt_lj

