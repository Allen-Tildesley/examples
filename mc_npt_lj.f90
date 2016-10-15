! mc_npt_lj.f90
! Monte Carlo, NPT ensemble
PROGRAM mc_npt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       energy_1, energy, energy_lrc, move, n, r, potovr
  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts isothermal-isobaric Monte Carlo at the given temperature and pressure
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! For the LJ potential we could use the known scaling of the separate parts
  ! with distances (i.e. with box scaling) to handle the volume move.
  ! However, this would require us to scale the cutoff distance with the box
  ! We do not do this here; instead we simply recalculate the potential energy,
  ! keeping r_cut fixed (in simulation units)

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in mc_module

  ! The logarithm of the box length is sampled uniformly

  ! Most important variables
  REAL :: box            ! box length
  REAL :: dr_max         ! maximum MC particle displacement
  REAL :: db_max         ! maximum MC box displacement
  REAL :: temperature    ! specified temperature
  REAL :: pressure       ! specified pressure
  REAL :: r_cut          ! potential cutoff distance
  REAL :: pot            ! total potential energy
  REAL :: vir            ! total virial
  REAL :: move_ratio     ! acceptance ratio of moves (to be averaged)
  REAL :: box_move_ratio ! acceptance ratio of box moves (to be averaged)
  REAL :: density        ! density (to be averaged)
  REAL :: pres_virial    ! virial pressure (to be averaged)
  REAL :: potential      ! potential energy per atom (to be averaged)

  TYPE(potovr) :: eng_old, eng_new ! Composite energy = pot & vir & overlap variables

  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL               :: box_scale, box_new, den_scale, delta, zeta1
  REAL               :: pot_lrc, vir_lrc
  REAL, DIMENSION(3) :: ri   ! position of atom i
  REAL, DIMENSION(3) :: zeta ! random numbers

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, pressure, r_cut, dr_max, db_max

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_npt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NPT ensemble'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  pressure    = 0.1
  r_cut       = 2.5
  dr_max      = 0.15
  db_max      = 0.025
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_npt_lj'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Pressure',                  pressure
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum displacement',      dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum box displacement',  db_max

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! first call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'  ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Simulation box length', box

  CALL allocate_arrays ( box, r_cut ) ! Allocate r

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! second call is to get r

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  eng_old = energy ( box, r_cut )
  IF ( eng_old%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_npt_lj'
  END IF
  pot = eng_old%pot
  vir = eng_old%vir
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: &
       &            'Move ratio', 'Box ratio', 'Density', 'Potential', 'Virial Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           ri(:)   = r(:,i)
           eng_old = energy_1 ( ri, i, box, r_cut )

           IF ( eng_old%ovr ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_npt_lj'
           END IF

           CALL RANDOM_NUMBER ( zeta )           ! Three uniform random numbers in range (0,1)
           zeta    = 2.0*zeta - 1.0              ! now in range (-1,+1)
           ri(:)   = ri(:) + zeta * dr_max / box ! Trial move to new position (box=1 units)
           ri(:)   = ri(:) - ANINT ( ri(:) )     ! Periodic boundary correction
           eng_new = energy_1 ( ri, i, box, r_cut )

           IF ( .NOT. eng_new%ovr ) THEN ! Consider non-overlapping configuration

              delta = ( eng_new%pot - eng_old%pot ) / temperature

              IF (  metropolis ( delta )  ) THEN ! Accept Metropolis test
                 pot = pot + eng_new%pot - eng_old%pot   ! Update potential energy
                 vir = vir + eng_new%vir - eng_old%vir   ! Update virial
                 CALL move ( i, ri )                     ! Update position
                 moves = moves + 1                       ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End consider non-overlapping configuration

        END DO ! End loop over atoms

        move_ratio = REAL(moves) / REAL(n)

        box_move_ratio = 0.0
        CALL RANDOM_NUMBER ( zeta1 )          ! Uniform random number in range (0,1)
        zeta1     = 2.0*zeta1 - 1.0           ! now in range (-1,+1)
        box_scale = EXP ( zeta1*db_max )      ! Sampling log(box) and log(vol) uniformly
        box_new   = box * box_scale           ! New box (in sigma units)
        den_scale = 1.0 / box_scale**3        ! Density scaling factor
        eng_new   = energy ( box_new, r_cut ) ! New energy, virial etc

        IF ( .NOT. eng_new%ovr ) THEN ! Consider non-overlapping configuration

           delta = ( (eng_new%pot-pot) + pressure * ( box_new**3 - box**3 )  ) / temperature &
                &   + REAL(n+1) * LOG(den_scale) ! Factor (n+1) consistent with box scaling

           IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
              pot     = eng_new%pot ! Update potential
              vir     = eng_new%vir ! Update virial
              box     = box_new     ! Update box
              box_move_ratio = 1.0  ! Increment move counter
           END IF ! reject Metropolis test

        END IF ! reject overlapping configuration

        CALL calculate ( )
        CALL blk_add ( [move_ratio,box_move_ratio,density,potential,pres_virial] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk        ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  eng_old = energy ( box, r_cut )
  IF ( eng_old%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_npt_lj'
  END IF
  pot = eng_old%pot
  vir = eng_old%vir
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  CALL calculate ( 'Final check' )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    density     = REAL(n) / box**3
    potential   = ( pot + pot_lrc ) / REAL ( n )
    pres_virial = density * temperature + ( vir + vir_lrc ) / box**3

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential energy', potential
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Virial pressure',  pres_virial
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density',          density
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_npt_lj

