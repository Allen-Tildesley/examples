! mc_npt_lj.f90
! Monte Carlo, NPT ensemble
PROGRAM mc_npt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis, random_translate_vector
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       potential_1, potential, move, n, r, &
       &                       potential_type, OPERATOR(+), OPERATOR(-)
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
  REAL :: box         ! Box length
  REAL :: dr_max      ! Maximum MC particle displacement
  REAL :: db_max      ! Maximum MC box displacement
  REAL :: temperature ! Specified temperature
  REAL :: pressure    ! Specified pressure
  REAL :: r_cut       ! Potential cutoff distance

  ! Quantities to be averaged
  REAL :: m_ratio ! Acceptance ratio of moves
  REAL :: v_ratio ! Acceptance ratio of volume moves
  REAL :: density ! Density
  REAL :: p_c     ! Pressure for simulated, cut, potential
  REAL :: p       ! Pressure for full potential with LRC
  REAL :: en_c    ! Internal energy per atom for simulated, cut, potential
  REAL :: en      ! Internal energy per atom for full potential with LRC

  ! Composite interaction = pot & vir & overlap variables
  TYPE(potential_type) :: system, system_new, atom_old, atom_new

  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL               :: box_scale, box_new, den_scale, delta, zeta
  REAL, DIMENSION(3) :: ri

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, pressure, r_cut, dr_max, db_max

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_npt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NPT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Simulation uses cut (but not shifted) potential'
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

  system = potential ( box, r_cut ) ! Initial energy and overlap check
  IF ( system%overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_npt_lj'
  END IF
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Volume ratio', &
       &            'Density', 'E/N (cut)', 'P (cut)', 'E/N (full)', 'P (full)' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           atom_old = potential_1 ( r(:,i), i, box, r_cut ) ! Old atom potential, virial etc

           IF ( atom_old%overlap ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_npt_lj'
           END IF

           ri(:) = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )                        ! Periodic boundary correction

           atom_new = potential_1 ( ri, i, box, r_cut ) ! New atom potential, virial etc

           IF ( .NOT. atom_new%overlap ) THEN ! Test for non-overlapping configuration

              delta = ( atom_new%pot_c - atom_old%pot_c ) / temperature ! Use cut (but not shifted) potential

              IF (  metropolis ( delta )  ) THEN ! Accept Metropolis test
                 system = system + atom_new - atom_old ! Update system values
                 CALL move ( i, ri )                   ! Update position
                 moves = moves + 1                     ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End test for non-overlapping configuration

        END DO ! End loop over atoms

        m_ratio = REAL(moves) / REAL(n)

        v_ratio = 0.0                   ! Zero volume move counter
        CALL RANDOM_NUMBER ( zeta )     ! Uniform random number in range (0,1)
        zeta      = 2.0*zeta - 1.0      ! now in range (-1,+1)
        box_scale = EXP ( zeta*db_max ) ! Sampling log(box) and log(vol) uniformly
        box_new   = box * box_scale     ! New box (in sigma units)
        den_scale = 1.0 / box_scale**3  ! Density scaling factor

        system_new = potential ( box_new, r_cut ) ! New system energy, virial etc

        IF ( .NOT. system_new%overlap ) THEN ! Test for non-overlapping configuration

           delta = ( system_new%pot_c - system%pot_c + pressure * ( box_new**3 - box**3 )  ) / temperature &
                &   + REAL(n+1) * LOG(den_scale) ! Factor (n+1) consistent with log(box) sampling

           IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
              system  = system_new ! Update system values
              box     = box_new    ! Update box
              v_ratio = 1.0        ! Set move counter
           END IF ! reject Metropolis test

        END IF ! End test for overlapping configuration

        CALL calculate ( )
        CALL blk_add ( [m_ratio,v_ratio,density,en_c,p_c,en,p] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk        ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  system = potential ( box, r_cut )
  IF ( system%overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_npt_lj'
  END IF
  CALL calculate ( 'Final check' )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    USE mc_module, ONLY : potential_lrc, pressure_lrc, pressure_delta
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates the properties of interest from system values
    ! In this example we simulate using the cut (but not shifted) potential
    ! Accordingly, < p_c > should match the input pressure and the values
    ! of < p_c >,  < en_c > and < density > should be consistent (for this potential)
    ! For comparison, long-range corrections are also applied to give
    ! estimates of < en > and < p > for the full (uncut) potential
    ! The value of the cut-and-shifted potential pot_s is not used, in this example
    
    density = REAL(n) / box**3                        ! Number density N/V
    en_c    = system%pot_c / REAL ( n )               ! PE/N for cut (but not shifted) potential
    en_c    = en_c + 1.5 * temperature                ! Add ideal gas contribution KE/N to give E_c/N
    en      = en_c + potential_lrc ( density, r_cut ) ! Add long-range contribution to give E/N estimate
    p_c     = system%vir / box**3                     ! Virial contribution to P_c
    p_c     = p_c + density * temperature             ! Add ideal gas contribution to P_c
    p       = p_c + pressure_lrc ( density, r_cut )   ! Add long-range contribution to give P
    p_c     = p_c + pressure_delta ( density, r_cut ) ! Add delta correction to P_c (not needed for P)

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density',    density
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (cut)',  en_c
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (cut)',    p_c
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (full)', en
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (full)',   p
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_npt_lj

