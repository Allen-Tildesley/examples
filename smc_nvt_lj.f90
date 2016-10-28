! smc_nvt_lj.f90
! Smart Monte Carlo, NVT ensemble
PROGRAM smc_nvt_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : random_normals, metropolis
  USE smc_module,       ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, force_1, r, r_old, zeta, v, move, n, &
       &                       potential_type

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
  REAL    :: box         ! Box length
  REAL    :: density     ! Density
  REAL    :: temperature ! Temperature (specified)
  INTEGER :: move_mode   ! Selects single- or multi-atom moves
  REAL    :: fraction    ! Fraction of atoms to move in multi-atom move
  REAL    :: dt          ! Time step
  REAL    :: r_cut       ! Potential cutoff distance

  ! Quantities to be averaged
  real :: m_ratio ! Acceptance ratio for moves
  REAL :: en_s    ! Internal energy per atom for simulated, cut-and-shifted, potential
  REAL :: en_f    ! Internal energy per atom for full potential with LRC
  REAL :: p_s     ! Pressure for simulated, cut-and-shifted, potential
  REAL :: p_f     ! Pressure for full potential with LRC
  real :: tc      ! Configurational temperature

  ! Composite interaction = forces & pot & cut & vir & lap & overlap variables
  TYPE(potential_type) :: total, total_old, partial_old, partial_new

  INTEGER :: blk, stp, nstep, nblock, ioerr
  INTEGER :: i, n_move
  REAL    :: v_rms, kin_old, kin_new, delta

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
  v_rms = sqrt ( temperature ) ! RMS value for velocity selection

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density

  CALL allocate_arrays ( box, r_cut )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! Second call gets r

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! First calculation of total system properties and check for overlap
  total = force ( box, r_cut )
  IF ( total%overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in smc_nvt_lj'
  END IF
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move Ratio', 'E/N (cut&shift)', 'P (cut&shift)', &
       &            'E/N (full)', 'P (full)', 'T (con)' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        SELECT CASE ( move_mode )

        CASE ( single_atom )

           n_move = 0
           DO i = 1, n ! Loop over atoms
              r_old(:,i)  = r(:,i)                    ! Store old position of this atom
              partial_old = force_1 ( i, box, r_cut ) ! Old force and pot etc for this atom

              IF ( partial_old%overlap ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in smc_nvt_lj'
              END IF

              CALL random_normals ( 0.0, v_rms, v(:,i) )           ! Choose 3 random momentum components
              kin_old     = 0.5*SUM(v(:,i)**2)                     ! Old kinetic energy of this atom
              v(:,i)      = v(:,i) + 0.5 * dt * partial_old%f(:,i) ! Kick half-step for one atom with old force
              r(:,i)      = r(:,i) + dt * v(:,i) / box             ! Drift step (positions in box=1 units)
              r(:,i)      = r(:,i) - ANINT ( r(:,i) )              ! Periodic boundaries (box=1 units)
              partial_new = force_1 ( i, box, r_cut )              ! New force and pot etc for this atom

              IF ( partial_new%overlap ) THEN ! Test for overlap
                 r(:,i) = r_old(:,i) ! Restore position: this move is rejected
              ELSE
                 v(:,i)  = v(:,i) + 0.5 * dt * partial_new%f(:,i) ! Kick half-step for one atom with new force
                 kin_new = 0.5*SUM(v(:,i)**2)                     ! New kinetic energy of this atom

                 delta = partial_new%pot - partial_old%pot ! Cut-and-shifted potential
                 delta = delta + kin_new - kin_old         ! Include kinetic energy change
                 delta = delta / temperature               ! Divide by temperature

                 IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                    total = total + partial_new - partial_old ! Update total values
                    n_move = n_move + 1                       ! Update move counter
                 ELSE
                    r(:,i) = r_old(:,i) ! Restore position: this move is rejected
                 END IF ! End accept Metropolis test

              END IF ! End test for overlap

           END DO ! End loop over atoms

           m_ratio = REAL(n_move) / REAL(n)

        CASE ( multi_atom )

           CALL RANDOM_NUMBER ( zeta )                         ! Select n uniform random numbers
           move = SPREAD ( zeta < fraction, dim=2, ncopies=3 ) ! Construct mask for moving atoms

           r_old     = r                         ! Store old positions
           total_old = total                     ! Store old totals
           CALL random_normals ( 0.0, v_rms, v ) ! Choose 3*n random momenta
           kin_old = 0.5*SUM(v**2)               ! Old kinetic energy
           WHERE ( move )
              v = v + 0.5 * dt * total%f(:,1:n) ! Kick half-step with old forces
              r = r + dt * v / box              ! Drift step (positions in box=1 units)
              r = r - ANINT ( r )               ! Periodic boundaries (box=1 units)
           END WHERE
           total = force ( box, r_cut ) ! New force and potential etc

           IF ( total%overlap ) THEN ! Test for overlap
              r       = r_old     ! Restore positions: this move is rejected
              total   = total_old ! Restore old totals
              m_ratio = 0.0       ! Set move counter
           ELSE
              WHERE ( move )
                 v = v + 0.5 * dt * total%f(:,1:n) ! Kick half-step with new forces
              END WHERE
              kin_new = 0.5*SUM(v**2) ! New kinetic energy

              delta = total%pot - total_old%pot ! Cut-and-shifted potential
              delta = delta + kin_new - kin_old ! Include kinetic energy change
              delta = delta / temperature       ! Divide by temperature

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 m_ratio  = 1.0 ! Set move counter
              ELSE
                 r       = r_old     ! Restore positions: this move is rejected
                 total   = total_old ! Restore old values
                 m_ratio = 0.0       ! Set move counter
              END IF ! End accept Metropolis test

           END IF ! End test for overlap

        END SELECT

        ! Calculate all variables for this step
        CALL calculate ( )
        CALL blk_add ( [m_ratio,en_s,p_s,en_f,p_f,tc] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  ! Double check book-keeping for totals, and final overlap
  total = force ( box, r_cut )
  IF ( total%overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in smc_nvt_lj'
  END IF
  CALL calculate ( 'Final check' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    USE smc_module, ONLY : potential_lrc, pressure_lrc
    IMPLICIT NONE
    CHARACTER (len=*), INTENT(in), OPTIONAL :: string

    en_s = total%pot / REAL ( n )   ! PE/N for cut-and-shifted potential
    en_s = en_s + 1.5 * temperature ! Add ideal gas contribution KE/N to give E_s/N

    en_f = total%cut / REAL ( n )                  ! PE/N for cut (but not shifted) potential
    en_f = en_f + 1.5 * temperature                ! Add ideal gas contribution KE/N
    en_f = en_f + potential_lrc ( density, r_cut ) ! Add long-range contribution to give E_f/N

    p_s = total%vir / box**3          ! Virial contribution to P
    p_s = p_s + density * temperature ! Add ideal gas contribution to give P_s

    p_f = p_s + pressure_lrc ( density, r_cut ) ! Add long-range contribution to give P_f

    tc = SUM ( total%f(:,1:n)**2 ) / total%lap

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (cut&shift)', en_s
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (cut&shift)',   p_s
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (full)',      en_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (full)',        p_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (con)',         tc
    END IF

  END SUBROUTINE calculate

END PROGRAM smc_nvt_lj

