! smc_nvt_lj.f90
! Smart Monte Carlo, NVT ensemble
PROGRAM smc_nvt_lj

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             !
  ! to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     !
  ! published by Oxford University Press ("the publishers").                                       !
  !                                                                                                !
  ! LICENCE                                                                                        !
  ! Creative Commons CC0 Public Domain Dedication.                                                 !
  ! To the extent possible under law, the authors have dedicated all copyright and related         !
  ! and neighboring rights to this software to the PUBLIC domain worldwide.                        !
  ! This software is distributed without any warranty.                                             !
  ! You should have received a copy of the CC0 Public Domain Dedication along with this software.  !
  ! If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               !
  !                                                                                                !
  ! DISCLAIMER                                                                                     !
  ! The authors and publishers make no warranties about the software, and disclaim liability       !
  ! for all uses of the software, to the fullest extent permitted by applicable law.               !
  ! The authors and publishers do not recommend use of this software for any purpose.              !
  ! It is made freely available, solely to clarify points made in the text. When using or citing   !
  ! the software, you should not imply endorsement by the authors or publishers.                   !
  !------------------------------------------------------------------------------------------------!

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
  ! The model is defined in smc_module

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               maths_module,     ONLY : random_normals, metropolis, lowercase
  USE               smc_module,       ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     force, force_1, r, r_old, zeta, v, move, n, potential_type

  IMPLICIT NONE

  ! Most important variables
  REAL    :: box         ! Box length
  REAL    :: temperature ! Temperature (specified)
  LOGICAL :: single_atom ! Selects single- or multi-atom moves
  REAL    :: fraction    ! Fraction of atoms to move in multi-atom move
  REAL    :: dt          ! Time step (effectively determines typical displacement)
  REAL    :: r_cut       ! Potential cutoff distance

  ! Composite interaction = forces & pot & cut & vir & lap & ovr variables
  TYPE(potential_type) :: total, total_old, partial_old, partial_new

  INTEGER :: blk, stp, nstep, nblock, ioerr
  INTEGER :: i, n_move
  REAL    :: v_rms, kin_old, kin_new, delta, m_ratio

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, r_cut, dt, single_atom, temperature, fraction

  WRITE ( unit=output_unit, fmt='(a)' ) 'smc_nvt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Smart Monte Carlo, constant-NVT ensemble'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 10000
  r_cut       = 2.5
  temperature = 1.0    ! Default temperature T
  dt          = 0.1    ! Together with v_rms=sqrt(T) determines typical displacement
  single_atom = .TRUE. ! .false. selects multi-atom mode, probably requiring smaller dt
  fraction    = 1.0    ! Only applicable in multi-atom mode

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in smc_nvt_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature',               temperature
  IF ( single_atom ) THEN
     WRITE ( unit=output_unit, fmt='(a,t40,a15)' ) 'Move mode is ', 'single-atom'
  ELSE
     WRITE ( unit=output_unit, fmt='(a,t40,a15)'   ) 'Move mode is ', 'multi-atom'
     WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Fraction of atoms moving', fraction
     IF ( fraction < 0.0 .OR. fraction > 1.0 ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Error: fraction out of range'
        STOP 'Error in smc_nvt_lj'
     END IF
  END IF
  v_rms = SQRT ( temperature ) ! RMS value for velocity selection
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Typical dr', v_rms*dt

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n) / box**3
  CALL allocate_arrays ( box, r_cut )
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! Second call gets r
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial forces, potential, etc plus overlap check
  total = force ( box, r_cut )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in smc_nvt_lj'
  END IF

  ! Initialize arrays for averaging and write column headings
  m_ratio = 0.0
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        IF ( single_atom ) THEN ! Single-atom moves

           n_move = 0
           DO i = 1, n ! Loop over atoms
              r_old(:,i)  = r(:,i)                    ! Store old position of this atom
              partial_old = force_1 ( i, box, r_cut ) ! Old force and pot etc for this atom

              IF ( partial_old%ovr ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in smc_nvt_lj'
              END IF

              CALL random_normals ( 0.0, v_rms, v(:,i) )           ! Choose 3 random momentum components
              kin_old     = 0.5*SUM(v(:,i)**2)                     ! Old kinetic energy of this atom
              v(:,i)      = v(:,i) + 0.5 * dt * partial_old%f(:,i) ! Kick half-step for one atom with old force
              r(:,i)      = r(:,i) + dt * v(:,i) / box             ! Drift step (positions in box=1 units)
              r(:,i)      = r(:,i) - ANINT ( r(:,i) )              ! Periodic boundaries (box=1 units)
              partial_new = force_1 ( i, box, r_cut )              ! New force and pot etc for this atom

              IF ( partial_new%ovr ) THEN ! Test for overlap
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

        ELSE ! Multi-atom moves

           CALL RANDOM_NUMBER ( zeta )                         ! Select N uniform random numbers
           move = SPREAD ( zeta < fraction, dim=1, ncopies=3 ) ! Construct mask for moving atoms

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

           IF ( total%ovr ) THEN ! Test for overlap
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

        END IF

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                        ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  ! Double check final overlap
  total = force ( box, r_cut )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in smc_nvt_lj'
  END IF

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box )

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE lrc_module,      ONLY : potential_lrc, pressure_lrc
    USE averages_module, ONLY : variable_type, msd
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(8) :: variables ! The 8 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: m_r, e_s, p_s, e_f, p_f, t_c, c_s, c_f
    REAL                :: vol, rho, fsq

    ! Preliminary calculations
    vol = box**3                    ! Volume
    rho = REAL(n) / vol             ! Density
    fsq = SUM ( total%f(:,1:n)**2 ) ! Total squared force

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move acceptance ratio
    m_r = variable_type ( nam = 'Move ratio', val = m_ratio, instant = .FALSE. )

    ! Internal energy per atom for simulated, cut-and-shifted, potential
    ! Ideal gas contribution plus total cut-and-shifted PE divided by N
    e_s = variable_type ( nam = 'E/N cut&shifted', val = 1.5*temperature + total%pot/REAL(n) )

    ! Internal energy per atom for full potential with LRC
    ! LRC plus ideal gas contribution plus total cut (but not shifted) PE divided by N
    e_f = variable_type ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + 1.5*temperature + total%cut/REAL(n) )

    ! Pressure for simulated, cut-and-shifted, potential
    ! Ideal gas contribution plus total virial divided by V
    p_s = variable_type ( nam = 'P cut&shifted', val = rho*temperature + total%vir/vol )

    ! Pressure for full potential with LRC
    ! LRC plus ideal gas contribution plus total virial divided by V
    p_f = variable_type ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*temperature + total%vir/vol )

    ! Configurational temperature
    ! Total squared force divided by total Laplacian
    t_c = variable_type ( nam = 'T config', val = fsq/total%lap )

    ! Heat capacity (excess, cut-and-shifted)
    ! Total PE divided by temperature and sqrt(N) to make result intensive
    ! Ideal gas contribution 1.5 added afterwards
    c_s = variable_type ( nam = 'Cv/N cut&shifted', val = total%pot/(temperature*SQRT(REAL(n))), &
         &                method = msd, add = 1.5, instant = .FALSE. )

    ! Heat capacity (excess, full)
    ! Total PE divided by temperature and sqrt(N) to make result intensive; LRC does not contribute
    ! Ideal gas contribution 1.5 added afterwards
    c_f = variable_type ( nam = 'Cv/N full', val = total%cut/(temperature*SQRT(REAL(n))), &
         &                method = msd, add = 1.5, instant = .FALSE. )

    ! Collect together for averaging
    variables = [ m_r, e_s, p_s, e_f, p_f, t_c, c_s, c_f ]

  END FUNCTION calc_variables

END PROGRAM smc_nvt_lj

