! md_chain_mts_lj.f90
! Molecular dynamics, multiple timesteps, chain molecule
PROGRAM md_chain_mts_lj

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

  ! Takes in a configuration of atoms in a linear chain (positions, velocities)
  ! NO periodic boundary conditions, no box
  ! Conducts molecular dynamics with springs and multiple timesteps
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Input configuration, output configuration, all calculations, and all results 
  ! are given in mass = 1 units, and in simulation units defined by the model 
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     zero_cm, force, spring, worst_bond, r, v, f, g, n, potential_type

  IMPLICIT NONE

  ! Most important variables
  INTEGER :: n_mts     ! Number of small steps per large step
  REAL    :: dt        ! Time step (smallest)
  REAL    :: bond      ! Bond length
  REAL    :: k_spring  ! Bond spring constant
  REAL    :: total_spr ! Total spring harmonic potential energy

  ! Composite interaction = pot & ovr variables
  TYPE(potential_type) :: total

  INTEGER :: blk, stp, nstep, nblock, stp_mts, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, dt, k_spring, n_mts

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_chain_mts_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE ensemble, chain molecule, multiple time steps'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  WRITE ( unit=output_unit, fmt='(a)' ) 'No periodic boundaries'
  CALL introduction

  ! Set sensible default run parameters for testing
  nblock   = 10
  nstep    = 100000
  dt       = 0.0002
  k_spring = 10000.0
  n_mts    = 10

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_chain_mts_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Bond spring constant',      k_spring
  WRITE ( unit=output_unit, fmt='(a,t40,i15  )' ) 'Multiple time step factor', n_mts
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Large time step',           dt*n_mts

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond ) ! First call is just to get n and bond
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',          n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Bond length (in sigma units)', bond
  CALL allocate_arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r, v ) ! Second call is to get r and v
  CALL zero_cm ! Set centre-of-mass position and velocity to zero
  WRITE ( unit=output_unit, fmt='(a,t40,es15.5)' ) 'Worst bond length deviation = ', worst_bond ( bond )

  ! Initial calculation of forces f, spring forces g, and potential energies
  CALL force ( total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in md_chain_mts_lj'
  END IF
  CALL spring ( k_spring, bond, total_spr )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! Single time step of length n_mts*dt

        v = v + 0.5 * REAL(n_mts) * dt * f  ! Kick half-step (long) with nonbonded forces f

        DO stp_mts = 1, n_mts ! Loop over n_mts steps of length dt
           v = v + 0.5 * dt * g                      ! Kick half-step (short) with spring forces g
           r = r + dt * v                            ! Drift step (short)
           CALL spring ( k_spring, bond, total_spr ) ! Evaluate spring forces g and potential
           v = v + 0.5 * dt * g                      ! Kick half-step (short) with spring forces g
        END DO ! End loop over n_mts steps of length dt

        CALL force ( total ) ! Evaluate nonbonded forces f and potential
        IF ( total%ovr ) THEN
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Error in md_chain_mts_lj'
        END IF

        v = v + 0.5 * REAL(n_mts) * dt * f ! Kick half-step (long) with nonbonded forces f

        ! End single time step of length n_mts*dt

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                        ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  WRITE ( unit=output_unit, fmt='(a,t40,es15.5)' ) 'Worst bond length deviation = ', worst_bond ( bond )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r, v ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE averages_module, ONLY : variable_type, msd, cke
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(5) :: variables ! The 5 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: e_f, t_k, r_g, c_f, conserved_msd
    REAL                :: kin, eng
    REAL, DIMENSION(3)  :: rcm
    REAL                :: rsq

    ! Preliminary calculations
    kin = 0.5*SUM(v**2)
    rcm = SUM ( r, dim=2 ) / REAL(n)                                 ! Centre of mass
    rsq = SUM ( ( r - SPREAD(rcm,dim=2,ncopies=n) ) ** 2 ) / REAL(n) ! Mean-squared distance from CM
    eng = kin+total%pot+total_spr                                    ! Total energy

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Internal energy per atom
    ! Total KE plus total LJ nonbonded energy plus total spring energy divided by N
    e_f = variable_type ( nam = 'E/N', val = eng/REAL(n) )

    ! Kinetic temperature
    ! Remove 6 degrees of freedom for conserved linear and angular momentum
    t_k = variable_type ( nam = 'T kinetic', val = 2.0*kin/REAL(3*n-6) )

    ! Radius of gyration
    r_g = variable_type ( nam = 'Rg', val = SQRT(rsq) )

    ! MSD of kinetic energy, intensive
    ! Use special method to convert to Cv/N
    c_f = variable_type ( nam = 'Cv/N', val = kin/SQRT(REAL(n)), method = cke, instant = .FALSE. )

    ! Mean-squared deviation of conserved energy
    conserved_msd = variable_type ( nam = 'Conserved MSD', val = eng/REAL(n), &
         &                          method = msd, e_format = .TRUE., instant = .FALSE. )

    ! Collect together for averaging
    variables = [ e_f, t_k, r_g, c_f, conserved_msd ]

  END FUNCTION calc_variables

END PROGRAM md_chain_mts_lj

