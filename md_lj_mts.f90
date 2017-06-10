! md_lj_mts.f90
! Molecular dynamics, NVE, multiple timesteps
PROGRAM md_lj_mts

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

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using velocity Verlet algorithm
  ! Uses no special neighbour lists, for clarity

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! This program just illustrates the idea of splitting the non-bonded interactions
  ! using a criterion based on distance, for use in a MTS scheme
  ! This would hardly ever be efficient for a simple potential of the LJ kind alone

  ! This program uses mass = 1 throughout
  ! Unlike most other example programs, positions are not divided by box length at all
  ! Input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     force, r, v, f, n, potential_type

  IMPLICIT NONE

  ! Most important variables
  REAL :: box     ! Box length
  REAL :: dt1     ! Time step (smallest)
  REAL :: lambda  ! Healing length for switch function

  INTEGER, PARAMETER        :: k_max = 3 ! Number of shells
  REAL,    DIMENSION(k_max) :: r_cut     ! Cutoff distance for each shell
  INTEGER, DIMENSION(k_max) :: n_mts     ! Successive ratios of number of steps for each shell
  REAL,    DIMENSION(k_max) :: dt        ! Timestep for each shell
  REAL,    DIMENSION(k_max) :: vol_shell ! Volume of each shell

  ! Composite interaction = pot & cut & vir & lap & ovr variables for each shell
  TYPE(potential_type), DIMENSION(k_max) :: total

  INTEGER            :: blk, stp1, stp2, stp3, nstep, nblock, k, ioerr
  REAL, DIMENSION(3) :: vcm

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

  NAMELIST /nml/ nblock, nstep, r_cut, lambda, dt1, n_mts

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_lj_mts'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE ensemble, multiple time steps'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL introduction

  ! Set sensible default run parameters for testing
  nblock = 10
  nstep  = 6250
  r_cut  = [ 2.4, 3.5, 4.0 ]
  n_mts  = [ 1, 4, 2 ]
  dt1    = 0.002
  lambda = 0.1

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_lj_mts'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'      ) 'Number of blocks',           nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'      ) 'Number of steps per block',  nstep
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.6))' ) 'Potential cutoff distances', r_cut(:)
  WRITE ( unit=output_unit, fmt='(a,t40,*(i15))'   ) 'Multiple step ratios',       n_mts(:)
  IF ( n_mts(1) /= 1 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)' ) 'n_mts(1) must be 1', n_mts(1)
     STOP 'Error in md_lj_mts'
  END IF
  IF ( ANY ( n_mts <= 0 ) ) THEN
     WRITE ( unit=error_unit, fmt='(a,*(i15))' ) 'n_mts values must be positive', n_mts
     STOP 'Error in md_lj_mts'
  END IF
  DO k = 1, k_max
     dt(k) = PRODUCT(n_mts(1:k))*dt1 ! Define time steps cumulatively
     IF ( k == 1 ) THEN
        vol_shell(k) = (4.0/3.0)*pi * r_cut(k)**3
     ELSE
        IF ( r_cut(k)-r_cut(k-1) < lambda ) THEN
           WRITE ( unit=error_unit, fmt='(a,3f15.6)' ) 'r_cut interval error', r_cut(k-1), r_cut(k), lambda
           STOP 'Error in md_lj_mts'
        END IF
        vol_shell(k) = (4.0/3.0)*pi * ( r_cut(k)**3 - r_cut(k-1)**3 )
     END IF
  END DO
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.6))' ) 'Time step for each shell', dt(:)
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.6))' ) 'Volume of each shell',     vol_shell(:)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)'    ) 'Switching length lambda',  lambda

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n) / box**3
  IF ( r_cut(k_max) > box/2.0  ) THEN
     WRITE ( unit=error_unit, fmt='(a,f15.6)') 'r_cut(k_max) too large ', r_cut(k_max)
     STOP 'Error in md_lj_mts'
  END IF
  CALL allocate_arrays ( r_cut )
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call to get r and v
  r(:,:) = r(:,:) - ANINT ( r(:,:) / box ) * box            ! Periodic boundaries
  vcm(:) = SUM ( v(:,:), dim=2 ) / REAL(n)                  ! Centre-of mass velocity
  v(:,:) = v(:,:) - SPREAD ( vcm(:), dim = 2, ncopies = n ) ! Set COM velocity to zero

  ! Calculate initial forces and pot, vir contributions for each shell
  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, total(k) )
     IF ( total(k)%ovr ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
        STOP 'Error in md_lj_mts'
     END IF
  END DO

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     ! The following set of nested loops is specific to k_max=3

     DO stp3 = 1, nstep ! Begin loop over steps

        ! Outer shell 3: a single step of size n_mts(3)*n_mts(2)*dt1

        CALL kick_propagator ( 0.5*dt(3), 3 ) ! Kick half-step (outer shell)

        DO stp2 = 1, n_mts(3) ! Middle shell 2: n_mts(3) steps of size n_mts(2)*dt1

           CALL kick_propagator ( 0.5*dt(2), 2 ) ! Kick half-step (middle shell)

           DO stp1 = 1, n_mts(2) ! Inner shell 1: n_mts(3)*n_mts(2) steps of size dt1

              CALL kick_propagator ( 0.5*dt(1), 1 ) ! Kick half-step (inner shell)

              CALL drift_propagator ( dt(1) )        ! Drift step
              r(:,:) = r(:,:) - ANINT ( r(:,:)/box ) * box ! Periodic boundaries

              CALL force ( box, r_cut, lambda, 1, total(1) )
              IF ( total(1)%ovr ) THEN
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
                 STOP 'Error in md_lj_mts'
              END IF

              CALL kick_propagator ( 0.5*dt(1), 1 ) ! Kick half-step (inner shell)

           END DO ! End inner shell 1

           CALL force ( box, r_cut, lambda, 2, total(2) )
           IF ( total(2)%ovr ) THEN ! Highly unlikely for middle shell
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
              STOP 'Highly unlikely error in md_lj_mts'
           END IF

           CALL kick_propagator ( 0.5*dt(2), 2 ) ! Kick half-step (middle shell)

        END DO ! End middle shell 2

        CALL force ( box, r_cut, lambda, 3, total(3) )
        IF ( total(3)%ovr ) THEN ! Extremely unlikely for outer shell
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Extremely unlikely error in md_lj_mts'
        END IF

        CALL kick_propagator ( 0.5*dt(3), 3 ) ! Kick half-step (outer shell)

        ! End outer shell 3

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                       ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk           ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, total(k) )
     IF ( total(k)%ovr ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
        STOP 'Error in md_lj_mts'
     END IF
  END DO

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r, v ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE kick_propagator ( t, k )
    IMPLICIT NONE
    REAL,    INTENT(in) :: t ! Timestep (typically half the current timestep)
    INTEGER, INTENT(in) :: k ! Force array shell number

    v(:,:) = v(:,:) + t * f(:,:,k)

  END SUBROUTINE kick_propagator

  SUBROUTINE drift_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Timestep (typically dt)

    r(:,:) = r(:,:) + t * v(:,:) 

  END SUBROUTINE drift_propagator

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE lrc_module,      ONLY : potential_lrc, pressure_lrc
    USE md_module,       ONLY : hessian
    USE averages_module, ONLY : variable_type, msd, cke
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(8) :: variables ! The 8 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: e_s, p_s, e_f, p_f, t_k, t_c, c_s, conserved_msd

    REAL :: kin, vol, rho, tmp, pot, cut, eng, vir, lap, fsq, hes

    ! Preliminary calculations
    kin = 0.5*SUM(v**2)                  ! Total kinetic energy
    tmp = 2.0 * kin / REAL ( 3*n - 3 )   ! Three degrees of freedom for conserved momentum
    vol = box**3                         ! Volume
    rho = REAL ( n ) / box**3            ! Density
    pot = SUM ( total(:)%pot )           ! Sum cut-and-shifted potential over shells
    cut = SUM ( total(:)%cut )           ! Sum cut (but not shifted) potential over shells
    vir = SUM ( total(:)%vir )           ! Sum virial over shells
    lap = SUM ( total(:)%lap )           ! Sum Laplacian over shells
    fsq = SUM ( SUM(f(:,:,:),dim=3)**2 ) ! Sum forces over shells before squaring
    hes = hessian ( box, r_cut(k_max) )  ! Hessian (not resolved into shells)
    eng = kin + pot                      ! Total energy (should be conserved)

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Kinetic temperature
    t_k = variable_type ( nam = 'T kinetic', val = tmp )

    ! Internal energy (cut-and-shifted) per atom
    ! Total KE plus cut-and-shifted PE divided by N
    e_s = variable_type ( nam = 'E/N cut&shifted', val = eng/REAL(n) ) 

    ! Internal energy (full, including LRC) per atom
    ! LRC plus total KE plus total cut (but not shifted) PE divided by N
    e_f = variable_type ( nam = 'E/N full', val = potential_lrc(rho,r_cut(k_max)) + (kin+cut)/REAL(n) )

    ! Pressure (cut-and-shifted)
    ! Ideal gas contribution plus total virial divided by V
    p_s = variable_type ( nam = 'P cut&shifted', val = rho*tmp + vir/vol )

    ! Pressure (full, including LRC)
    ! LRC plus ideal gas contribution plus total virial divided by V
    p_f = variable_type ( nam = 'P full', val = pressure_lrc(rho,r_cut(k_max)) + rho*tmp + vir/vol )

    ! Configurational temperature
    ! Total squared force divided by Laplacian with small Hessian correction
    t_c = variable_type ( nam = 'T config', val = fsq/(lap-2.0*(hes/fsq)) )

    ! MSD of kinetic energy, intensive
    ! Use special method to convert to Cv/N
    c_s = variable_type ( nam = 'Cv/N cut&shifted', val = kin/SQRT(REAL(n)), method = cke, instant = .FALSE. )

    ! Mean-squared deviation of conserved energy
    conserved_msd = variable_type ( nam = 'Conserved MSD', val = eng/REAL(n), &
         &                          method = msd, e_format = .TRUE., instant = .FALSE. )

    ! Collect together for averaging
    variables = [ e_s, p_s, e_f, p_f, t_k, t_c, c_s, conserved_msd ]

  END FUNCTION calc_variables

END PROGRAM md_lj_mts

