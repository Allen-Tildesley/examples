! md_lj_mts.f90
! Molecular dynamics, NVE, multiple timesteps
PROGRAM md_lj_mts

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, r, v, f, n, potential_type

  IMPLICIT NONE

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

  ! Most important variables
  REAL :: box     ! Box length
  REAL :: dt      ! Time step (smallest)
  REAL :: lambda  ! Healing length for switch function

  INTEGER, PARAMETER        :: k_max = 3 ! Number of shells
  REAL,    DIMENSION(k_max) :: r_cut     ! Cutoff distance for each shell
  INTEGER, DIMENSION(k_max) :: n_mts     ! Successive ratios of number of steps for each shell

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  ! Composite interaction = pot & cut & vir & lap & ovr variables for each shell
  TYPE(potential_type), DIMENSION(k_max) :: total

  INTEGER            :: blk, stp1, stp2, stp3, nstep, nblock, k, ioerr
  REAL               :: pairs
  REAL, DIMENSION(3) :: vcm

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

  NAMELIST /nml/ nblock, nstep, r_cut, lambda, dt, n_mts

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_lj_mts'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE ensemble, multiple time steps'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL introduction
  CALL time_stamp

  ! Set sensible default run parameters for testing
  nblock = 10
  nstep  = 6000
  r_cut  = [ 2.4, 3.5, 4.0 ]
  n_mts  = [ 1, 4, 2 ]
  dt     = 0.002
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
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.6))' ) 'Potential cutoff distances', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,*(i15))'   ) 'Multiple step ratios', n_mts(:)
  IF ( n_mts(1) /= 1 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)' ) 'n_mts(1) must be 1', n_mts(1)
     STOP 'Error in md_lj_mts'
  END IF
  IF ( ANY ( n_mts <= 0 ) ) THEN
     WRITE ( unit=error_unit, fmt='(a,*(i15))' ) 'n_mts values must be positive', n_mts
     STOP 'Error in md_lj_mts'
  END IF
  DO k = 1, k_max
     WRITE ( unit=output_unit, fmt='(a,i1,t40,f15.6)' ) 'Time step for shell ', k, PRODUCT(n_mts(1:k))*dt
  END DO

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

  ! Diagnostic information to assess efficiency
  DO k = 1, k_max
     IF ( k == 1 ) THEN
        pairs = r_cut(k)**3
     ELSE
        pairs = r_cut(k)**3 - r_cut(k-1)**3
        IF ( r_cut(k)-r_cut(k-1) < lambda ) THEN
           WRITE ( unit=error_unit, fmt='(a,3f15.6)' ) 'r_cut interval error', r_cut(k-1), r_cut(k), lambda
           STOP 'Error in md_lj_mts'
        END IF
     END IF
     pairs = REAL((n*(n-1))/2) * (4.0/3.0)*pi * pairs / box**3
     WRITE ( unit=output_unit, fmt='(a,i1,t40,i15)' ) 'Estimated pairs in shell ', k, NINT ( pairs )
  END DO

  ! Calculate initial forces and pot, vir contributions for each shell
  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, total(k) )
     IF ( total(k)%ovr ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
        STOP 'Error in md_lj_mts'
     END IF
  END DO
  CALL calculate ( 'Initial values' )
  
  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( variables )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     ! The following set of nested loops is specific to k_max=3

     DO stp3 = 1, nstep ! Begin loop over steps

        ! Outer shell 3: a single step of size n_mts(3)*n_mts(2)*dt

        CALL kick_propagator ( 0.5*n_mts(3)*n_mts(2)*dt, 3 ) ! Kick half-step (outer shell)

        DO stp2 = 1, n_mts(3) ! Middle shell 2: n_mts(3) steps of size n_mts(2)*dt

           CALL kick_propagator ( 0.5*n_mts(2)*dt, 2 ) ! Kick half-step (middle shell)

           DO stp1 = 1, n_mts(2) ! Inner shell 1: n_mts(3)*n_mts(2) steps of size dt

              CALL kick_propagator ( 0.5*dt, 1 ) ! Kick half-step (inner shell)

              CALL drift_propagator ( dt )                 ! Drift step
              r(:,:) = r(:,:) - ANINT ( r(:,:)/box ) * box ! Periodic boundaries

              CALL force ( box, r_cut, lambda, 1, total(1) )
              IF ( total(1)%ovr ) THEN
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
                 STOP 'Error in md_lj_mts'
              END IF

              CALL kick_propagator ( 0.5*dt, 1 ) ! Kick half-step (inner shell)

           END DO ! End inner shell 1

           CALL force ( box, r_cut, lambda, 2, total(2) )
           IF ( total(2)%ovr ) THEN ! Highly unlikely for middle shell
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
              STOP 'Highly unlikely error in md_lj_mts'
           END IF

           CALL kick_propagator ( 0.5*n_mts(2)*dt, 2 ) ! Kick half-step (middle shell)

        END DO ! End middle shell 2

        CALL force ( box, r_cut, lambda, 3, total(3) )
        IF ( total(3)%ovr ) THEN ! Extremely unlikely for outer shell
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Extremely unlikely error in md_lj_mts'
        END IF

        CALL kick_propagator ( 0.5*n_mts(3)*n_mts(2)*dt, 3 ) ! Kick half-step (outer shell)

        ! End outer shell 3

        ! Calculate and accumulate variables for this step
        CALL calculate ( )
        CALL blk_add ( variables )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                       ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk           ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ! Output run averages

  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, total(k) )
     IF ( total(k)%ovr ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
        STOP 'Error in md_lj_mts'
     END IF
  END DO
  CALL calculate ( 'Final values' )
  CALL time_stamp

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

  SUBROUTINE calculate ( string ) 
    USE md_module,       ONLY : potential_lrc, pressure_lrc, hessian
    USE averages_module, ONLY : write_variables, msd, cke
    IMPLICIT NONE
    CHARACTER (len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: e_s, p_s, e_f, p_f, t_k, t_c, c_s, conserved_msd

    REAL :: kin, vol, rho, tmp, pot, cut, vir, lap, fsq, hes

    ! Preliminary calculations
    kin = 0.5*SUM(v**2)                      ! Total kinetic energy
    tmp = 2.0 * kin / REAL ( 3*n - 3 )       ! Three degrees of freedom for conserved momentum
    vol = box**3                             ! Volume
    rho = REAL ( n ) / box**3                ! Density
    pot = SUM ( total(:)%pot )               ! Sum cut-and-shifted potential over shells
    cut = SUM ( total(:)%cut )               ! Sum cut (but not shifted) potential over shells
    vir = SUM ( total(:)%vir )               ! Sum virial over shells
    lap = SUM ( total(:)%lap )               ! Sum Laplacian over shells
    fsq = SUM ( SUM ( f(:,:,:), dim=3 )**2 ) ! Sum forces over shells before squaring
    hes = hessian ( box, r_cut(k_max) )      ! Hessian (not resolved into shells)

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
    e_s = variable_type ( nam = 'E/N cut&shifted', val = (kin+pot)/REAL(n) ) 

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
    c_s = variable_type ( nam = 'Cv/N cut&shifted', val = kin/SQRT(REAL(n)), method = cke )

    ! Mean-squared deviation of conserved energy (not divided by N)
    conserved_msd = variable_type ( nam = 'Conserved MSD', val = kin+pot, method = msd )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ e_s, p_s, e_f, p_f, t_k, t_c, c_s, conserved_msd ]

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       CALL write_variables ( variables(1:6) ) ! Don't write out MSD variables
    END IF

  END SUBROUTINE calculate

END PROGRAM md_lj_mts

