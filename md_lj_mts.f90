! md_lj_mts.f90
! Molecular dynamics, NVE, multiple timesteps
PROGRAM md_lj_mts

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, r, v, f, n

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
  REAL :: density ! Density
  REAL :: dt      ! Time step (smallest)
  REAL :: lambda  ! Healing length for switch function

  INTEGER, PARAMETER        :: k_max = 3 ! Number of shells
  REAL,    DIMENSION(k_max) :: r_cut     ! Cutoff distance for each shell
  REAL,    DIMENSION(k_max) :: pot       ! Total cut-and-shifted potential energy for each shell
  REAL,    DIMENSION(k_max) :: cut       ! Total cut (but not shifted) potential energy for each shell
  REAL,    DIMENSION(k_max) :: vir       ! Total virial for each shell
  REAL,    DIMENSION(k_max) :: lap       ! Total Laplacian for each shell
  INTEGER, DIMENSION(k_max) :: n_mts     ! Successive ratios of number of steps for each shell

  ! Quantities to be averaged
  REAL :: en_s ! Internal energy (cut-and-shifted) per atom
  REAL :: p_s  ! Pressure (cut-and-shifted)
  REAL :: en_f ! Internal energy (full, including LRC) per atom
  REAL :: p_f  ! Pressure (full, including LRC)
  REAL :: tk   ! Kinetic temperature
  REAL :: tc   ! Configurational temperature

  INTEGER :: blk, stp1, stp2, stp3, nstep, nblock, k, ioerr
  REAL    :: pairs
  LOGICAL :: overlap

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

  NAMELIST /nml/ nblock, nstep, r_cut, lambda, dt, n_mts

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_lj_mts'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE ensemble, multiple time steps'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock = 10
  nstep  = 1000
  r_cut  = [ 2.4, 3.5, 4.0 ]
  n_mts  = [ 1, 4, 2 ]
  dt     = 0.002
  lambda = 0.1

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_lj_mts'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'      ) 'Number of blocks',           nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'      ) 'Number of steps per block',  nstep
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.5))' ) 'Potential cutoff distances', r_cut

  DO k = 1, k_max
     IF ( k == 1 ) THEN
        pairs = r_cut(k)**3
     ELSE
        pairs = r_cut(k)**3 - r_cut(k-1)**3
        IF ( r_cut(k)-r_cut(k-1) < lambda ) THEN
           WRITE ( unit=error_unit, fmt='(a,3f15.5)' ) 'r_cut interval error', r_cut(k-1), r_cut(k), lambda
           STOP 'Error in md_lj_mts'
        END IF
     END IF
     pairs = REAL(n*(n-1)/2) * (4.0/3.0)*pi * pairs / box**3
     WRITE ( unit=output_unit, fmt='(a,i1,t40,i15)' ) 'Estimated pairs in shell ', k, NINT ( pairs )
  END DO

  WRITE ( unit=output_unit, fmt='(a,t40,*(i15))'  ) 'Multiple step ratios', n_mts(:)
  IF ( n_mts(1) /= 1 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)' ) 'n_mts(1) must be 1', n_mts(1)
     STOP 'Error in md_lj_mts'
  END IF
  IF ( ANY ( n_mts <= 0 ) ) THEN
     WRITE ( unit=error_unit, fmt='(a,*(i15))' ) 'n_mts values must be positive', n_mts
     STOP 'Error in md_lj_mts'
  END IF
  DO k = 1, k_max
     WRITE ( unit=output_unit, fmt='(a,i1,t40,f15.5)' ) 'Time step for shell ', k, PRODUCT(n_mts(1:k))*dt
  END DO

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box ** 3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density
  IF ( r_cut(k_max) > box/2.0  ) THEN
     WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_cut(k_max) too large ', r_cut(k_max)
     STOP 'Error in md_lj_mts'
  END IF

  CALL allocate_arrays ( r_cut )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call to get r and v

  r(:,:) = r(:,:) - ANINT ( r(:,:) / box ) * box ! Periodic boundaries

  ! Calculate initial forces and pot, vir contributions for each shell
  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, pot(k), cut(k), vir(k), lap(k), overlap )
     IF ( overlap ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
        STOP 'Error in md_lj_mts'
     END IF
  END DO
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'E/N (cut&shift)', 'P (cut&shift)', &
       &            'E/N (full)', 'P (full)', 'T (kin)', 'T (con)' ] )

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

              CALL force ( box, r_cut, lambda, 1, pot(1), cut(1), vir(1), lap(1), overlap )
              IF ( overlap ) THEN
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
                 STOP 'Error in md_lj_mts'
              END IF

              CALL kick_propagator ( 0.5*dt, 1 ) ! Kick half-step (inner shell)

           END DO ! End inner shell 1

           CALL force ( box, r_cut, lambda, 2, pot(2), cut(2), vir(2), lap(2), overlap )
           IF ( overlap ) THEN ! Highly unlikely for middle shell
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
              STOP 'Highly unlikely error in md_lj_mts'
           END IF

           CALL kick_propagator ( 0.5*n_mts(2)*dt, 2 ) ! Kick half-step (middle shell)

        END DO ! End middle shell 2

        CALL force ( box, r_cut, lambda, 3, pot(3), cut(3), vir(3), lap(3), overlap )
        IF ( overlap ) THEN ! Extremely unlikely for outer shell
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Extremely unlikely error in md_lj_mts'
        END IF

        CALL kick_propagator ( 0.5*n_mts(3)*n_mts(2)*dt, 3 ) ! Kick half-step (outer shell)

        ! End outer shell 3

        CALL calculate ( )
        CALL blk_add ( [en_s,p_s,en_f,p_f,tk,tc] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk           ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, pot(k), cut(k), vir(k), lap(k), overlap )
     IF ( overlap ) THEN
        WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
        STOP 'Error in md_lj_mts'
     END IF
  END DO
  CALL calculate ( 'Final values' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r, v )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE kick_propagator ( t, k )
    IMPLICIT NONE
    REAL,    INTENT(in) :: t ! Timestep (typically half the current timestep)
    INTEGER, INTENT(in) :: k ! Force array shell

    v(:,:) = v(:,:) + t * f(:,:,k)

  END SUBROUTINE kick_propagator

  SUBROUTINE drift_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Timestep (typically dt)

    r(:,:) = r(:,:) + t * v(:,:) 

  END SUBROUTINE drift_propagator

  SUBROUTINE calculate ( string ) 
    USE md_module, ONLY : potential_lrc, pressure_lrc, hessian
    IMPLICIT NONE
    CHARACTER (len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates variables of interest and (optionally) writes them out

    REAL :: fsq, beta, kin

    kin = 0.5*SUM(v**2)
    tk  = 2.0 * kin / REAL ( 3*(n-1) )

    ! Sum potential terms and virial over shells to get totals
    en_s = ( SUM(pot) + kin ) / REAL ( n ) 
    en_f = ( SUM(cut) + kin ) / REAL ( n ) + potential_lrc ( density, r_cut(k_max) )
    p_s  = density * tk + SUM(vir) / box**3
    p_f  = p_s + pressure_lrc ( density, r_cut(k_max) )

    fsq  = SUM ( SUM(f,dim=3)**2 ) ! Sum forces over shells before squaring
    beta = ( SUM(lap) / fsq ) - 2.0*hessian ( box, r_cut(k_max) ) / (fsq**2) ! include 1/N Hessian correction
    tc   = 1.0 / beta

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (cut&shift)', en_s
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (cut&shift)',   p_s
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (full)',      en_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (full)',        p_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (kin)',         tk
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (con)',         tc
    END IF

  END SUBROUTINE calculate

END PROGRAM md_lj_mts

