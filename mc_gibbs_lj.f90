! mc_gibbs_lj.f90
! Monte Carlo, Gibbs ensemble
PROGRAM mc_gibbs_lj
  
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE maths_module,     ONLY : metropolis, random_integer, random_translate_vector
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       potential_1, potential, move, swap, n, r, potential_type

  IMPLICIT NONE

  ! Takes in a pair of configurations of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Gibbs ensemble Monte Carlo at the given temperature, total volume and total N
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Note that long-range corrections are not included in the acceptance/rejection
  ! of creation and destruction moves, but are added to the simulation averages

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in mc_module

  ! Most important variables
  REAL, dimension(2) :: box ! Box lengths
  REAL :: dr_max            ! Maximum MC displacement
  REAL :: temperature       ! Specified temperature
  REAL :: r_cut             ! Potential cutoff distance

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  ! Composite interaction = pot & vir & ovr variables
  TYPE(potential_type), dimension(2) :: total
  TYPE(potential_type)               :: partial_old, partial_new

  INTEGER            :: blk, stp, i, nstep, nblock
  INTEGER            :: try, ntry, m_try, m_acc, u_try, u_acc, d_try, d_acc, v_try, v_acc, ioerr
  REAL               :: prob_move, prob_create, delta, zeta, m_ratio, u_ratio, d_ratio, v_ratio
  REAL, DIMENSION(3) :: ri

  CHARACTER(len=5), PARAMETER :: cnf1_prefix = 'cnf1.'
  CHARACTER(len=5), PARAMETER :: cnf2_prefix = 'cnf2.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, r_cut, dr_max

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_gibbs_lj'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, Gibbs ensemble'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  temperature = 1.2
  r_cut       = 2.5
  dr_max      = 0.15

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_zvt_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',              nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',     nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature',                   temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Potential cutoff distance',     r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum displacement',          dr_max

  ! Read in initial configurations and allocate necessary arrays
  CALL read_cnf_atoms ( cnf1_prefix//inp_tag, n1, box(1) ) ! First call is just to get n and box
  CALL read_cnf_atoms ( cnf2_prefix//inp_tag, n2, box(2) ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n1
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box(1)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n1) / box(1)**3
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n2
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box(2)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n2) / box(2)**3
  CALL allocate_arrays ( box, r_cut ) ! Allocate r
  CALL read_cnf_atoms ( cnf1_prefix//inp_tag, n1, box(1), r(:,1:n1) )       ! Second call is to get r
  CALL read_cnf_atoms ( cnf2_prefix//inp_tag, n2, box(2), r(:,n1+1:n1+n2) ) ! Second call is to get r
  r1(:,1:n1)       = r(:,1:n1)       / box(1)  ! Convert positions to box units
  r1(:,n1+1:n1+n2) = r(:,n1+1:n1+n2) / box(2)  ! Convert positions to box units
  r(:,:)           = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial energy and overlap check
  total(1) = potential ( 1,    n1,    box(1), r_cut ) 
  total(2) = potential ( n1+1, n1+n2, box(2), r_cut ) 
  IF ( total(1)%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration 1'
     STOP 'Error in mc_gibbs_lj'
  END IF
  IF ( total(2)%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration 2'
     STOP 'Error in mc_gibbs_lj'
  END IF
  CALL calculate ( 'Initial values' )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( output_unit, variables )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        m_try = 0
        m_acc = 0
        u_try = 0
        u_acc = 0
        d_try = 0
        d_acc = 0
        v_try = 0
        v_acc = 0

        DO i = 1, n1 ! Loop over atoms in system 1
              partial_old = potential_1 ( 1, n1, r(:,i), i, box(1), r_cut ) ! Old atom potential, virial etc

              IF ( partial_old%ovr ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in mc_gibbs_lj'
              END IF

              ri(:) = random_translate_vector ( dr_max/box(1), r(:,i) ) ! Trial move to new position (in box=1 units)
              ri(:) = ri(:) - ANINT ( ri(:) )                           ! Periodic boundary correction

              partial_new = potential_1 ( 1, n1, ri, i, box(1), r_cut ) ! New atom potential, virial etc

              IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

                 delta = partial_new%pot - partial_old%pot ! Use cut (but not shifted) potential
                 delta = delta / temperature               ! Divide by temperature

                 IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                    CALL move ( i, ri )                       ! Update position
                    total(1) = total(1) + partial_new - partial_old ! Update total values
                    m_acc = m_acc + 1                         ! Increment move counter
                 END IF ! End accept Metropolis test

              END IF ! End test for overlapping configuration
           END DO ! End loop over atoms in system 1

           DO i = n1+1, n1+n2 ! Loop over atoms in system 2
              partial_old = potential_1 ( n1+1, n1+n2, r(:,i), i, box(2), r_cut ) ! Old atom potential, virial etc

              IF ( partial_old%ovr ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in mc_gibbs_lj'
              END IF

              ri(:) = random_translate_vector ( dr_max/box(2), r(:,i) ) ! Trial move to new position (in box=1 units)
              ri(:) = ri(:) - ANINT ( ri(:) )                           ! Periodic boundary correction

              partial_new = potential_1 ( n1+1, n1+n2, ri, i, box(2), r_cut ) ! New atom potential, virial etc

              IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

                 delta = partial_new%pot - partial_old%pot ! Use cut (but not shifted) potential
                 delta = delta / temperature               ! Divide by temperature

                 IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                    CALL move ( i, ri )                       ! Update position
                    total(2) = total(2) + partial_new - partial_old ! Update total values
                    m_acc = m_acc + 1                         ! Increment move counter
                 END IF ! End accept Metropolis test

              END IF ! End test for overlapping configuration
           END DO ! End loop over atoms in system 2

           ELSE IF ( zeta < prob_move + prob_create ) THEN ! Try create

              c_try = c_try + 1

              IF ( n+1 > SIZE(r,dim=2) ) then
                 WRITE ( unit=error_unit, fmt='(a,2i5)') 'n has grown too large', n+1, SIZE(r,dim=2)
                 STOP 'Error in mc_zvt_lj'
              END IF

              CALL RANDOM_NUMBER ( ri ) ! Three uniform random numbers in range (0,1)
              ri = ri - 0.5             ! Now in range (-0.5,+0.5) for box=1 units

              partial_new = potential_1 ( ri, n+1, box, r_cut ) ! New atom potential, virial, etc

              IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

                 delta = partial_new%pot / temperature                    ! Use cut (not shifted) potential
                 delta = delta - LOG ( activity * box**3 / REAL ( n+1 ) ) ! Activity term for creation

                 IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                    CALL create ( ri )          ! Create new particle
                    total = total + partial_new ! Update total values
                    c_acc = c_acc + 1           ! Increment creation move counter
                 END IF ! End accept Metropolis test

              END IF ! End test for overlapping configuration

           ELSE ! Try destroy

              d_try = d_try + 1

              i = random_integer ( 1, n ) ! Choose particle at random

              partial_old = potential_1 ( r(:,i), i, box, r_cut ) ! Old atom potential, virial, etc

              IF ( partial_old%ovr ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap found on particle removal'
                 STOP 'Error in mc_zvt_lj'
              END IF

              delta = -partial_old%pot / temperature                      ! Use cut (not shifted) potential
              delta = delta - LOG ( REAL ( n ) / ( activity * box**3 )  ) ! Activity term for destruction

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 CALL destroy ( i )          ! Destroy chosen particle
                 total = total - partial_old ! Update total values
                 d_acc = d_acc + 1           ! Increment destruction move counter
              END IF ! End accept Metropolis test

           END IF ! End choice of move type

        END DO ! End loop over tries of different kinds

        m_ratio = 0.0
        c_ratio = 0.0
        d_ratio = 0.0
        IF ( m_try > 0 ) m_ratio = REAL(m_acc) / REAL(m_try)
        IF ( c_try > 0 ) c_ratio = REAL(c_acc) / REAL(c_try)
        IF ( d_try > 0 ) d_ratio = REAL(d_acc) / REAL(d_try)

        ! Calculate and accumulate variables for this step
        CALL calculate ( )
        CALL blk_add ( variables )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )                                  ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk                   ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r(:,1:n)*box ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit ) ! Output run averages

  CALL calculate ( 'Final values' )

  ! Double-check book-keeping of totals and overlap
  total = potential ( box, r_cut )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_zvt_lj'
  END IF
  CALL calculate ( 'Final check' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r(:,1:n)*box ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    USE mc_module,       ONLY : potential_lrc, pressure_lrc, pressure_delta, force_sq
    USE averages_module, ONLY : write_variables
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    ! In this example we simulate using the cut (but not shifted) potential
    ! The values of < p_c >,  < e_c > and < density > should be consistent (for this potential)
    ! For comparison, long-range corrections are also applied to give
    ! estimates of < e_f > and < p_f > for the full (uncut) potential
    ! The value of the cut-and-shifted potential is not used, in this example

    TYPE(variable_type) :: m_r, c_r, d_r, density, e_c, p_c, e_f, p_f, t_c
    REAL                :: fsq, vol, rho

    ! Preliminary calculations (m_ratio, total etc are known already)
    fsq = force_sq ( box, r_cut )
    vol = box**3
    rho = REAL(n) / vol

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %msd: indicating if mean squared deviation required
    ! If not set below, %msd adopts its default value of .false.
    ! The %msd and %nam components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move, creation, and destruction acceptance ratios

    IF ( PRESENT ( string ) ) THEN ! The ratio is meaningless in this case
       m_r = variable_type ( nam = 'Move ratio',    val = 0.0 )
       c_r = variable_type ( nam = 'Create ratio',  val = 0.0 )
       d_r = variable_type ( nam = 'Destroy ratio', val = 0.0 )
    ELSE
       m_r = variable_type ( nam = 'Move ratio',    val = m_ratio )
       c_r = variable_type ( nam = 'Create ratio',  val = c_ratio )
       d_r = variable_type ( nam = 'Destroy ratio', val = d_ratio )
    END IF

    ! Density
    density = variable_type ( nam = 'Density', val = rho )

    ! Internal energy per atom for simulated, cut, potential
    ! Ideal gas contribution plus cut (but not shifted) PE divided by N
    e_c = variable_type ( nam = 'E/N cut', val = 1.5*temperature + total%pot/REAL(n) )

    ! Internal energy per atom for full potential with LRC
    ! LRC plus ideal gas contribution plus cut (but not shifted) PE divided by N
    e_f = variable_type ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + 1.5*temperature + total%pot/REAL(n) )

    ! Pressure for simulated, cut, potential
    ! delta correction plus ideal gas contribution plus total virial divided by V
    p_c = variable_type ( nam = 'P cut', val = pressure_delta(rho,r_cut) + rho*temperature + total%vir/vol )

    ! Pressure for full potential with LRC
    ! LRC plus ideal gas contribution plus total virial divided by V 
    p_f = variable_type ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*temperature + total%vir/vol )

    ! Configurational temperature
    ! Total squared force divided by total Laplacian
    t_c = variable_type ( nam = 'T config', val = fsq/total%lap )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ m_r, c_r, d_r, density, e_c, p_c, e_f, p_f, t_c ]

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       CALL write_variables ( output_unit, variables(4:) ) ! Don't write out move ratios
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_gibbs_lj
