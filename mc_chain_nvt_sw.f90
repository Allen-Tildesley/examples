! mc_chain_nvt_sw.f90
! Monte Carlo, single chain, NVT, square wells
PROGRAM mc_chain_nvt_sw
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       regrow, crank, pivot, qcount, weight, n, r

  IMPLICIT NONE

  ! Takes in a configuration of atom positions in a linear chain
  ! NO periodic boundary conditions, no box
  ! Conducts Monte Carlo, NVT ensemble using various moves
  ! such as CBMC regrowth, pivot, crankshaft
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Input configuration, output configuration, all calculations, and all results 
  ! are given in simulation units defined by the model:
  ! atomic core diameter sigma = 1, well-depth epsilon=1
  ! Potential energy is -q, where q is the total number of attractive well interactions
  ! This is a negative integer, but for convenience we refer to q as energy
  ! Configurational weights are calculated on the basis of the athermal interactions
  ! only, i.e. weight=1 if no overlap, weight=0 if any overlap

  ! Most important variables
  REAL    :: temperature    ! temperature (specified, in units of well depth)
  REAL    :: bond           ! bond length (in units of core diameter)
  REAL    :: range          ! range of attractive well (specified, same units)
  INTEGER :: m_max          ! maximum atoms in regrow
  INTEGER :: k_max          ! number of random tries per atom in regrow
  REAL    :: crank_max      ! maximum move angle in crankshaft
  REAL    :: crank_fraction ! fraction of atoms to try in crankshaft per step
  INTEGER :: n_crank        ! number of crankshaft tries per step
  REAL    :: pivot_max      ! maximum move angle in pivot
  REAL    :: pivot_fraction ! fraction of atoms to try in pivot per step
  INTEGER :: n_pivot        ! number of pivot tries per step
  INTEGER :: q              ! total attractive interaction (negative of energy)
  INTEGER :: q_max          ! max value of q seen in simulation (minimum energy)
  INTEGER :: nq             ! Maximum anticipated energy

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  ! Histograms and tables
  REAL, DIMENSION(:), ALLOCATABLE :: h ! Histogram of q values (0:nq)
  REAL, DIMENSION(:), ALLOCATABLE :: s ! Entropy histogram used in acceptance (0:nq)

  INTEGER :: blk, stp, nstep, nblock, ioerr, try, n_acc
  LOGICAL :: accepted
  REAL    :: r_ratio, c_ratio, p_ratio

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=4), PARAMETER :: his_prefix = 'his.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, m_max, k_max, crank_max, crank_fraction, pivot_max, pivot_fraction, &
       &         temperature, range

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_chain_nvt_sw'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble, chain molecule, square wells'
  CALL introduction
  CALL time_stamp

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock         = 10    ! Number of blocks
  nstep          = 10000 ! Number of sweeps per block
  m_max          = 3     ! Maximum atoms in regrow
  k_max          = 32    ! Number of random tries per atom in regrow
  crank_max      = 0.5   ! Maximum move angle in crankshaft
  crank_fraction = 0.5   ! Fraction of atoms to try in crank moves
  pivot_max      = 0.5   ! Maximum move angle in pivot
  pivot_fraction = 0.2   ! Fraction of atoms to try in pivot moves
  temperature    = 1.0   ! Temperature (in units of well depth)
  range          = 1.5   ! Range of attractive well

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_chain_nvt_sw'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',                nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',       nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Max atoms in regrow',             m_max
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Random tries per atom in regrow', k_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Max move angle in crankshaft',    crank_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Crank fraction',                  crank_fraction
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Max move angle in pivot',         pivot_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Pivot fraction',                  pivot_fraction
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature/well depth',          temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Attractive well range',           range
  IF ( range < 1.0 ) THEN
     WRITE ( unit=output_unit, fmt='(a)' ) 'Warning, range < core diameter (1.0)'
  END IF

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond ) ! First call is just to get n and bond
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',          n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Bond length (in sigma units)', bond
  CALL allocate_arrays
  nq = 6*n ! Anticipated maximum number of pair interactions within range
  ALLOCATE ( h(0:nq), s(0:nq) )
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r ) ! Second call gets r

  ! Set number of crankshaft and pivot moves per step
  n_crank = NINT(crank_fraction*REAL(n))
  n_pivot = NINT(pivot_fraction*REAL(n))
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of crankshaft tries per step', n_crank
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of pivot tries per step',      n_pivot

  ! Initialize Boltzmann exponents, s(q), which are just values of E/kT
  s = [ ( -REAL(q)/temperature, q = 0, nq ) ]

  ! Initial energy calculation plus overlap check
  IF ( weight() == 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_chain_nvt_sw'
  END IF
  q = qcount ( range )
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Initial energy', q
  q_max = q ! Max q seen so far
  CALL calculate ( 'Initial values' )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( variables )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin
     h(:) = 0.0

     DO stp = 1, nstep ! Begin loop over steps

        r_ratio = 0.0
        CALL regrow ( s, m_max, k_max, bond, range, q, accepted )
        IF ( accepted ) r_ratio = 1.0
        CALL update_histogram

        p_ratio = 0.0
        IF ( n_pivot > 0 ) THEN
           n_acc = 0
           DO try = 1, n_pivot
              CALL pivot ( s, pivot_max, range, q, accepted )
              IF ( accepted ) n_acc = n_acc + 1
              CALL update_histogram
           END DO
           p_ratio = REAL(n_acc) / REAL(n_pivot)
        END IF

        c_ratio = 0.0
        IF ( n_crank > 0 ) THEN
           n_acc = 0
           DO try = 1, n_crank
              CALL crank ( s, crank_max, range, q, accepted )
              IF ( accepted ) n_acc = n_acc + 1
              CALL update_histogram
           END DO
           c_ratio = REAL(n_acc) / REAL(n_crank)
        END IF

        ! Calculate and accumulate variables for this step
        CALL calculate ( )
        CALL blk_add ( variables )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                     ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk         ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r ) ! Save configuration
     CALL write_histogram ( his_prefix//sav_tag )             ! Save histogram

  END DO ! End loop over blocks

  CALL run_end ! Output run averages

  CALL calculate ( 'Final values' )

  ! Final energy calculation and overlap check
  IF ( weight() == 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_chain_nvt_sw'
  END IF
  q = qcount ( range ) ! Count all non-bonded interactions
  CALL calculate ( 'Final check' )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r )
  CALL time_stamp

  CALL deallocate_arrays
  DEALLOCATE ( h, s )
  CALL conclusion

CONTAINS

  SUBROUTINE calculate ( string )
    USE averages_module, ONLY : write_variables, msd
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: r_r, c_r, p_r, e_x, r_g, c_x
    REAL, DIMENSION(3)  :: rcm
    real                :: rsq

    ! Preliminary calculations
    rcm = SUM ( r, dim=2 ) / REAL(n)                                 ! Centre of mass
    rsq = SUM ( ( r - SPREAD(rcm,dim=2,ncopies=n) ) ** 2 ) / REAL(n) ! Mean-squared distance from COM

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Acceptance ratios for regrowth, crankshaft, and pivot moves

    IF ( PRESENT ( string ) ) THEN ! The ratios are meaningless in this case
       r_r = variable_type ( nam = 'Regrow ratio', val = 0.0 )
       c_r = variable_type ( nam = 'Crank ratio',  val = 0.0 )
       p_r = variable_type ( nam = 'Pivot ratio',  val = 0.0 )
    ELSE
       r_r = variable_type ( nam = 'Regrow ratio', val = r_ratio )
       c_r = variable_type ( nam = 'Crank ratio',  val = c_ratio )
       p_r = variable_type ( nam = 'Pivot ratio',  val = p_ratio )
    END IF

    ! Total potential energy of system  (extensive, i.e. not divided by N)
    ! PE=negative of the number of square-well contacts
    e_x = variable_type ( nam = 'PE whole-chain', val = -REAL(q) )

    ! Radius of gyration
    r_g = variable_type ( nam = 'Rg', val = SQRT(rsq) )

    ! Heat Capacity (excess, without ideal gas contribution, extensive)
    ! MSD of total PE / T
    ! PE=negative of the number of square-well contacts
    c_x = variable_type ( nam = 'Cv(ex) whole-chain', val = -REAL(q)/temperature, method = msd )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ r_r, c_r, p_r, e_x, r_g, c_x ]

    IF ( PRESENT ( string ) ) THEN ! Output required
       WRITE ( unit=output_unit, fmt='(a)' ) string
       CALL write_variables ( variables(4:5) ) ! Not acceptance ratios or heat capacity
    END IF

  END SUBROUTINE calculate

  SUBROUTINE update_histogram
    IMPLICIT NONE

    IF ( q > nq .OR. q < 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'q out of range ', q, nq
       STOP 'Error in update_histogram'
    END IF

    h(q)  = h(q) + 1.0
    q_max = MAX ( q, q_max )

  END SUBROUTINE update_histogram

  SUBROUTINE write_histogram ( filename )
    USE, INTRINSIC :: iso_fortran_env, ONLY : iostat_end, iostat_eor
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: filename

    INTEGER :: q, ioerr, his_unit

    OPEN ( newunit=his_unit, file=filename, status='replace', action='write', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)' ) 'Error opening ', filename, ioerr
       STOP 'Error in write_histogram'
    END IF

    DO q = 0, q_max
       WRITE ( unit=his_unit, fmt='(i10,es20.8)') q, h(q)
    END DO

    CLOSE ( unit=his_unit )

  END SUBROUTINE write_histogram

END PROGRAM mc_chain_nvt_sw

