! mc_chain_nvt_sw.f90
! Monte Carlo, single chain, NVT, square wells
PROGRAM mc_chain_nvt_sw
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_module,        ONLY : model_description, allocate_arrays, deallocate_arrays, &
       &                       regrow, cranks, pivots, write_histogram, verbose, qcount, weight, &
       &                       n, nq, range, bond, r, h, s, wl

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
  ! Energy is -q, where q is the total number of attractive well interactions
  ! This is a negative integer, but for convenience we refer to q as energy
  ! Configurational weights are calculated on the basis of the athermal interactions
  ! only, i.e. weight=1 if no overlap, weight=0 if any overlap

  ! Most important variables
  REAL    :: temperature    ! temperature (specified, in units of well depth)
  REAL    :: regrow_ratio   ! acceptance ratio for regrowth moves
  REAL    :: crank_ratio    ! acceptance ratio for crankshaft moves
  REAL    :: pivot_ratio    ! acceptance ratio for pivot moves
  INTEGER :: m_max          ! maximum atoms in regrow
  INTEGER :: k_max          ! number of random tries per atom in regrow
  REAL    :: crank_max      ! maximum move angle in crankshaft
  REAL    :: crank_fraction ! fraction of atoms to try in crankshaft
  INTEGER :: n_crank        ! number of atoms to try in crankshaft
  REAL    :: pivot_max      ! maximum move angle in pivot
  REAL    :: pivot_fraction ! fraction of atoms to try in pivot
  INTEGER :: n_pivot        ! number of atoms to try in pivot
  INTEGER :: q              ! total attractive interaction (negative of energy)

  INTEGER :: blk, stp, nstep, nblock, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.', his_prefix = 'his.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, m_max, k_max, crank_max, crank_fraction, pivot_max, pivot_fraction, &
       &         temperature, range

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_chain_nvt_sw'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble, chain molecule, square wells'
  CALL model_description ( output_unit )
  CALL time_stamp ( output_unit )
  wl = .FALSE. ! not using the Wang-Landau method

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing

  nblock         = 50    ! number of blocks
  nstep          = 10000 ! number of sweeps per block
  m_max          = 3     ! maximum atoms in regrow
  k_max          = 32    ! number of random tries per atom in regrow
  crank_max      = 0.5   ! maximum move angle in crankshaft
  crank_fraction = 0.5   ! fraction of atoms to try in crank moves
  pivot_max      = 0.5   ! maximum move angle in pivot
  pivot_fraction = 0.2   ! fraction of atoms to try in pivot moves
  temperature    = 1.0   ! temperature (in units of well depth)
  range          = 0.5   ! range of attractive well

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_chain_nvt_sw'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',                nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',       nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Max atoms in regrow',             m_max
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Random tries per atom in regrow', k_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Max move angle in crankshaft',    crank_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Crank fraction',                  crank_fraction
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Max move angle in pivot',         pivot_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Pivot fraction',                  pivot_fraction
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature/well depth',          temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Attractive well range',           range

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond ) ! First call is just to get n and bond
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',          n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Bond length (in sigma units)', bond
  n_crank = NINT(crank_fraction*REAL(n))
  n_pivot = NINT(pivot_fraction*REAL(n))
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of atoms in crankshaft',  n_crank
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of atoms in pivot',       n_pivot

  CALL allocate_arrays

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r ) ! Second call gets r

  ! initialize Boltzmann exponents, s(q), which are just values of E/kT
  s = [ (-REAL(q)/temperature, q = 0, nq) ]

  verbose = .TRUE.
  IF ( weight() == 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_chain_nvt_sw'
  END IF
  q = qcount()
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Initial energy', q
  verbose = .FALSE.

  CALL run_begin ( [ CHARACTER(len=15) :: 'Regrow ratio', 'Crank ratio', 'Pivot ratio', 'Energy' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin
     h = 0.0

     DO stp = 1, nstep ! Begin loop over steps

        CALL regrow ( m_max, k_max, q, regrow_ratio )
        CALL pivots ( n_pivot, pivot_max, q, pivot_ratio )
        CALL cranks ( n_crank, crank_max, q, crank_ratio )

        ! Calculate all variables for this step
        CALL blk_add ( [regrow_ratio,crank_ratio,pivot_ratio,REAL(q)] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk         ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r ) ! save configuration
     CALL write_histogram ( his_prefix//sav_tag )             ! save histogram

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Final pair interactions', q

  verbose = .TRUE.
  IF ( weight() == 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_chain_nvt_sw'
  END IF
  q = qcount() ! count all non-bonded interactions
  verbose = .FALSE.
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Check final pair interactions', q

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays

END PROGRAM mc_chain_nvt_sw

