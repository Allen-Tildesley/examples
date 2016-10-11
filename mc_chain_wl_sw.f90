! mc_chain_wl_sw.f90
! Monte Carlo, single chain, Wang-Landau, square wells
PROGRAM mc_chain_wl_sw
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       regrow, crank, pivot, qcount, weight, &
       &                       n, r

  IMPLICIT NONE

  ! Takes in a configuration of atom positions in a linear chain
  ! NO periodic boundary conditions, no box
  ! Conducts Monte Carlo, Wang-Landau method, using various moves
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

  ! The model is the same one studied by
  ! MP Taylor, W Paul, K Binder, J Chem Phys, 131, 114907 (2009)

  ! Most important variables
  REAL    :: bond           ! bond length (in units of core diameter)
  REAL    :: range          ! range of attractive well (specified, same units)
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
  REAL    :: flatness       ! flatness criterion
  LOGICAL :: flat           ! Flag indicating that histogram is flat
  INTEGER :: nstep          ! Steps per block (at end of which we do flatness check)
  INTEGER :: stage          ! Counts stages in entropy function refinement
  INTEGER :: nstage         ! Determines termination point
  REAL    :: ds             ! Entropy increment for Wang-Landau method
  INTEGER :: nq             ! Maximum anticipated energy

  REAL, DIMENSION(:), ALLOCATABLE :: h ! Histogram of q values (0:nq)
  REAL, DIMENSION(:), ALLOCATABLE :: s ! Entropy histogram used in acceptance (0:nq)

  INTEGER :: blk, stp, ioerr, try, n_acc, q_max
  logical :: accepted

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.', his_prefix = 'his.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with stage number

  NAMELIST /nml/ nstep, nstage, m_max, k_max, crank_max, crank_fraction, pivot_max, pivot_fraction, &
       &         range, flatness

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_chain_wl_sw'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, Wang-Landau method, chain molecule, square wells'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing

  nstage         = 20    ! 2**(-20) approx 10**(-6) for smallest modification factor
  nstep          = 10000 ! number of steps per block
  m_max          = 3     ! maximum atoms in regrow
  k_max          = 32    ! number of random tries per atom in regrow
  crank_max      = 0.5   ! maximum move angle in crankshaft
  crank_fraction = 0.5   ! fraction of atoms to try in crank moves
  pivot_max      = 0.5   ! maximum move angle in pivot
  pivot_fraction = 0.2   ! fraction of atoms to try in pivot moves
  range          = 1.3   ! range of attractive well
  flatness       = 0.8   ! histogram flatness criterion

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_chain_wl_sw'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of stages of refinement',  nstage
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',       nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Flatness criterion',              flatness
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Max atoms in regrow',             m_max
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Random tries per atom in regrow', k_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Max move angle in crankshaft',    crank_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Crank fraction',                  crank_fraction
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Max move angle in pivot',         pivot_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Pivot fraction',                  pivot_fraction
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Attractive well range',           range

  IF ( range < 1.0 ) THEN
     WRITE ( unit=output_unit, fmt='(a)' ) 'Warning, range < core diameter (1.0)'
  END IF

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond ) ! First call is just to get n and bond
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',          n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Bond length (in sigma units)', bond
  n_crank = NINT(crank_fraction*REAL(n))
  n_pivot = NINT(pivot_fraction*REAL(n))
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of atoms in crankshaft',  n_crank
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of atoms in pivot',       n_pivot

  CALL allocate_arrays
  nq = 6*n ! Anticipated maximum number of pair interactions within range
  ALLOCATE ( h(0:nq), s(0:nq) )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r ) ! Second call gets r

  IF ( weight() == 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_chain_wl_sw'
  END IF
  q = qcount ( range )
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Initial energy', q
  q_max = q ! Max q seen so far

  CALL run_begin ( [ CHARACTER(len=15) :: 'Regrow ratio', 'Crank ratio', 'Pivot ratio' ] )

  stage = 1
  blk   = 0
  ds    = 1.0 ! Initial entropy refinement term
  flat  = .FALSE.
  s(:)  = 0.0
  h(:)  = 0.0

  DO ! Begin loop over blocks

     IF ( stage > nstage ) EXIT ! Run is finished

     IF  ( flat ) THEN ! Start of next stage
        h(:) = 0.0
        s(1:q_max) = s(1:q_max) - MINVAL ( s(1:q_max) ) ! Reset baseline for entropy
        flat = .FALSE.
     END IF

     blk = blk + 1
     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL regrow ( s, m_max, k_max, bond, range, q, accepted )
        IF ( accepted ) THEN
           regrow_ratio = 1.0
        ELSE
           regrow_ratio = 0.0
        END IF
        CALL update_histogram

        n_acc = 0
        DO try = 1, n_pivot
           CALL pivot ( s, pivot_max, range, q, accepted )
           IF ( accepted ) n_acc = n_acc + 1
           CALL update_histogram
        END DO
        IF ( n_pivot > 0 ) THEN
           pivot_ratio = REAL(n_acc) / real(n_pivot)
        ELSE
           pivot_ratio = 0.0
        END IF

        n_acc = 0
        DO try = 1, n_crank
           CALL crank ( s, crank_max, range, q, accepted )
           IF ( accepted ) n_acc = n_acc + 1
           CALL update_histogram
        END DO
        IF ( n_crank > 0 ) THEN
           crank_ratio = REAL(n_acc) / REAL(n_crank)
        ELSE
           crank_ratio = 0.0
        END IF

        ! Calculate all variables for this step
        CALL blk_add ( [regrow_ratio,crank_ratio,pivot_ratio] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )

     flat = histogram_flat ( flatness ) ! Check for flatness

     IF ( flat ) THEN ! End of this stage
        IF ( nstage < 1000 ) WRITE(sav_tag,'(i3.3)') stage       ! number configuration by stage
        CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r ) ! save configuration
        CALL write_histogram ( his_prefix//sav_tag )             ! save histogram
        WRITE ( unit=output_unit, fmt='(a,i2,a,2i5)' ) 'End of stage ', stage, ' q diagnostics ', q_max, COUNT(h(:)>0.5)
        stage = stage + 1 ! Ready for next stage
        ds    = ds / 2.0  ! Entropy change reduction
     END IF

  END DO ! End loop over blocks

  CALL deallocate_arrays
  DEALLOCATE ( h, s )
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE update_histogram

    ! This routine simply updates the probability histogram of (negative) energies
    ! It also updates the entropy histogram by ds

    IF ( q > nq .OR. q < 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'q out of range ', q, nq
       STOP 'Error in update_histogram'
    END IF

    h(q) = h(q) + 1.0
    s(q) = s(q) + ds
    q_max = MAX ( q, q_max )

  END SUBROUTINE update_histogram

  FUNCTION histogram_flat ( flatness ) RESULT ( flat )
    LOGICAL            :: flat     ! Returns a signal that the histogram is "sufficiently flat"
    REAL,   INTENT(in) :: flatness ! Specified degree of flatness to pass the test

    ! We only look at the histogram up to the maximum q visited during the run
    ! The flatness parameter should be in (0,1), higher values corresponding to flatter histograms

    REAL :: avg

    IF ( flatness <= 0.0 .OR. flatness >= 1.0 ) THEN ! Ensure sensible value
       WRITE ( unit=error_unit, fmt='(a,f15.5)' ) 'Flatness error ', flatness
       STOP 'Error in histogram_flat'
    END IF

    avg  = SUM ( h(1:q_max) ) / REAL(q_max)

    IF ( avg < 0.0 ) THEN ! This should never happen, unless h somehow overflows
       WRITE ( unit=error_unit, fmt='(a,*(es20.8))' ) 'Error in h ', h(1:q_max)
       STOP 'Error in histogram_flat'
    END IF

    flat = REAL(MINVAL(h(1:q_max))) > flatness*avg

  END FUNCTION histogram_flat

  SUBROUTINE write_histogram ( filename )
    USE, INTRINSIC :: iso_fortran_env, ONLY : iostat_end, iostat_eor

    CHARACTER(len=*), INTENT(in) :: filename

    INTEGER :: q, ioerr, his_unit

    OPEN ( newunit=his_unit, file=filename, status='replace', action='write', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)' ) 'Error opening ', filename, ioerr
       STOP 'Error in write_histogram'
    END IF

    DO q = 0, q_max
       WRITE ( unit=his_unit, fmt='(i10,2es20.8)') q, h(q), s(q)
    END DO

    CLOSE ( unit=his_unit )

  END SUBROUTINE write_histogram

END PROGRAM mc_chain_wl_sw

