! mc_chain_wl_sw.f90
! Monte Carlo, single chain, Wang-Landau, square wells
PROGRAM mc_chain_wl_sw

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

  ! In this program, the flatness of the h(q) histogram is checked at the end of every block
  ! A stage consists of many blocks, during which the entropy is adjusted by a value ds
  ! The stage is deemed to be finished when h(q) is sufficiently flat
  ! At the end of each stage, the histograms are printed and reset, and ds is divided by 2

  ! The model is the same one studied by
  ! MP Taylor, W Paul, K Binder, J Chem Phys, 131, 114907 (2009)

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     regrow, crank, pivot, qcount, weight, n, r

  IMPLICIT NONE

  ! Most important variables
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
  INTEGER :: q_min          ! min value of q seen in simulation (maximum energy)
  INTEGER :: q_max          ! max value of q seen in simulation (minimum energy)
  REAL    :: flatness       ! flatness criterion
  LOGICAL :: flat           ! Flag indicating that histogram is flat
  INTEGER :: nstep          ! Steps per block (at end of which we do flatness check)
  INTEGER :: stage          ! Counts stages in entropy function refinement
  INTEGER :: nstage         ! Determines termination point for run
  REAL    :: ds             ! Entropy increment for Wang-Landau method
  INTEGER :: nq             ! Maximum anticipated q

  ! Histograms and tables
  REAL, DIMENSION(:), ALLOCATABLE :: h ! Histogram of q values (0:nq)
  REAL, DIMENSION(:), ALLOCATABLE :: g ! Histogram of radius of gyration (0:nq)
  REAL, DIMENSION(:), ALLOCATABLE :: s ! Entropy histogram used in acceptance (0:nq)
  REAL, DIMENSION(:), ALLOCATABLE :: e ! Energy histogram

  INTEGER :: blk, stp, ioerr, try, n_acc
  LOGICAL :: accepted
  REAL    :: r_ratio, c_ratio, p_ratio

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=4), PARAMETER :: his_prefix = 'his.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with stage number

  NAMELIST /nml/ nstep, nstage, m_max, k_max, crank_max, crank_fraction, pivot_max, pivot_fraction, &
       &         range, flatness

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_chain_wl_sw'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, Wang-Landau method, chain molecule, square wells'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nstage         = 20    ! 2**(-20) = approx 10**(-6) for smallest modification factor
  nstep          = 10000 ! number of steps per block
  m_max          = 3     ! maximum atoms in regrow
  k_max          = 32    ! number of random tries per atom in regrow
  crank_max      = 0.5   ! maximum move angle in crankshaft
  crank_fraction = 0.5   ! fraction of atoms to try in crank moves
  pivot_max      = 0.5   ! maximum move angle in pivot
  pivot_fraction = 0.2   ! fraction of atoms to try in pivot moves
  range          = 1.5   ! range of attractive well
  flatness       = 0.9   ! histogram flatness criterion

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_chain_wl_sw'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of stages of refinement',  nstage
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',       nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Flatness criterion',              flatness
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Max atoms in regrow',             m_max
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Random tries per atom in regrow', k_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Max move angle in crankshaft',    crank_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Crank fraction',                  crank_fraction
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Max move angle in pivot',         pivot_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Pivot fraction',                  pivot_fraction
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
  ALLOCATE ( h(0:nq), g(0:nq), s(0:nq), e(0:nq) )
  e = [ ( -REAL(q), q = 0, nq ) ] ! Energy values
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r ) ! Second call gets r

  ! Set number of crankshaft and pivot moves per step
  n_crank = NINT(crank_fraction*REAL(n))
  n_pivot = NINT(pivot_fraction*REAL(n))
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of crankshaft tries per step', n_crank
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of pivot tries per step',      n_pivot

  ! Initial energy calculation plus overlap check
  IF ( weight() == 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_chain_wl_sw'
  END IF
  q = qcount ( range )
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Initial energy', q
  q_min = q ! Min q seen so far
  q_max = q ! Max q seen so far

  ! Initialize arrays for averaging and write column headings
  r_ratio = 0.0
  c_ratio = 0.0
  p_ratio = 0.0
  h(:) = 0.0
  CALL run_begin ( calc_variables() )

  stage = 1
  blk   = 0
  ds    = 1.0 ! Initial entropy refinement term
  flat  = .FALSE.
  s(:)  = 0.0
  h(:)  = 0.0
  g(:)  = 0.0

  DO ! Begin loop over blocks

     IF ( stage > nstage ) EXIT ! Run is finished

     IF  ( flat ) THEN ! Start of next stage
        h(:) = 0.0
        g(:) = 0.0
        flat = .FALSE.
     END IF

     blk = blk + 1
     CALL blk_begin

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
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk ) ! Output block averages

     flat = histogram_flat ( flatness ) ! Check for flatness

     IF ( flat ) THEN ! End of this stage
        WRITE ( unit=output_unit, fmt='(a)' ) REPEAT('-',63)
        WRITE ( unit=output_unit, fmt='(a,i3,a,3i3)' ) 'stage', stage, &
             & ' q_min q_max q_count =', q_min, q_max, COUNT(h(q_min:q_max)>0.5)
        WRITE ( unit=output_unit, fmt='(a)' ) REPEAT('-',63)
        IF ( nstage < 1000 ) WRITE(sav_tag,'(i3.3)') stage       ! Number configuration by stage
        CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r ) ! Save configuration
        CALL write_histogram ( his_prefix//sav_tag )             ! Save histogram
        stage = stage + 1 ! Ready for next stage
        ds    = ds / 2.0  ! Entropy change reduction
     END IF

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages, although these are really not relevant
  CALL write_results                ! Specimen results at selected temperatures

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r ) ! Write out final configuration

  CALL deallocate_arrays
  DEALLOCATE ( h, g, s, e )
  CALL conclusion

CONTAINS

  FUNCTION calc_variables () RESULT ( variables )
    USE averages_module, ONLY : variable_type
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(3) :: variables ! The 3 variables listed below

    ! This function returns all variables of interest in an array, for use in the main program

    TYPE(variable_type) :: r_r, c_r, p_r

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Acceptance ratios for regrowth, crankshaft, and pivot moves
    r_r = variable_type ( nam = 'Regrow ratio', val = r_ratio, instant = .FALSE. )
    c_r = variable_type ( nam = 'Crank ratio',  val = c_ratio, instant = .FALSE. )
    p_r = variable_type ( nam = 'Pivot ratio',  val = p_ratio, instant = .FALSE. )

    ! Collect together for averaging
    variables = [ r_r, c_r, p_r ]

  END FUNCTION calc_variables

  SUBROUTINE update_histogram
    IMPLICIT NONE

    ! This routine updates the probability histogram of (negative) energies
    ! It also updates the entropy histogram by ds
    ! and a histogram of values of radius of gyration

    REAL, DIMENSION(3) :: r_cm
    REAL               :: r_g

    IF ( q > nq .OR. q < 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'q out of range ', q, nq
       STOP 'Error in update_histogram'
    END IF

    r_cm = SUM ( r, dim=2 ) / REAL(n) ! Centre of mass
    r_g  = SQRT ( SUM ( ( r - SPREAD(r_cm,dim=2,ncopies=n) ) ** 2 ) / REAL(n) )

    h(q) = h(q) + 1.0
    s(q) = s(q) + ds
    g(q) = g(q) + r_g

    q_min = MIN ( q, q_min )
    q_max = MAX ( q, q_max )

  END SUBROUTINE update_histogram

  FUNCTION histogram_flat ( flatness ) RESULT ( flat )
    IMPLICIT NONE
    LOGICAL            :: flat     ! Returns a signal that the histogram is "sufficiently flat"
    REAL,   INTENT(in) :: flatness ! Specified degree of flatness to pass the test

    ! We only look at the histogram in the range of q visited during the run
    ! The flatness parameter should be in (0,1), higher values corresponding to flatter histograms

    REAL :: avg

    IF ( flatness <= 0.0 .OR. flatness >= 1.0 ) THEN ! Ensure sensible value
       WRITE ( unit=error_unit, fmt='(a,f15.6)' ) 'Flatness error ', flatness
       STOP 'Error in histogram_flat'
    END IF

    avg = SUM ( h(q_min:q_max) ) / REAL(q_max+1-q_min)

    IF ( avg < 0.0 ) THEN ! This should never happen, unless h somehow overflows
       WRITE ( unit=error_unit, fmt='(a,*(es20.8))' ) 'Error in h ', h(q_min:q_max)
       STOP 'Error in histogram_flat'
    END IF

    flat = REAL(MINVAL(h(q_min:q_max))) > flatness*avg

  END FUNCTION histogram_flat

  SUBROUTINE write_histogram ( filename )
    USE, INTRINSIC :: iso_fortran_env, ONLY : iostat_end, iostat_eor
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: filename

    ! This routine writes out the histograms at the end of each stage
    ! Note that h and g will be reset to zero at the start of the next stage
    ! so it is OK to normalize them here
    ! Also we reset the baseline for entropy to avoid the numbers getting too large
    ! We set s(q_min)=0.0 for convenience; relative to this, most other entropies
    ! will be negative, but this is not a problem.

    INTEGER :: q, ioerr, his_unit
    REAL    :: norm

    OPEN ( newunit=his_unit, file=filename, status='replace', action='write', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)' ) 'Error opening ', filename, ioerr
       STOP 'Error in write_histogram'
    END IF

    ! Normalize radius of gyration entries
    WHERE ( h > 0.5 ) g = g / h

    ! Normalize h, converting it to a set of probabilities
    norm = SUM(h)
    h    = h / norm

    ! Reset baseline for entropy
    norm           = s(q_min)
    s(q_min:q_max) = s(q_min:q_max) - norm

    DO q = q_min, q_max
       WRITE ( unit=his_unit, fmt='(3f15.6,es20.8)') e(q), h(q), g(q), s(q)
    END DO

    CLOSE ( unit=his_unit )

  END SUBROUTINE write_histogram

  SUBROUTINE write_results

    IMPLICIT NONE

    ! Calculates some specimen results after the simulation is finished
    ! The same calculation can be done afterwards using the data written out by write_histogram
    ! An example program wl_hist.f90 is provided to illustrate this

    REAL, DIMENSION(:), ALLOCATABLE :: t_vals ! Array of temperatures
    REAL, DIMENSION(:), ALLOCATABLE :: boltz  ! Array of Boltzmann factors (including DOS)

    INTEGER :: i, qs_max
    REAL    :: t, s_max, norm, g_avg, e_avg, e_msd

    ALLOCATE ( boltz(q_min:q_max) )
    ALLOCATE ( t_vals(10) ) ! Ideally we could omit this as Fortran 2003 should automatically allocate
    t_vals = [ 0.15, 0.18, 0.2, 0.22, 0.25, 0.3, 0.5, 1.0, 2.0, 5.0 ]

    WRITE ( unit=output_unit, fmt='(a)'   ) 'Specimen results'
    WRITE ( unit=output_unit, fmt='(4a15)') 'T', 'Rg', 'PE/N', 'Cv(ex)/N'

    DO i = 1, SIZE(t_vals) ! Loop over selected temperatures
       t = t_vals(i)

       ! Locate maximum Boltzmann factor (helps avoid overflows)
       qs_max = MAXLOC ( s(q_min:q_max) - e(q_min:q_max) / t, dim=1 )
       s_max  = s(qs_max + q_min - 1)

       ! Compute Boltzmann weights including density of states
       boltz = s(q_min:q_max) - s_max - e(q_min:q_max) / t
       boltz = EXP ( boltz )

       ! Calculate averages
       norm  = SUM ( boltz )
       g_avg = SUM ( boltz * g(q_min:q_max) ) / norm
       e_avg = SUM ( boltz * e(q_min:q_max) ) / norm
       e_msd = SUM ( boltz * (e(q_min:q_max) - e_avg )**2  ) / norm
       e_avg = e_avg / REAL(n)            ! Energy per atom
       e_msd = e_msd / ( REAL(n) * t**2 ) ! Heat capacity per atom

       WRITE ( unit=output_unit, fmt='(4f15.6)' ) t, g_avg, e_avg, e_msd

    END DO ! End loop over selected temperatures

    DEALLOCATE ( boltz )

  END SUBROUTINE write_results

END PROGRAM mc_chain_wl_sw

