! qmc_pi_lj.f90
! Quantum Monte Carlo, path-integral method
PROGRAM qmc_pi_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE maths_module,     ONLY : metropolis, random_translate_vector
  USE qmc_module,       ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       potential_1, spring_1, potential, spring, n, p, r, potential_type

  IMPLICIT NONE

  ! Takes in a set of configurations of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts path-integral Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1
  ! The importance of quantum effects is specified through the reduced de Boer length
  ! lambda = (hbar/sqrt(mass*epsilon))/sigma which takes typical values of
  ! 0.01 for Xe, 0.03 for Ar, and 0.095 for Ne.
  ! This means that the quantum spring constant may be expressed k_spring = P*(T/lambda)**2
  ! where T stands for the reduced temperature kB*T/epsilon

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in qmc_module

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: lambda      ! de Boer length
  REAL :: dr_max      ! Maximum MC displacement
  REAL :: temperature ! Specified temperature
  REAL :: r_cut       ! Potential cutoff distance

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  ! Composite interaction = pot & ovr variables
  TYPE(potential_type) :: total, partial_old, partial_new

  INTEGER            :: blk, stp, i, k, nstep, nblock, moves, ioerr
  REAL               :: total_spr, partial_old_spr, partial_new_spr, delta, k_spring, m_ratio
  REAL, DIMENSION(3) :: rik

  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav'    ! May be overwritten with block number
  CHARACTER(len=6)            :: cnf_prefix = 'cnf##.' ! Will have polymer id inserted
  CHARACTER(len=2)            :: k_tag                 ! Will contain polymer id

  NAMELIST /nml/ nblock, nstep, p, temperature, r_cut, dr_max, lambda

  WRITE ( unit=output_unit, fmt='(a)' ) 'qmc_pi_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Path-integral Monte Carlo, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Simulation uses cut (but not shifted) potential'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  r_cut       = 2.5
  dr_max      = 0.15
  p           = 4
  lambda      = 0.1

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in qmc_pi_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Size of each ring polymer', p
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum displacement',      dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'de Boer length',            lambda
  IF ( p < 2 .OR. p > 99 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'p must lie between 2 and 99', p
     STOP 'Error in qmc_pi_lj'
  END IF
  k_spring = REAL(p) * ( temperature / lambda ) ** 2
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Quantum spring constant', k_spring

  ! Read in initial configuration and allocate necessary arrays
  ! Read each polymer index from a unique file; we assume that n and box are all the same!
  cnf_prefix(4:5) = '01' ! Read first LJ configuration to get basic parameters n and box
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box )
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n) / box**3
  CALL allocate_arrays ( box, r_cut ) ! Allocate r
  DO k = 1, p ! Loop over ring polymer indices
     WRITE(k_tag,fmt='(i2.2)') k ! Convert into character form
     cnf_prefix(4:5) = k_tag     ! Insert into configuration filename
     CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r(:,:,k) ) ! Read r for polymer index k
  END DO ! End loop over ring polymer indices
  r(:,:,:) = r(:,:,:) / box                ! Convert all positions to box units
  r(:,:,:) = r(:,:,:) - ANINT ( r(:,:,:) ) ! Periodic boundaries (box=1 units)

  ! Calculate classical LJ and quantum spring potential energies & check overlap
  total = potential ( box, r_cut )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in qmc_pi_lj'
  END IF
  total_spr = spring ( box, k_spring )
  CALL calculate ( 'Initial values' )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( variables )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           DO k = 1, p ! Begin loop over ring polymer indices

              partial_old = potential_1 ( r(:,i,k), i, k, box, r_cut ) ! Old atom classical potential etc

              IF ( partial_old%ovr ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in qmc_pi_lj'
              END IF

              partial_old_spr = spring_1 ( r(:,i,k), i, k, box, k_spring ) ! Old atom quantum potential

              rik(:) = random_translate_vector ( dr_max/box, r(:,i,k) ) ! Trial move to new position (in box=1 units) 
              rik(:) = rik(:) - ANINT ( rik(:) )                        ! Periodic boundary correction

              partial_new = potential_1 ( rik, i, k, box, r_cut ) ! New atom classical potential etc

              IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

                 partial_new_spr = spring_1 ( rik, i, k, box, k_spring ) ! New atom quantum potential

                 delta =         partial_new%pot - partial_old%pot ! Change in classical cut (but not shifted) potential
                 delta = delta + partial_new_spr - partial_old_spr ! Add change in quantum potential
                 delta = delta / temperature                       ! Divide by temperature

                 IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                    total     = total     + partial_new     - partial_old     ! Update classical total values
                    total_spr = total_spr + partial_new_spr - partial_old_spr ! Update quantum system potential
                    r(:,i,k)  = rik                                           ! Update position
                    moves     = moves + 1                                     ! Increment move counter
                 END IF ! End accept Metropolis test

              END IF ! End test for non-overlapping configuration

           END DO ! End loop over ring polymer indices

        END DO ! End loop over atoms

        m_ratio = REAL(moves) / REAL(n)

        ! Calculate and accumulate variables for this step
        CALL calculate ( )
        CALL blk_add ( variables )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                                  ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk                  ! Number configuration by block
     DO k = 1, p ! Loop over ring polymer indices
        WRITE(k_tag,fmt='(i2.2)') k                                        ! Convert into character form
        cnf_prefix(4:5) = k_tag                                            ! Insert into configuration filename
        CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r(:,:,k)*box ) ! Save to unique file
     END DO ! End loop over ring polymer indices

  END DO ! End loop over blocks

  CALL run_end ! Output run averages

  CALL calculate ( 'Final values' )

  ! Final double-check on book-keeping for totals, and overlap
  total = potential ( box, r_cut )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in qmc_pi_lj'
  END IF
  total_spr = spring ( box, k_spring )
  CALL calculate ( 'Final check' )

  DO k = 1, p ! Loop over ring polymer indices
     WRITE(k_tag,fmt='(i2.2)') k                                        ! Convert into character form
     cnf_prefix(4:5) = k_tag                                            ! Insert into configuration filename
     CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r(:,:,k)*box ) ! Write to unique file
  END DO ! End loop over ring polymer indices

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE calculate ( string )
    USE lrc_module,      ONLY : potential_lrc
    USE averages_module, ONLY : write_variables, variable_type
    IMPLICIT NONE
    CHARACTER (len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    ! In this example we simulate using the cut (but not shifted) potential
    ! The values of < e_c > and < density > should be consistent (for this potential)
    ! For comparison, long-range corrections are also applied to give
    ! an estimate of < e_f > for the full (uncut) potential
    ! The value of the cut-and-shifted potential is not used, in this example

    TYPE(variable_type) :: m_r, e_c, e_f
    REAL                :: vol, rho, kin

    ! Preliminary calculations
    vol = box**3                    ! Volume
    rho = REAL(n) / vol             ! Density
    kin = 1.5 * n * p * temperature ! Average kinetic energy for NP-atom system

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Acceptance ratio of moves

    IF ( PRESENT ( string ) ) THEN ! The ratio is meaningless in this case
       m_r = variable_type ( nam = 'Move ratio', val = 0.0 )
    ELSE
       m_r = variable_type ( nam = 'Move ratio', val = m_ratio )
    END IF

    ! Internal energy per atom for simulated, cut, potential
    ! Total (cut but not shifted) PE already divided by factor P
    ! plus total classical KE for NP-atom system MINUS total spring potential
    ! all divided by N
    e_c = variable_type ( nam = 'E/N cut', val = (kin+total%pot-total_spr)/REAL(n) )

    ! Internal energy per atom for full potential with LRC
    ! LRC plus total (cut but not shifted) PE already divided by factor P
    ! plus total classical KE for NP-atom system MINUS total spring potential
    ! all divided by N
    e_f = variable_type ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + (kin+total%pot-total_spr)/REAL(n) )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ m_r, e_c, e_f ]

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       CALL write_variables ( variables(2:) ) ! Don't write out move ratio
    END IF

  END SUBROUTINE calculate

END PROGRAM qmc_pi_lj

