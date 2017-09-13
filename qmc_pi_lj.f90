! qmc_pi_lj.f90
! Quantum Monte Carlo, path-integral method
PROGRAM qmc_pi_lj

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

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               maths_module,     ONLY : metropolis, random_translate_vector
  USE               qmc_module,       ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     potential_1, spring_1, potential, spring, n, p, r, potential_type

  IMPLICIT NONE

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: lambda      ! de Boer length
  REAL :: dr_max      ! Maximum MC displacement
  REAL :: dc_max      ! Maximum centre-of-mass MC displacement
  REAL :: temperature ! Specified temperature
  REAL :: r_cut       ! Potential cutoff distance

  ! Composite interaction = pot & ovr variables
  TYPE(potential_type) :: total, partial_old, partial_new

  INTEGER            :: blk, stp, i, k, km, kp, nstep, nblock, r_moves, c_moves, ioerr
  REAL               :: total_spr, partial_old_spr, partial_new_spr, delta, k_spring
  REAL               :: r_ratio, c_ratio
  REAL               :: rad_g
  REAL, DIMENSION(3) :: rik, dc

  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav'    ! May be overwritten with block number
  CHARACTER(len=6)            :: cnf_prefix = 'cnf##.' ! Will have polymer id inserted
  CHARACTER(len=2)            :: k_tag                 ! Will contain polymer id

  NAMELIST /nml/ nblock, nstep, p, temperature, r_cut, dr_max, dc_max, lambda

  WRITE ( unit=output_unit, fmt='(a)' ) 'qmc_pi_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Path-integral Monte Carlo, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Simulation uses cut (but not shifted) potential'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 10000
  temperature = 0.701087
  r_cut       = 2.5
  dr_max      = 0.05
  dc_max      = 0.1
  p           = 4
  lambda      = 0.092 ! from Ne parameters epsilon=36.8K, sigma=0.2789nm, m=20.18u

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
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum COM displacement',  dc_max
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

  ! Initialize arrays for averaging and write column headings
  r_ratio = 0.0
  c_ratio = 0.0
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        r_moves = 0
        c_moves = 0

        DO i = 1, n ! Begin loop over atoms

           ! Centre of mass move

           dc = 0.0
           dc = random_translate_vector ( dc_max/box, dc ) ! Random displacement
           partial_old = potential_type ( 0.0, 0.0, .FALSE. )
           partial_new = potential_type ( 0.0, 0.0, .FALSE. )
           DO k = 1, p ! Loop over ring polymer indices
              partial_old = partial_old + potential_1 ( r(:,i,k), i, k, box, r_cut ) ! Old atom classical potential etc
              rik(:) = r(:,i,k) + dc
              rik(:) = rik(:) - ANINT ( rik(:) )                                     ! Periodic boundary correction
              partial_new = partial_new + potential_1 ( rik, i, k, box, r_cut )      ! New atom classical potential etc
           END DO

           IF ( partial_old%ovr ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in qmc_pi_lj'
           END IF

           IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

              delta = partial_new%pot - partial_old%pot ! Change in classical cut (but not shifted) potential
              delta = delta / temperature               ! Divide by temperature

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 total     = total + partial_new - partial_old          ! Update classical total values
                 r(:,i,:)  = r(:,i,:) + SPREAD ( dc, dim=2, ncopies=p ) ! Update position
                 r(:,i,:)  = r(:,i,:) - ANINT ( r(:,i,:) )              ! Periodic boundary corrections
                 c_moves   = c_moves + 1                                ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End test for non-overlapping configuration

           ! Individual atom moves

           DO k = 1, p ! Begin loop over ring polymer indices

              partial_old = potential_1 ( r(:,i,k), i, k, box, r_cut ) ! Old atom classical potential etc

              IF ( partial_old%ovr ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in qmc_pi_lj'
              END IF

              ! Get neighbour indices
              kp = k+1
              IF ( kp > p ) kp = 1
              km = k-1
              if ( km < 1 ) km = p
              partial_old_spr = spring_1 ( r(:,i,k), r(:,i,km), box, k_spring ) &
                   &          + spring_1 ( r(:,i,k), r(:,i,kp), box, k_spring ) ! Old atom quantum potential

              rik(:) = random_translate_vector ( dr_max/box, r(:,i,k) ) ! Trial move to new position (in box=1 units) 
              rik(:) = rik(:) - ANINT ( rik(:) )                        ! Periodic boundary correction

              partial_new = potential_1 ( rik, i, k, box, r_cut ) ! New atom classical potential etc

              IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

                 partial_new_spr = spring_1 ( rik, r(:,i,km), box, k_spring ) &
                      &          + spring_1 ( rik, r(:,i,kp), box, k_spring ) ! New atom quantum potential

                 delta =         partial_new%pot - partial_old%pot ! Change in classical cut (but not shifted) potential
                 delta = delta + partial_new_spr - partial_old_spr ! Add change in quantum potential
                 delta = delta / temperature                       ! Divide by temperature

                 IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                    total     = total     + partial_new     - partial_old     ! Update classical total values
                    total_spr = total_spr + partial_new_spr - partial_old_spr ! Update quantum system potential
                    r(:,i,k)  = rik                                           ! Update position
                    r_moves   = r_moves + 1                                   ! Increment move counter
                 END IF ! End accept Metropolis test

              END IF ! End test for non-overlapping configuration

           END DO ! End loop over ring polymer indices

        END DO ! End loop over atoms

        r_ratio = REAL(r_moves) / REAL(n*p)
        c_ratio = REAL(c_moves) / REAL(n)

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                                  ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk                  ! Number configuration by block
     DO k = 1, p ! Loop over ring polymer indices
        WRITE(k_tag,fmt='(i2.2)') k                                        ! Convert into character form
        cnf_prefix(4:5) = k_tag                                            ! Insert into configuration filename
        CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r(:,:,k)*box ) ! Save to unique file
     END DO ! End loop over ring polymer indices

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  ! Final double-check on book-keeping for totals, and overlap
  total = potential ( box, r_cut )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in qmc_pi_lj'
  END IF
  total_spr = spring ( box, k_spring )

  DO k = 1, p ! Loop over ring polymer indices
     WRITE(k_tag,fmt='(i2.2)') k                                        ! Convert into character form
     cnf_prefix(4:5) = k_tag                                            ! Insert into configuration filename
     CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r(:,:,k)*box ) ! Write to unique file
  END DO ! End loop over ring polymer indices

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE lrc_module,      ONLY : potential_lrc, pressure_lrc
    USE averages_module, ONLY : variable_type
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(7) :: variables ! The 7 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    ! In this example we simulate using the cut (but not shifted) potential
    ! but we only report results which have had the long-range corrections applied
    ! The value of the cut-and-shifted potential is not used, in this example

    TYPE(variable_type) :: r_r, c_r, e_f, p_f, e_q, k_q, r_g
    REAL                :: vol, rho, kin, kin_q

    ! Preliminary calculations
    vol   = box**3                    ! Volume
    rho   = REAL(n) / vol             ! Density
    kin   = 1.5 * n * p * temperature ! Average kinetic energy for NP-atom system
    kin_q = kin - total_spr           ! Quantum estimator for kinetic energy
    rad_g = rad_gyr ( r )

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Acceptance ratio of atomic moves
    r_r = variable_type ( nam = 'Atomic move ratio', val = r_ratio, instant = .FALSE. )

    ! Acceptance ratio of centre-of-mass moves
    c_r = variable_type ( nam = 'COM move ratio', val = c_ratio, instant = .FALSE. )

    ! Internal energy per atom for full potential with LRC
    ! LRC plus total (cut but not shifted) PE already divided by factor P
    ! plus KE estimator: total classical KE for NP-atom system MINUS total spring potential
    ! all divided by N
    e_f = variable_type ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + (kin_q+total%pot)/REAL(n) )

    ! Kinetic energy per atom, just for interest's sake
    k_q = variable_type ( nam = 'KE/N', val = kin_q/REAL(n) )

    ! Pressure for full potential with LRC
    ! LRC plus ideal gas contribution plus total virial divided by V
    kin_q = kin_q / 1.5 ! Convert KE estimator to kinetic energy part of pressure
    p_f = variable_type ( nam = 'P full', val = pressure_lrc(rho,r_cut) + (kin_q+total%vir)/vol )

    ! Quantum spring energy per atom, just for interest's sake
    e_q = variable_type ( nam = 'Espring/N', val = total_spr/REAL(n) )

    ! Quantum polymer radius of gyration, just for interest's sake
    r_g = variable_type ( nam = 'Radius of gyration', val = rad_g )

    ! Collect together for averaging
    variables = [ r_r, c_r, e_f, p_f, e_q, k_q, r_g ]

  END FUNCTION calc_variables

  FUNCTION rad_gyr ( r ) RESULT ( r_g )
    IMPLICIT NONE
    REAL, DIMENSION(:,:,:), INTENT(in) :: r   ! Coordinate array
    REAL                               :: r_g ! RMS radius of gyration

    ! Calculate average radius of gyration of polymers
    ! The formula we use involves a single sweep over atoms, and is origin-independent
    ! To implement periodic boundaries, we take the origin on atom 1

    INTEGER            :: i, k
    REAL, DIMENSION(3) :: r_ik, r_cm
    REAL               :: r_sq

    r_g = 0.0 ! Zero function accumulator

    DO i = 1, n ! Loop over polymers
       r_cm = 0.0
       r_sq = 0.0
       DO k = 2, p ! Loop over atoms in polymer, omitting first one
          r_ik = r(:,i,k) - r(:,i,1) ! Position relative to atom 1
          r_ik = r_ik - ANINT(r_ik)  ! Position with PBC applied (box = 1 units)
          r_cm = r_cm + r_ik         ! Increment centre-of-mass accumulator
          r_sq = r_sq + SUM(r_ik**2) ! Increment squared distance accumulator
       END DO ! End loop over atoms in polymer, omitting first one
       r_cm = r_cm / REAL(p)      ! Normalize centre-of-mass vector
       r_sq = r_sq / REAL(p)      ! Normalize mean-squared distances
       r_sq = r_sq - SUM(r_cm**2) ! Definition of squared radius of gyration
       IF ( r_sq < 0.0 ) THEN
          WRITE(*,*) 'Warning! ', r_sq
          r_sq = 0.0
       END IF
       r_g = r_g + SQRT(r_sq)    ! Accumulate root-mean-square radius of gyration
    END DO ! End loop over polymers

    r_g = box * r_g / REAL(n) ! Average RMS Rg in sigma=1 units

  END FUNCTION rad_gyr

END PROGRAM qmc_pi_lj

