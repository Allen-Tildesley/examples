! mc_npt_hs.f90
! Monte Carlo simulation, NPT ensemble, hard spheres
PROGRAM mc_npt_hs

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

  ! Takes in a configuration of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at given NPT (the temperature is irrelevant)
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! We take kT=1 throughout defining the unit of energy
  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model,
  ! in this case, for hard spheres, sigma = 1

  ! The logarithm of the box length is sampled uniformly

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               maths_module,     ONLY : metropolis, random_translate_vector
  USE               mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     overlap_1, overlap, n, r

  IMPLICIT NONE

  ! Most important variables
  REAL :: box      ! Box length
  REAL :: dr_max   ! Maximum MC displacement
  REAL :: db_max   ! Maximum MC box displacement
  REAL :: pressure ! Specified pressure

  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL, DIMENSION(3) :: ri
  REAL               :: box_scale, box_new, den_scale, delta, zeta, m_ratio, v_ratio

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, dr_max, db_max, pressure

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_npt_hs'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NPT'
  CALL introduction

  CALL RANDOM_INIT ( .FALSE., .TRUE. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock   = 10
  nstep    = 10000
  dr_max   = 0.15
  db_max   = 0.005
  pressure = 4.0

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_npt_hs'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Pressure',                  pressure
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum displacement',      dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum box displacement',  db_max

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  CALL allocate_arrays ! Allocates r
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r ) ! Second call to get r
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial calculation and overlap check
  IF ( overlap ( box ) ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_npt_hs'
  END IF

  ! Initialize arrays for averaging and write column headings
  m_ratio = 0.0
  v_ratio = 0.0
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           ri(:) = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )                        ! Periodic boundary correction

           IF ( .NOT. overlap_1 ( ri, i, box ) ) THEN ! Accept
              r(:,i) = ri(:)     ! Update position
              moves  = moves + 1 ! Increment move counter
           END IF ! End accept

        END DO ! End loop over atoms

        m_ratio = REAL(moves) / REAL(n)

        v_ratio = 0.0                   ! Zero volume move counter
        CALL RANDOM_NUMBER ( zeta )     ! Uniform random number in range (0,1)
        zeta      = 2.0*zeta - 1.0      ! Now in range (-1,+1)
        box_scale = EXP ( zeta*db_max ) ! Sampling log(box) and log(vol) uniformly
        box_new   = box * box_scale     ! New box (in sigma units)
        den_scale = 1.0 / box_scale**3  ! Density scaling factor

        IF ( .NOT. overlap ( box_new ) ) THEN ! Test for non-overlapping configuration

           delta = pressure * ( box_new**3 - box**3 ) ! PV term (temperature = 1.0 )
           delta = delta + REAL(n+1) * LOG(den_scale) ! Factor (n+1) consistent with log(box) sampling

           IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
              box     = box_new  ! Update box
              v_ratio = 1.0      ! Set move counter
           END IF ! reject Metropolis test

        END IF ! End test for overlapping configuration

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                        ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk            ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  ! Final overlap check and pressure calculation
  IF ( overlap ( box ) ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_npt_hs'
  END IF

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE averages_module, ONLY : variable_type
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(3) :: variables ! The 3 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: m_r, v_r, density
    REAL                :: vol, rho

    ! Preliminary calculations (m_ratio, v_ratio, box etc already known)
    vol = box**3        ! Volume
    rho = REAL(n) / vol ! Density

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move acceptance ratios
    m_r = variable_type ( nam = 'Move ratio',   val = m_ratio, instant = .FALSE. )
    v_r = variable_type ( nam = 'Volume ratio', val = v_ratio, instant = .FALSE. )

    ! Density
    density = variable_type ( nam = 'Density', val = rho )

    ! Collect together for averaging
    variables = [ m_r, v_r, density ]

  END FUNCTION calc_variables

END PROGRAM mc_npt_hs
