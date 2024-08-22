! mc_nvt_sc.f90
! Monte Carlo, NVT ensemble, linear hard molecules
PROGRAM mc_nvt_sc

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

  ! Takes in a configuration of linear molecules (positions and orientations)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo for hard particles (the temperature is irrelevant)
  ! Uses no special neighbour lists
  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in reduced units kT=1

  ! Despite the program name, there is nothing here specific to spherocylinders
  ! The model is defined in mc_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor, &
       &                                    COMPILER_VERSION, COMPILER_OPTIONS
  
  USE config_io_module, ONLY : read_cnf_mols, write_cnf_mols
  USE averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : random_rotate_vector, random_translate_vector
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       overlap_1, overlap, n, r, e

  IMPLICIT NONE

  ! Most important variables
  REAL :: box     ! box length (in units where sigma=1)
  REAL :: dr_max  ! maximum MC displacement
  REAL :: de_max  ! maximum MC rotation
  REAL :: eps_box ! pressure scaling parameter

  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL, DIMENSION(3) :: ri, ei
  REAL               :: m_ratio

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, dr_max, de_max, eps_box

  WRITE ( unit=output_unit, fmt='(a)'   ) 'mc_nvt_sc'
  WRITE ( unit=output_unit, fmt='(2a)'  ) 'Compiler: ', COMPILER_VERSION()
  WRITE ( unit=output_unit, fmt='(2a/)' ) 'Options:  ', COMPILER_OPTIONS()
  
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT, hard linear molecules'
  CALL introduction

  CALL RANDOM_INIT ( .FALSE., .TRUE. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock  = 10
  nstep   = 10000
  dr_max  = 0.05
  de_max  = 0.05
  eps_box = 0.001

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_sc'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',           nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',  nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum displacement',       dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum rotation',           de_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Pressure scaling parameter', eps_box

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Box (in sigma units)', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',              REAL(n) / box**3
  CALL allocate_arrays
  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! Second call to get r and e
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial pressure and order calculation and overlap check
  IF ( overlap ( box ) ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_sc'
  END IF

  ! Initialize arrays for averaging and write column headings
  m_ratio = 0.0
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           ri(:) = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )                        ! Periodic boundary correction
           ei(:) = random_rotate_vector ( de_max, e(:,i) )        ! Trial move to new orientation

           IF ( .NOT. overlap_1 ( ri, ei, i, box ) ) THEN ! Accept
              r(:,i) = ri(:)     ! Update position
              e(:,i) = ei(:)     ! Update orientation
              moves  = moves + 1 ! Increment move counter
           END IF ! End accept

        END DO ! End loop over atoms

        m_ratio = REAL(moves) / REAL(n)

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                          ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk              ! Number configuration by block
     CALL write_cnf_mols ( cnf_prefix//sav_tag, n, box, r*box, e ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  ! Final overlap check
  IF ( overlap ( box ) ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_sc'
  END IF

  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE averages_module, ONLY : variable_type
    USE maths_module,    ONLY : nematic_order
    USE mc_module,       ONLY : n_overlap
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(3) :: variables ! The 3 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: m_r, p, order
    REAL                :: vol, rho, vir, ord

    ! Preliminary calculations (m_ratio, eps_box, box etc are already known)
    vol = box**3                                                   ! Volume
    rho = REAL(n) / vol                                            ! Density
    vir = REAL ( n_overlap ( box/(1.0+eps_box) ) ) / (3.0*eps_box) ! Virial
    ord = nematic_order ( e )                                      ! Order

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move acceptance ratio
    m_r = variable_type ( nam = 'Move ratio', val = m_ratio, instant = .FALSE. )

    ! Pressure in units kT/sigma**3
    ! Ideal gas contribution plus total virial divided by V
    p = variable_type ( nam = 'P', val = rho + vir/vol )

    ! Orientational order parameter
    order = variable_type ( nam = 'Nematic order', val = ord )

    ! Collect together for averaging
    variables = [ m_r, p, order ]

  END FUNCTION calc_variables

END PROGRAM mc_nvt_sc
