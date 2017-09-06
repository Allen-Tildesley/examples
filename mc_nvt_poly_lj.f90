! mc_nvt_poly_lj.f90
! Monte Carlo, NVT ensemble, polyatomic molecule
PROGRAM mc_nvt_poly_lj

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

  ! Takes in a configuration of polyatomic molecules (positions and quaternions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model (including the cutoff distance) is defined in mc_module

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_mols, write_cnf_mols
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               maths_module,     ONLY : metropolis, random_rotate_quaternion, random_translate_vector, q_to_a
  USE               mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     potential_1, potential, n, na, db, r, e, d, potential_type

  IMPLICIT NONE

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: dr_max      ! Maximum MC translational displacement
  REAL :: de_max      ! Maximum MC rotational displacement
  REAL :: temperature ! Specified temperature

  ! Composite interaction = pot & vir & ovr variables
  TYPE(potential_type) :: total, partial_old, partial_new

  INTEGER               :: blk, stp, i, a, nstep, nblock, moves, ioerr
  REAL                  :: delta, m_ratio
  REAL, DIMENSION(3)    :: ri
  REAL, DIMENSION(0:3)  :: ei
  REAL, DIMENSION(3,3)  :: ai
  REAL, DIMENSION(na,3) :: di

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, dr_max, de_max

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_nvt_poly_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble, polyatomic molecule'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 20000
  temperature = 1.0
  dr_max      = 0.05
  de_max      = 0.05

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_poly_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum r displacement',    dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum e displacement',    de_max

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of molecules',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density of molecules',  REAL(n) / box**3
  CALL allocate_arrays ( box )
  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! Second call is to get r and e
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Calculate all bond vectors
  DO i = 1, n ! Loop over molecules
     ai = q_to_a ( e(:,i) ) ! Rotation matrix for i
     DO a = 1, na ! Loop over all atoms
        d(:,a,i) = MATMUL ( db(:,a), ai ) ! NB: equivalent to ai_T*db, ai_T=transpose of ai
     END DO ! End loop over all atoms
  END DO ! End loop over molecules

  ! Initial energy and overlap check
  total = potential ( box )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_poly_lj'
  END IF

  ! Initialize arrays for averaging and write column headings
  m_ratio = 0.0
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           partial_old = potential_1 ( r(:,i), d(:,:,i), i, box ) ! Old molecule potential, virial, etc

           IF ( partial_old%ovr ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_nvt_poly_lj'
           END IF

           ri = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
           ri = ri - ANINT ( ri )                              ! Periodic boundary correction
           ei = random_rotate_quaternion ( de_max, e(:,i) )    ! Trial rotation
           ai = q_to_a ( ei ) ! Rotation matrix for i
           DO a = 1, na ! Loop over all atoms
              di(:,a) = MATMUL ( db(:,a), ai ) ! NB: equivalent to ai_T*db, ai_T=transpose of ai
           END DO ! End loop over all atoms

           partial_new = potential_1 ( ri, di, i, box ) ! New molecule potential, virial etc

           IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

              delta = partial_new%pot - partial_old%pot ! Use cut-and-shifted potential
              delta = delta / temperature               ! Divide by temperature

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 total    = total + partial_new - partial_old ! Update potential energy
                 r(:,i)   = ri                                ! Update position
                 e(:,i)   = ei                                ! Update quaternion
                 d(:,:,i) = di                                ! Update bond vectors
                 moves    = moves + 1                         ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End test for non-overlapping configuration

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

  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE averages_module, ONLY : variable_type
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(3) :: variables ! The 3 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    ! In this example we simulate using the shifted-force potential only
    ! The values of < p_sf >, < e_sf > and density should be consistent (for this potential)
    ! There are no long-range or delta corrections

    TYPE(variable_type) :: m_r, e_sf, p_sf
    REAL                :: vol, rho

    ! Preliminary calculations
    vol = box**3        ! Volume
    rho = REAL(n) / vol ! Density

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move acceptance ratio
    m_r = variable_type ( nam = 'Move ratio', val = m_ratio, instant = .FALSE. )

    ! Internal energy per molecule (shifted-force potential)
    ! Ideal gas contribution (assuming nonlinear molecules) plus total PE divided by N
    e_sf = variable_type ( nam = 'E/N shifted force', val = 3.0*temperature + total%pot/REAL(n) )

    ! Pressure (shifted-force potential)
    ! Ideal gas contribution plus total virial divided by V
    p_sf = variable_type ( nam = 'P shifted force', val = rho*temperature + total%vir/vol )

    ! Collect together for averaging
    variables = [ m_r, e_sf, p_sf ]

  END FUNCTION calc_variables

END PROGRAM mc_nvt_poly_lj

