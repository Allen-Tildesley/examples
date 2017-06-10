! mc_chain_nvt_cbmc_lj.f90
! Monte Carlo, single chain, NVT, CBMC
PROGRAM mc_chain_nvt_cbmc_lj

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
  ! Conducts Monte Carlo, NVT ensemble using CBMC regrowth moves
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Input configuration, output configuration, all calculations, and all results 
  ! are given in simulation units defined by the model.
  ! E.g. for Lennard-Jones, atomic diameter sigma = 1, well-depth epsilon=1
  ! Configurational weights are calculated on the basis of the nonbonded interactions

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in mc_module

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, regrow, n, r

  IMPLICIT NONE

  ! Most important variables
  REAL    :: temperature ! Temperature (specified, in units of well depth)
  REAL    :: bond        ! Intramolecular bond length
  REAL    :: k_spring    ! Strength of intramolecular bond springs
  INTEGER :: m_max       ! Maximum atoms in regrow
  INTEGER :: k_max       ! Number of random tries per atom in regrow

  INTEGER :: blk, stp, nstep, nblock, ioerr
  LOGICAL :: accepted
  REAL    :: m_ratio

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, m_max, k_max, temperature, k_spring

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_chain_nvt_cbmc_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble, CBMC, chain molecule'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Simulation uses full nonbonded potential (no cutoff)'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10     ! Number of blocks
  nstep       = 100000 ! Number of steps per block
  m_max       = 3      ! Maximum atoms in regrow
  k_max       = 32     ! Number of random tries per atom in regrow
  temperature = 1.0    ! Temperature (in units of well depth)
  k_spring    = 400.0  ! Strength of intramolecular bond springs (same units)

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_chain_nvt_cbmc_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',                nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',       nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Max atoms in regrow',             m_max
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Random tries per atom in regrow', k_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature/well depth',          temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Bond spring strength',            k_spring

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond ) ! First call is just to get n and bond
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',          n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Bond length (in sigma units)', bond
  CALL allocate_arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r ) ! Second call gets r

  ! Initialize arrays for averaging and write column headings
  m_ratio = 0.0
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL regrow ( temperature, m_max, k_max, bond, k_spring, accepted )
        IF ( accepted ) THEN
           m_ratio = 1.0
        ELSE
           m_ratio = 0.0
        END IF

        ! Calculate and accumulate quantities for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                     ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk         ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r )

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE mc_module,       ONLY : potential, spring_pot, potential_type
    USE averages_module, ONLY : variable_type, msd
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(4) :: variables ! The 4 variables listed below

    ! This function returns all variables of interest in an array, for use in the main program

    TYPE(variable_type)  :: m_r, e_x, r_g, c_x
    TYPE(potential_type) :: total
    REAL, DIMENSION(3)   :: rcm
    REAL                 :: spr, rsq

    ! Preliminary calculations
    total = potential ( ) ! Nonbonded potential with overlap flag
    IF ( total%ovr ) THEN ! Overlap test (might happen with initial configuration)
       WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
       STOP 'Error in mc_chain_nvt_cbmc_lj/calculate'
    END IF ! End overlap test
    spr = spring_pot ( bond, k_spring )                              ! Total spring potential energy
    rcm = SUM ( r, dim=2 ) / REAL(n)                                 ! Centre of mass
    rsq = SUM ( ( r - SPREAD(rcm,dim=2,ncopies=n) ) ** 2 ) / REAL(n) ! Mean-squared distance from CM

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move acceptance ratio
    m_r = variable_type ( nam = 'Regrow ratio', val = m_ratio, instant = .FALSE. )

    ! Total potential energy per atom (excess, without ideal gas contribution)
    ! Total PE of bond springs plus total LJ PE (not cut, nor shifted) divided by N
    e_x = variable_type ( nam = 'PE/N', val = (spr+total%pot)/REAL(n) )

    ! Radius of gyration
    r_g = variable_type ( nam = 'Rg', val = SQRT(rsq) )

    ! Heat Capacity per atom (excess, without ideal gas contribution)
    ! MSD of PE / (sqrt(N)*T)
    ! Total PE of bond springs plus total LJ PE (not cut, nor shifted), divided by T
    c_x = variable_type ( nam = 'Cv(ex)/N', val = (spr+total%pot)/(SQRT(REAL(n))*temperature), &
         &                method = msd, instant = .FALSE. )

    ! Collect together for averaging
    variables = [ m_r, e_x, r_g, c_x ]

  END FUNCTION calc_variables

END PROGRAM mc_chain_nvt_cbmc_lj
