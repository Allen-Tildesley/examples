! mc_chain_nvt_cbmc_lj.f90
! Monte Carlo, single chain, NVT, CBMC
PROGRAM mc_chain_nvt_cbmc_lj
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, regrow, n, r

  IMPLICIT NONE

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

  ! Most important variables
  REAL    :: temperature ! Temperature (specified, in units of well depth)
  REAL    :: bond        ! Intramolecular bond length
  REAL    :: k_spring    ! Strength of intramolecular bond springs
  INTEGER :: m_max       ! Maximum atoms in regrow
  INTEGER :: k_max       ! Number of random tries per atom in regrow

  ! Quantities for averaging
  REAL :: m_ratio ! Acceptance ratio for regrowth moves
  REAL :: pe      ! Total nonbonded + bond spring potential energy
  REAL :: r_g     ! Radius of gyration

  INTEGER :: blk, stp, nstep, nblock, ioerr
  LOGICAL :: accepted

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, m_max, k_max, temperature, k_spring

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_chain_nvt_cbmc_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble, CBMC, chain molecule'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Simulation uses full nonbonded potential (no cutoff)'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing

  nblock         = 10    ! Number of blocks
  nstep          = 10000 ! Number of steps per block
  m_max          = 3     ! Maximum atoms in regrow
  k_max          = 32    ! Number of random tries per atom in regrow
  temperature    = 1.0   ! Temperature (in units of well depth)
  k_spring       = 400.0 ! Strength of intramolecular bond springs (same units)

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_chain_nvt_cbmc_lj'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',                nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',       nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Max atoms in regrow',             m_max
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Random tries per atom in regrow', k_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature/well depth',          temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Bond spring strength',            k_spring

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond ) ! First call is just to get n and bond
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',          n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Bond length (in sigma units)', bond

  CALL allocate_arrays

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, bond, r ) ! Second call gets r

  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Regrow ratio', 'PE', 'Rg' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL regrow ( temperature, m_max, k_max, bond, k_spring, accepted )
        IF ( accepted ) THEN
           m_ratio = 1.0
        ELSE
           m_ratio = 0.0
        END IF

        CALL calculate()
        CALL blk_add ( [m_ratio,pe,r_g] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk         ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    USE mc_module, ONLY : potential, spring_pot, potential_type
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    REAL, DIMENSION(3)   :: r_cm
    TYPE(potential_type) :: total

    total = potential ( ) ! Calculate nonbonded potential with overlap flag

    IF ( total%ovr ) THEN ! Overlap test (might happen with initial configuration)
       WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
       STOP 'Error in mc_chain_nvt_cbmc_lj/calculate'
    END IF ! End overlap test

    pe   = total%pot + spring_pot ( bond, k_spring )
    r_cm = SUM ( r, dim=2 ) / REAL(n) ! Centre of mass
    r_g  = SQRT ( SUM ( ( r - SPREAD(r_cm,dim=2,ncopies=n) ) ** 2 ) / REAL(n) )

    IF ( PRESENT ( string ) ) THEN ! output required
       WRITE ( unit=output_unit, fmt='(a)'          ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'PE', pe
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Rg', r_g
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_chain_nvt_cbmc_lj

