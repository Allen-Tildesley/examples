! mc_chain_nvt_cbmc_lj.f90
! Monte Carlo, single chain, NVT, CBMC
PROGRAM mc_chain_nvt_cbmc_lj
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       regrow, energy, spring_pot, n, r

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
  REAL    :: pot_bonds   ! Bond spring potential energy per atom (for averaging)
  REAL    :: pot         ! Nonbonded potential energy
  REAL    :: potential   ! Nonbonded potential energy per atom (for averaging)
  REAL    :: move_ratio  ! Acceptance ratio for regrowth moves (for averaging)
  REAL    :: r_gyration  ! Radius of gyration (for averaging)
  INTEGER :: m_max       ! Maximum atoms in regrow
  INTEGER :: k_max       ! Number of random tries per atom in regrow

  INTEGER :: blk, stp, nstep, nblock, ioerr
  LOGICAL :: overlap

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, m_max, k_max, temperature, k_spring

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_chain_nvt_cbmc_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble, CBMC, chain molecule'
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

  CALL energy ( overlap, pot )
  IF ( overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_chain_nvt_cbmc_lj'
  END IF

  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Regrow ratio', 'Nonbonded Pot', 'Spring Pot', 'Rad Gyration' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL regrow ( temperature, m_max, k_max, bond, k_spring, pot, move_ratio )

        ! Calculate all variables for this step
        CALL calculate()
        CALL blk_add ( [move_ratio,potential,pot_bonds,r_gyration] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk         ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    REAL, DIMENSION(3) :: r_cm

    potential  = pot / REAL ( n )
    pot_bonds  = spring_pot ( bond, k_spring ) / REAL(n)
    r_cm       = SUM ( r, dim=2 ) / REAL(n) ! Centre of mass
    r_gyration = SQRT ( SUM ( ( r - SPREAD(r_cm,dim=2,ncopies=n) ) ** 2 ) / REAL(n) )

    IF ( PRESENT ( string ) ) THEN ! output required
       WRITE ( unit=output_unit, fmt='(a)'          ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Nonbonded Pot', potential
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Spring Pot',    pot_bonds
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Rad Gyration',  r_gyration
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_chain_nvt_cbmc_lj

