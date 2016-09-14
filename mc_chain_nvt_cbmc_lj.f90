! mc_chain_nvt_cbmc_lj.f90
! Monte Carlo, single chain, NVT, CBMC
PROGRAM mc_chain_nvt_cbmc_lj
  !
  ! TODO MPA provide code
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_module,        ONLY : model_description, allocate_arrays, deallocate_arrays, &
       &                       regrow, energy, &
       &                       n, bond, r

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
  REAL    :: temperature ! temperature (in units of well depth)
  REAL    :: k_spring    ! strength of intramolecular bond springs
  REAL    :: spring_pot  ! bond spring potential energy (for averaging)
  REAL    :: pot         ! nonbonded potential energy
  REAL    :: potential   ! nonbonded potential energy per atom (for averaging)
  REAL    :: move_ratio  ! acceptance ratio for regrowth moves
  INTEGER :: m_max       ! maximum atoms in regrow
  INTEGER :: k_max       ! number of random tries per atom in regrow

  INTEGER :: blk, stp, nstep, nblock, ioerr
  LOGICAL :: overlap

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, m_max, k_max, temperature, k_spring

  WRITE ( unit=output_unit, fmt='(a)' ) 'mc_chain_nvt_cbmc_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble, CBMC, chain molecule'
  CALL model_description ( output_unit )
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing

  nblock         = 10    ! number of blocks
  nstep          = 10000 ! number of sweeps per block
  m_max          = 3     ! maximum atoms in regrow
  k_max          = 32    ! number of random tries per atom in regrow
  temperature    = 1.0   ! temperature (in units of well depth)
  k_spring       = 10.0  ! strength of intramolecular bond springs

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
  potential  = pot / REAL ( n )
  spring_pot = 0.5 * k_spring * SUM((r(:,1:n-1)-r(:,2:n))**2) / REAL(n)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial potential energy', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial spring energy',    spring_pot

  CALL run_begin ( [ CHARACTER(len=15) :: 'Regrow ratio', 'Energy', 'Spring' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL regrow ( temperature, m_max, k_max, k_spring, pot, move_ratio )
        potential  = pot / REAL ( n )
        spring_pot = 0.5 * k_spring * SUM((r(:,1:n-1)-r(:,2:n))**2) / REAL(n)
        CALL blk_add ( [move_ratio,potential,spring_pot] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk         ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, bond, r ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  potential  = pot / REAL(n)
  spring_pot = 0.5 * k_spring * SUM((r(:,1:n-1)-r(:,2:n))**2) / REAL(n)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Final potential energy', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Final spring energy',    spring_pot

  CALL energy ( overlap, pot )
  IF ( overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_chain_nvt_cbmc_lj'
  END IF
  potential = pot / REAL ( n )
  WRITE ( unit=output_unit, fmt='(a)'           ) 'Final check'
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final potential energy', potential

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, bond, r )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays

END PROGRAM mc_chain_nvt_cbmc_lj

