! qmc_pi_lj.f90
! Quantum Monte Carlo, path-integral, LJ atoms
PROGRAM qmc_pi_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE utility_module, ONLY : metropolis, read_cnf_atoms, write_cnf_atoms, time_stamp, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add
  USE qmc_pi_lj_module, ONLY : initialize, finalize, energy_cl_1, energy_qu_1, energy_cl, energy_qu, move, &
       &                     n, p, r, ne
  IMPLICIT NONE

  ! Takes in a set of configurations of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts path-integral Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1
  ! The importance of quantum effects is specified through the reduced de Boer length
  ! lambda = (hbar/sqrt(mass*epsilon))/sigma which takes typical values of
  ! 0.01 for Xe, 0.03 for Ar, and 0.095 for Ne.
  ! This means that the quantum spring constant may be expressed k_spring = P*T**2/lambda**2
  ! where T is the reduced temperature kB*T/epsilon

  ! Most important variables
  REAL :: sigma        ! atomic diameter (in units where box=1)
  REAL :: box          ! box length (in units where sigma=1)
  REAL :: lambda       ! de Boer length (in units where sigma=1)
  REAL :: density      ! reduced density n*sigma**3/box**3
  REAL :: dr_max       ! maximum MC displacement
  REAL :: temperature  ! specified temperature
  REAL :: r_cut        ! potential cutoff distance
  REAL :: pot_cl       ! classical potential energy
  REAL :: pot_qu       ! quantum potential energy
  REAL :: move_ratio   ! acceptance ratio of moves (to be averaged)
  REAL :: potential_cl ! classical potential energy per atom (LJ sigma=1 units, to be averaged)
  REAL :: potential_qu ! quantum potential energy per atom (LJ sigma=1 units, to be averaged)
  REAL :: energy       ! total energy per atom (LJ sigma=1 units, to be averaged)

  LOGICAL            :: overlap
  INTEGER            :: blk, stp, i, k, nstep, nblock, moves, ioerr
  REAL               :: pot_cl_old, pot_cl_new, pot_qu_old, pot_qu_new, kin, delta, k_spring
  REAL, DIMENSION(3) :: rik  ! position of atom i in polymer k
  REAL, DIMENSION(3) :: zeta ! random numbers

  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number
  CHARACTER(len=6)            :: cnf_prefix = 'cnf##.' ! will have polymer id inserted
  CHARACTER(len=2)            :: k_tag           ! will contain polymer id

  NAMELIST /nml/ nblock, nstep, p, temperature, r_cut, dr_max, lambda

  WRITE( unit=output_unit, fmt='(a)' ) 'qmc_pi_lj'
  WRITE( unit=output_unit, fmt='(a)' ) 'Path-integral Monte Carlo, constant-NVT, Lennard-Jones'
  WRITE( unit=output_unit, fmt='(a)' ) 'Results in units epsilon = sigma = 1'
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  r_cut       = 2.5
  dr_max      = 0.15
  p           = 4
  lambda      = 0.1
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in qmc_pi_lj'
  END IF

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of ring polymers',   p
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum displacement',      dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'de Boer length',            lambda

  IF ( p < 2 .OR. p > 99 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'p must lie between 2 and 99', p
     STOP 'Error in qmc_pi_lj'
  END IF

  k_spring = REAL(p) * ( temperature/lambda )**2
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'quantum spring constant', k_spring

  WRITE(k_tag,fmt='(i2.2)') 1 ! convert into character form
  cnf_prefix(4:5) = k_tag     ! insert into configuration filename

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! first call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box (in sigma units)', box
  sigma   = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Reduced density', density

  ! Convert run and potential parameters to box units
  sigma  = sigma / box
  r_cut  = r_cut / box
  dr_max = dr_max / box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'sigma (in box units)', sigma
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'r_cut (in box units)', r_cut
  IF ( r_cut > 0.5 ) THEN
     WRITE ( unit=error_unit, fmt='(a,f15.5)') 'r_cut too large ', r_cut
     STOP 'Error in qmc_pi_lj'
  END IF

  CALL initialize ! Allocate r

  ! read each ring polymer from a unique file
  DO k = 1, p
     WRITE(k_tag,fmt='(i2.2)') k ! convert into character form
     cnf_prefix(4:5) = k_tag     ! insert into configuration filename
     CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r(:,:,k) )
  END DO

  ! Convert to box units
  r(:,:,:) = r(:,:,:) / box
  r(:,:,:) = r(:,:,:) - ANINT ( r(:,:,:) ) ! Periodic boundaries

  ! Calculate classical LJ potential energy
  CALL energy_cl ( sigma, r_cut, overlap, pot_cl )
  IF ( overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in qmc_pi_lj'
  END IF

  ! Calculate quantum spring potential energy
  CALL energy_qu ( k_spring, sigma, pot_qu )
  kin   = 1.5 * n * p * temperature  ! Kinetic energy
  potential_cl = pot_cl / REAL(n)
  potential_qu = pot_qu / REAL(n)
  energy = ( kin + pot_cl - pot_qu ) / REAL(n) ! total energy per atom
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial classical potential energy', potential_cl
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial quantum potential energy',   potential_qu
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Initial energy',                     energy

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Pot-classical', 'Pot-quantum', 'Energy' ] )


  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           DO k = 1, p ! Begin loop over ring polymers

              CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
              zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

              rik(:) = r(:,i,k)
              CALL energy_cl_1 ( rik, i, k, ne, sigma, r_cut, overlap, pot_cl_old )
              IF ( overlap ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in qmc_pi_lj'
              END IF
              CALL energy_qu_1 ( rik, i, k, ne, k_spring, sigma, pot_qu_old )

              rik(:) = rik(:) + zeta * dr_max   ! trial move to new position
              rik(:) = rik(:) - ANINT ( rik(:) ) ! periodic boundary correction

              CALL energy_cl_1 ( rik, i, k, ne, sigma, r_cut, overlap, pot_cl_new )

              IF ( .NOT. overlap ) THEN ! consider non-overlapping configuration
                 CALL energy_qu_1 ( rik, i, k, ne, k_spring, sigma, pot_qu_new )
                 delta = ( pot_cl_new + pot_qu_new - pot_cl_old - pot_qu_old ) / temperature
                 IF ( metropolis ( delta ) ) THEN ! accept Metropolis test
                    pot_cl = pot_cl + pot_cl_new - pot_cl_old ! update classical LJ potential energy
                    pot_qu = pot_qu + pot_qu_new - pot_qu_old ! update quantum spring potential energy
                    CALL move ( i, k, rik )                   ! update position
                    moves  = moves + 1                        ! increment move counter
                 END IF ! reject Metropolis test
              END IF ! reject overlapping configuration

           END DO ! End loop over ring polymers

        END DO ! End loop over atoms

        ! Calculate all variables for this step
        move_ratio = REAL(moves) / REAL(n)
        potential_cl = pot_cl / REAL(n)
        potential_qu = pot_qu / REAL(n)
        energy = ( kin + pot_cl - pot_qu ) / REAL(n) ! total energy per atom
        CALL blk_add ( [move_ratio,potential_cl,potential_qu,energy] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk        ! number configuration by block

     ! save each ring polymer to a unique file
     DO k = 1, p
        WRITE(k_tag,fmt='(i2.2)') k ! convert into character form
        cnf_prefix(4:5) = k_tag     ! insert into configuration filename
        CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r(:,:,k)*box ) ! save configuration
     END DO

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  potential_cl = pot_cl / REAL(n)
  potential_qu = pot_qu / REAL(n)
  energy = ( kin + pot_cl - pot_qu ) / REAL(n) ! total energy per atom
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final classical potential energy', potential_cl
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final quantum potential energy',   potential_qu
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final energy',                     energy
  CALL energy_cl ( sigma, r_cut, overlap, pot_cl )
  IF ( overlap ) THEN ! this should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in qmc_pi_lj'
  END IF
  CALL energy_qu ( k_spring, sigma, pot_qu )
  WRITE ( unit=output_unit, fmt='(a)'           ) 'Final check'
  potential_cl = pot_cl / REAL(n)
  potential_qu = pot_qu / REAL(n)
  energy = ( kin + pot_cl - pot_qu ) / REAL(n) ! total energy per atom
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final classical potential energy', potential_cl
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final quantum potential energy',   potential_qu
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final energy',                     energy

  ! write each ring polymer to a unique file
  DO k = 1, p
     WRITE(k_tag,fmt='(i2.2)') k ! convert into character form
     cnf_prefix(4:5) = k_tag     ! insert into configuration filename
     CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r(:,:,k)*box )
  END DO
  CALL time_stamp ( output_unit )

  CALL finalize

END PROGRAM qmc_pi_lj

