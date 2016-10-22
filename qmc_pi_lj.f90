! qmc_pi_lj.f90
! Quantum Monte Carlo, path-integral method
PROGRAM qmc_pi_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis, random_translate_vector
  USE qmc_module,       ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       potential_1, spring_1, potential, spring, n, p, r, &
       &                       potential_type, OPERATOR(+), OPERATOR(-)

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
  REAL :: box         ! box length
  REAL :: lambda      ! de Boer length
  REAL :: density     ! density
  REAL :: dr_max      ! maximum MC displacement
  REAL :: temperature ! specified temperature
  REAL :: r_cut       ! potential cutoff distance

  ! Quantities to be averaged
  REAL :: m_ratio ! acceptance ratio of moves
  REAL :: en_c    ! Internal energy per atom for simulated, cut, potential
  REAL :: en      ! Internal energy per atom for full potential with LRC

  ! Composite interaction = pot & overlap variables
  TYPE(potential_type) :: system, atom_old, atom_new

  INTEGER            :: blk, stp, i, k, nstep, nblock, moves, ioerr
  REAL               :: system_spring, atom_old_spring, atom_new_spring, delta, k_spring
  REAL, DIMENSION(3) :: rik

  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav'       ! may be overwritten with block number
  CHARACTER(len=6)            :: cnf_prefix = 'cnf##.' ! will have polymer id inserted
  CHARACTER(len=2)            :: k_tag                 ! will contain polymer id

  NAMELIST /nml/ nblock, nstep, p, temperature, r_cut, dr_max, lambda

  WRITE ( unit=output_unit, fmt='(a)' ) 'qmc_pi_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Path-integral Monte Carlo, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Simulation uses cut (but not shifted) potential'
  CALL introduction ( output_unit )
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

  k_spring = REAL(p) * ( temperature / lambda ) ** 2
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Quantum spring constant', k_spring

  ! Read first polymer configuration to get basic parameters n and box
  WRITE(k_tag,fmt='(i2.2)') 1 ! convert into character form
  cnf_prefix(4:5) = k_tag     ! insert into configuration filename
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box )
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density

  CALL allocate_arrays ( box, r_cut ) ! Allocate r

  ! Read each ring polymer from a unique file; we assume that n and box are compatible!
  DO k = 1, p
     WRITE(k_tag,fmt='(i2.2)') k ! convert into character form
     cnf_prefix(4:5) = k_tag     ! insert into configuration filename
     CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r(:,:,k) )
  END DO

  r(:,:,:) = r(:,:,:) / box                ! Convert positions to box units
  r(:,:,:) = r(:,:,:) - ANINT ( r(:,:,:) ) ! Periodic boundaries (box=1 units)

  ! Calculate classical LJ and quantum spring potential energies
  system = potential ( box, r_cut )
  IF ( system%overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in qmc_pi_lj'
  END IF
  system_spring = spring ( box, k_spring )
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'E/N (cut)', 'E/N (full)' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           DO k = 1, p ! Begin loop over ring polymers

              atom_old = potential_1 ( r(:,i,k), i, k, box, r_cut ) ! Old atom classical potential etc

              IF ( atom_old%overlap ) THEN ! should never happen
                 WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                 STOP 'Error in qmc_pi_lj'
              END IF

              atom_old_spring = spring_1 ( r(:,i,k), i, k, box, k_spring ) ! Old atom quantum potential

              rik(:) = random_translate_vector ( dr_max/box, r(:,i,k) ) ! Trial move to new position (in box=1 units) 
              rik(:) = rik(:) - ANINT ( rik(:) )                        ! Periodic boundary correction

              atom_new = potential_1 ( rik, i, k, box, r_cut ) ! New atom classical potential etc

              IF ( .NOT. atom_new%overlap ) THEN ! Test for non-overlapping configuration

                 atom_new_spring = spring_1 ( rik, i, k, box, k_spring ) ! New atom quantum potential

                 delta = ( atom_new%pot_c + atom_new_spring - atom_old%pot_c - atom_old_spring ) / temperature

                 IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                    system        = system        + atom_new        - atom_old        ! Update classical system values
                    system_spring = system_spring + atom_new_spring - atom_old_spring ! Update quantum system potential
                    r(:,i,k) = rik                                                    ! Update position
                    moves    = moves + 1                                              ! Increment move counter
                 END IF ! End accept Metropolis test

              END IF ! End test for non-overlapping configuration

           END DO ! End loop over ring polymers

        END DO ! End loop over atoms

        m_ratio = REAL(moves) / REAL(n)

        ! Calculate all variables for this step
        CALL calculate ( )
        CALL blk_add ( [m_ratio,en_c,en] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,fmt='(i3.3)') blk ! number configuration by block

     ! Save each ring polymer to a unique file
     DO k = 1, p
        WRITE(k_tag,fmt='(i2.2)') k ! convert into character form
        cnf_prefix(4:5) = k_tag     ! insert into configuration filename
        CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r(:,:,k)*box )
     END DO

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  system = potential ( box, r_cut )
  IF ( system%overlap ) THEN ! this should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in qmc_pi_lj'
  END IF
  system_spring = spring ( box, k_spring )
  CALL calculate ( 'Final check' )

  ! Write each ring polymer to a unique file
  DO k = 1, p
     WRITE(k_tag,fmt='(i2.2)') k ! convert into character form
     cnf_prefix(4:5) = k_tag     ! insert into configuration filename
     CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r(:,:,k)*box )
  END DO
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS
  
  SUBROUTINE calculate ( string )
    USE qmc_module, ONLY : potential_lrc
    IMPLICIT NONE
    CHARACTER (len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates the properties of interest from system values
    ! and optionally writes them out (e.g. at the start and end of the run)
    ! In this example we simulate using the cut (but not shifted) potential
    ! The values of < en_c > and < density > should be consistent (for this potential)
    ! For comparison, long-range corrections are also applied to give
    ! an estimate of < en > for the full (uncut) potential
    ! The value of the cut-and-shifted potential pot_s is not used, in this example

    real :: kin

    kin   = 1.5 * n * p * temperature                        ! Kinetic energy
    en_c  = ( kin + system%pot_c - system_spring ) / REAL(n) ! Total energy per atom (cut but not shifted)
    en    = en_c + potential_lrc ( density, r_cut )          ! Add long-range contribution to PE/N

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (cut)',  en_c
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (full)', en
    END IF

  END SUBROUTINE calculate

END PROGRAM qmc_pi_lj

