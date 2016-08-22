! mc_nvt_poly_lj.f90
! Monte Carlo, NVT ensemble, polyatomic molecule, LJ atoms
PROGRAM mc_nvt_poly_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module,  ONLY : read_cnf_mols, write_cnf_mols
  USE averages_module,   ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,    ONLY : metropolis, random_rotate_quaternion
  USE mc_poly_lj_module, ONLY : allocate_arrays, deallocate_arrays, energy_1, energy, q_to_d, &
       &                        n, na, r, e, d, ne

  IMPLICIT NONE

  ! Takes in a configuration of polyatomic molecules (positions and quaternions)
  ! For this example we specialize to triatomic molecules, natom=3
  ! Shifted LJ potential with no long range corrections
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1

  ! Most important variables
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: dr_max      ! maximum MC translational displacement
  REAL :: de_max      ! maximum MC rotational displacement
  REAL :: temperature ! specified temperature
  REAL :: r_cut       ! potential cutoff distance
  REAL :: rm_cut      ! molecule-molecule cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: vir         ! total virial
  REAL :: move_ratio  ! acceptance ratio of moves (to be averaged)
  REAL :: pressure    ! pressure (LJ sigma=1 units, to be averaged)
  REAL :: potential   ! potential energy per molecule (LJ sigma=1 units, to be averaged)

  LOGICAL            :: overlap
  INTEGER            :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL               :: pot_old, pot_new, vir_old, vir_new, delta

  REAL, DIMENSION(3)       :: ri   ! position of atom i
  REAL, DIMENSION(0:3)     :: ei   ! orientation of atom i
  INTEGER, PARAMETER       :: natom = 3
  REAL, DIMENSION(3,natom) :: di, db ! bond vectors
  REAL, DIMENSION(3)       :: zeta ! random numbers

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, r_cut, dr_max, de_max

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_nvt_poly_lj'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT, polyatomic, Lennard-Jones'
  WRITE( unit=output_unit, fmt='(a)' ) 'Results in units epsilon = sigma = 1'
  CALL time_stamp ( output_unit )

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  r_cut       = 2.5
  dr_max      = 0.15
  de_max      = 0.1
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_nvt_poly_lj'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature',               temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum r displacement',    dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Maximum e displacement',    de_max

  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box ) ! first call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of molecules',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box (in sigma units)', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Reduced density', density

  ! Body-fixed bond vectors in sigma units
  na = natom
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Number of atoms per molecule', na
  db(:,1) = [ 0.0, 0.0,  1.0/SQRT(3.0) ]
  db(:,2) = [-0.5, 0.0, -0.5/SQRT(3.0) ]
  db(:,3) = [ 0.5, 0.0, -0.5/SQRT(3.0) ]
  rm_cut = r_cut + 2.0 * SQRT ( MAXVAL ( SUM(db**2,dim=1) ) )
  DO i = 1, 3
     WRITE ( unit=output_unit, fmt='(a,i1,t40,3f15.5)' ) 'Body-fixed bond vector ', i, db(:,i)
  END DO
  WRITE( unit=output_unit, fmt='(a,t40,f15.5)') 'Molecule-molecule cutoff distance', rm_cut

  CALL allocate_arrays ( box, rm_cut )

  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! second call is to get r and e

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
  db     = db / box                  ! bond vectors

  ! Convert quaternions to space-fixed bond vectors
  DO i = 1, n
     d(:,:,i) = q_to_d ( e(:,i), db )
  END DO

  CALL energy ( box, r_cut, rm_cut, overlap, pot, vir )
  IF ( overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_poly_lj'
  END IF
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial potential energy (sigma units)', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)') 'Initial pressure (sigma units)',         pressure

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Potential', 'Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)

           ri(:)   = r(:,i)   ! copy old position
           ei(:)   = e(:,i)   ! copy old quaternion
           di(:,:) = d(:,:,i) ! copy old bond vectors
           CALL  energy_1 ( ri, di, i, ne, box, r_cut, rm_cut, overlap, pot_old, vir_old )
           IF ( overlap ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_nvt_poly_lj'
           END IF
           ri(:) = ri(:) + zeta * dr_max / box  ! trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )      ! periodic boundary correction
           ei = random_rotate_quaternion ( de_max, ei ) ! trial rotation
           di = q_to_d ( ei, db ) ! new space-fixed bond vectors (in box=1 units)
           CALL  energy_1 ( ri, di, i, ne, box, r_cut, rm_cut, overlap, pot_new, vir_new )

           IF ( .NOT. overlap ) THEN ! consider non-overlapping configuration
              delta = ( pot_new - pot_old ) / temperature
              IF ( metropolis ( delta ) ) THEN      ! accept Metropolis test
                 pot      = pot + pot_new - pot_old ! update potential energy
                 vir      = vir + vir_new - vir_old ! update virial
                 r(:,i)   = ri                      ! update position
                 e(:,i)   = ei                      ! update quaternion
                 d(:,:,i) = di                      ! update bond vectors
                 moves  = moves + 1                 ! increment move counter
              END IF ! reject Metropolis test
           END IF ! reject overlapping configuration

        END DO ! End loop over atoms

        ! Calculate all variables for this step
        move_ratio = REAL(moves) / REAL(n)
        potential  = pot / REAL(n)
        pressure   = density * temperature + vir / box**3
        CALL blk_add ( [move_ratio,potential,pressure] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk              ! number configuration by block
     CALL write_cnf_mols ( cnf_prefix//sav_tag, n, box, r*box, e ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final potential energy (sigma units)', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final pressure (sigma units)',         pressure

  CALL energy ( box, r_cut, rm_cut, overlap, pot, vir )
  IF ( overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_poly_lj'
  END IF
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE ( unit=output_unit, fmt='(a)'           ) 'Final check'
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final potential energy (sigma units)', potential
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Final pressure (sigma units)',         pressure

  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays

END PROGRAM mc_nvt_poly_lj

