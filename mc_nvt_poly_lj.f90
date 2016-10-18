! mc_nvt_poly_lj.f90
! Monte Carlo, NVT ensemble, polyatomic molecule
PROGRAM mc_nvt_poly_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_mols, write_cnf_mols
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis, random_rotate_quaternion
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       energy_1, energy, q_to_d, n, na, r, e, d, potential_type

  IMPLICIT NONE

  ! Takes in a configuration of polyatomic molecules (positions and quaternions)
  ! For this example we specialize to triatomic molecules, natom=3
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r and bond vectors db are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in mc_module

  ! Most important variables
  REAL :: box         ! box length
  REAL :: density     ! density
  REAL :: dr_max      ! maximum MC translational displacement
  REAL :: de_max      ! maximum MC rotational displacement
  REAL :: temperature ! specified temperature
  REAL :: r_cut       ! potential cutoff distance
  REAL :: rm_cut      ! molecule-molecule cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: vir         ! total virial
  REAL :: move_ratio  ! acceptance ratio of moves (to be averaged)
  REAL :: pres_virial ! virial pressure (to be averaged)
  REAL :: potential   ! potential energy per molecule (to be averaged)

  TYPE(potential_type) :: eng_old, eng_new ! Composite energy = pot & vir & overlap variables
  INTEGER      :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL         :: delta

  INTEGER, PARAMETER       :: natom = 3
  REAL, DIMENSION(3)       :: ri     ! position of molecule i
  REAL, DIMENSION(0:3)     :: ei     ! orientation (quaternion) of molecule i
  REAL, DIMENSION(3,natom) :: di, db ! bond vectors
  REAL, DIMENSION(3)       :: zeta   ! random numbers

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, temperature, r_cut, dr_max, de_max

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_nvt_poly_lj'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble, polyatomic molecule'
  CALL introduction ( output_unit )
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
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of molecules',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density

  ! Body-fixed bond vectors in sigma units
  na = natom ! sets the number of atoms per molecule in the module
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Number of atoms per molecule', na
  db(:,1) = [ 0.0, 0.0,  1.0/SQRT(3.0) ]
  db(:,2) = [-0.5, 0.0, -0.5/SQRT(3.0) ]
  db(:,3) = [ 0.5, 0.0, -0.5/SQRT(3.0) ]
  rm_cut = r_cut + 2.0 * SQRT ( MAXVAL ( SUM(db**2,dim=1) ) ) ! centre-centre cutoff distance
  DO i = 1, 3
     WRITE ( unit=output_unit, fmt='(a,i1,t40,3f15.5)' ) 'Body-fixed bond vector ', i, db(:,i)
  END DO
  WRITE( unit=output_unit, fmt='(a,t40,f15.5)') 'Molecule-molecule cutoff distance', rm_cut

  CALL allocate_arrays ( box, rm_cut )

  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! second call is to get r and e

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
  db     = db / box                  ! Convert bond vectors to box units

  DO i = 1, n
     d(:,:,i) = q_to_d ( e(:,i), db ) ! Convert quaternions to space-fixed bond vectors
  END DO

  eng_old = energy ( box, r_cut, rm_cut )
  IF ( eng_old%overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_poly_lj'
  END IF
  pot = eng_old%pot
  vir = eng_old%vir
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'Potential', 'Virial Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           ri(:)   = r(:,i)   ! copy old position
           ei(:)   = e(:,i)   ! copy old quaternion
           di(:,:) = d(:,:,i) ! copy old bond vectors
           eng_old = energy_1 ( ri, di, i, box, r_cut, rm_cut )

           IF ( eng_old%overlap ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_nvt_poly_lj'
           END IF

           CALL RANDOM_NUMBER ( zeta ) ! three uniform random numbers in range (0,1)
           zeta = 2.0*zeta - 1.0       ! now in range (-1,+1)
           ri(:) = ri(:) + zeta * dr_max / box  ! trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )      ! periodic boundary correction
           ei = random_rotate_quaternion ( de_max, ei ) ! trial rotation
           di = q_to_d ( ei, db ) ! new space-fixed bond vectors (in box=1 units)
           eng_new = energy_1 ( ri, di, i, box, r_cut, rm_cut )

           IF ( .NOT. eng_new%overlap ) THEN ! consider non-overlapping configuration
              delta = ( eng_new%pot - eng_old%pot ) / temperature
              IF ( metropolis ( delta ) ) THEN ! accept Metropolis test
                 pot      = pot + eng_new%pot - eng_old%pot ! update potential energy
                 vir      = vir + eng_new%vir - eng_old%vir ! update virial
                 r(:,i)   = ri                              ! update position
                 e(:,i)   = ei                              ! update quaternion
                 d(:,:,i) = di                              ! update bond vectors
                 moves    = moves + 1                       ! increment move counter
              END IF ! reject Metropolis test
           END IF ! reject overlapping configuration

        END DO ! End loop over atoms

        ! Calculate all variables for this step
        move_ratio  = REAL(moves) / REAL(n)
        CALL calculate ( )
        CALL blk_add ( [move_ratio,potential,pres_virial] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk              ! number configuration by block
     CALL write_cnf_mols ( cnf_prefix//sav_tag, n, box, r*box, e ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  eng_old = energy ( box, r_cut, rm_cut )
  IF ( eng_old%overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_poly_lj'
  END IF
  pot = eng_old%pot
  vir = eng_old%vir
  CALL calculate ( 'Final check' )

  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    potential   = pot / REAL ( n )
    pres_virial = density * temperature + vir / box**3

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential energy', potential
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Virial pressure',  pres_virial
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_nvt_poly_lj

