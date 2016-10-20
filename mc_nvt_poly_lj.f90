! mc_nvt_poly_lj.f90
! Monte Carlo, NVT ensemble, polyatomic molecule
PROGRAM mc_nvt_poly_lj

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_mols, write_cnf_mols
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE maths_module,     ONLY : metropolis, random_rotate_quaternion, random_translate_vector
  USE mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       potential_1, potential, q_to_d, n, na, r, e, d, db, potential_type

  IMPLICIT NONE

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
  ! The model is defined in mc_module

  ! Most important variables
  REAL :: box         ! box length
  REAL :: density     ! density
  REAL :: dr_max      ! maximum MC translational displacement
  REAL :: de_max      ! maximum MC rotational displacement
  REAL :: temperature ! specified temperature
  REAL :: r_cut       ! potential cutoff distance
  REAL :: pot         ! total potential energy
  REAL :: vir         ! total virial

  ! Quantities to be averaged
  REAL :: move_ratio ! acceptance ratio of moves
  REAL :: p_cutsh    ! pressure (cut & shifted potential)
  REAL :: en_cutsh   ! internal energy per molecule (cut & shifted potential)

  ! Composite interaction = pot & vir & overlap variables
  TYPE(potential_type) :: system, molecule_old, molecule_new

  INTEGER :: blk, stp, i, nstep, nblock, moves, ioerr
  REAL    :: delta

  REAL, DIMENSION(3)    :: ri
  REAL, DIMENSION(0:3)  :: ei
  REAL, DIMENSION(3,na) :: di

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

  CALL allocate_arrays ( box, r_cut )

  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! second call is to get r and e

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  DO i = 1, n ! Loop over molecules
     d(:,:,i) = q_to_d ( e(:,i), db ) ! Convert quaternions to space-fixed atom vectors
  END DO ! End loop over molecules

  system = potential ( box, r_cut )
  IF ( system%overlap ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in mc_nvt_poly_lj'
  END IF
  pot = system%pot
  vir = system%vir
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Move ratio', 'E/N (cut&shift)', 'P (cut&shift)' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        moves = 0

        DO i = 1, n ! Begin loop over atoms

           molecule_old = potential_1 ( r(:,i), d(:,:,i), i, box, r_cut ) ! Old molecule potential, virial, etc

           IF ( molecule_old%overlap ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_nvt_poly_lj'
           END IF

           ri = random_translate_vector ( dr_max/box, r(:,i) ) ! Trial move to new position (in box=1 units)
           ri = ri - ANINT ( ri )                              ! Periodic boundary correction
           ei = random_rotate_quaternion ( de_max, e(:,i) )    ! Trial rotation
           di = q_to_d ( ei, db )                              ! New space-fixed atom vectors

           molecule_new = potential_1 ( ri, di, i, box, r_cut ) ! New molecule potential, virial etc

           IF ( .NOT. molecule_new%overlap ) THEN ! Test for non-overlapping configuration

              delta = ( molecule_new%pot - molecule_old%pot ) / temperature

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 pot      = pot + molecule_new%pot - molecule_old%pot ! Update potential energy
                 vir      = vir + molecule_new%vir - molecule_old%vir ! Update virial
                 r(:,i)   = ri                                        ! Update position
                 e(:,i)   = ei                                        ! Update quaternion
                 d(:,:,i) = di                                        ! Update atom vectors
                 moves    = moves + 1                                 ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End test for non-overlapping configuration

        END DO ! End loop over atoms

        move_ratio  = REAL(moves) / REAL(n)

        ! Calculate all variables for this step
        CALL calculate ( )
        CALL blk_add ( [move_ratio,en_cutsh,p_cutsh] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk              ! number configuration by block
     CALL write_cnf_mols ( cnf_prefix//sav_tag, n, box, r*box, e ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL calculate ( 'Final values' )

  system = potential ( box, r_cut )
  IF ( system%overlap ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in mc_nvt_poly_lj'
  END IF
  pot = system%pot
  vir = system%vir
  CALL calculate ( 'Final check' )

  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE calculate ( string )
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: string

    en_cutsh = pot / REAL ( n )                ! PE for cut & shifted potential
    en_cutsh = en_cutsh + 3.0 * temperature    ! Add ideal gas contribution KE/N assuming nonlinear molecules 
    p_cutsh  = vir / box**3                    ! Virial contribution to P
    p_cutsh  = p_cutsh + density * temperature ! Add ideal gas contribution to P

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (cut&shift)', en_cutsh
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (cut&shift)',   p_cutsh
    END IF

  END SUBROUTINE calculate

END PROGRAM mc_nvt_poly_lj

