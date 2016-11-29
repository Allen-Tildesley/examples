! md_nve_hs.f90
! Molecular dynamics, NVE ensemble, hard spheres
PROGRAM md_nve_hs

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add, variable_type
  USE md_module,        ONLY : allocate_arrays, deallocate_arrays, update, overlap, collide, &
       &                       n, r, v, coltime, partner, lt, gt

  IMPLICIT NONE

  ! Takes in a hard-sphere configuration (positions and velocities)
  ! Checks for overlaps    
  ! Conducts molecular dynamics simulation
  ! Uses no special neighbour lists
  ! ... so is restricted to small number of atoms
  ! Assumes that collisions can be predicted by looking at 
  ! nearest neighbour particles in periodic boundaries
  ! ... so is unsuitable for low densities

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are stored divided by the box length
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in units sigma = 1, mass = 1

  ! Most important variables
  REAL :: box        ! Box length (in units where sigma=1)
  REAL :: vir        ! Total collisional virial
  REAL :: kin        ! Kinetic energy
  REAL :: temp_kinet ! Temperature (conserved)
  REAL :: t          ! Time

  ! Quantities to be averaged
  TYPE(variable_type), DIMENSION(:), ALLOCATABLE :: variables

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  INTEGER            :: i, j, k, ncoll, coll, blk, nblock, ioerr
  REAL               :: tij, vir_sum
  REAL, DIMENSION(3) :: vcm

  NAMELIST /nml/ nblock, ncoll

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_nve_hs'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE, hard spheres'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Results in units sigma = 1, mass = 1'
  CALL time_stamp

  ! Set sensible default run parameters for testing
  nblock = 10
  ncoll  = 10000

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_nve_hs'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Number of blocks',               nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Number of collisions per block', ncoll

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Box (in sigma units)', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',              REAL(n) / box**3
  CALL allocate_arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v
  r(:,:) = r(:,:) / box                                     ! Convert positions to box=1 units
  r(:,:) = r(:,:) - ANINT ( r(:,:) )                        ! Periodic boundaries
  vcm(:) = SUM ( v(:,:), dim=2 ) / REAL(n)                  ! Centre-of mass velocity
  v(:,:) = v(:,:) - SPREAD ( vcm(:), dim = 2, ncopies = n ) ! Set COM velocity to zero
  kin        = 0.5 * SUM ( v**2 )
  temp_kinet = 2.0 * kin / REAL ( 3*(n-1) )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature', temp_kinet

  ! Initial overlap check
  IF ( overlap ( box ) ) THEN
     WRITE ( unit=error_unit, fmt='(a)' ) 'Particle overlap in initial configuration'
     STOP 'Error in md_nve_hs'
  END IF

  ALLOCATE ( variables(2) )

  ! Initial search for collision partners >i
  coltime(:) = HUGE(1.0)
  partner(:) = n
  DO i = 1, n
     CALL update ( i, box, gt ) 
  END DO

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( variables )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin
     vir_sum = 0.0
     t       = 0.0

     DO coll = 1, ncoll ! Begin loop over collisions

        i   = MINLOC ( coltime, dim=1 ) ! Locate minimum collision time
        j   = partner(i)                ! Collision partner
        tij = coltime(i)                ! Time to collision

        t          = t + tij                     ! Advance time by tij
        coltime(:) = coltime(:) - tij            ! Reduce times to next collision by tij
        r(:,:)     = r(:,:) + tij * v(:,:) / box ! Advance all positions by tij (box=1 units)
        r(:,:)     = r(:,:) - ANINT ( r(:,:) )   ! Apply periodic boundaries

        CALL collide ( i, j, box, vir ) ! Compute collision dynamics

        vir_sum = vir_sum + vir

        DO k = 1, n
           IF ( ( k == i ) .OR. ( partner(k) == i ) .OR. ( k == j ) .OR. ( partner(k) == j ) ) THEN
              CALL update ( k, box, gt ) ! Search for partners >k
           ENDIF
        END DO

        CALL update ( i, box, lt ) ! search for partners <i
        CALL update ( j, box, lt ) ! search for partners <j

     END DO ! End loop over collisions

     ! Calculate and accumulate variables for this step
     CALL calculate
     CALL blk_add ( variables )
     CALL blk_end ( blk )                                           ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ! Output run averages

  WRITE ( unit=output_unit, fmt='(a,t40,2i5)' ) 'Final colliding pair', i, j

  IF ( overlap ( box ) ) STOP 'Particle overlap in final configuration'

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )
  CALL time_stamp

  CALL deallocate_arrays

CONTAINS

  SUBROUTINE calculate
    IMPLICIT NONE

    ! This routine calculates all variables of interest
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: coll_rate, p_coll
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

    ! Collision rate per particle
    coll_rate = variable_type ( nam = 'Collision rate', val = 2.0*REAL(ncoll)/t/REAL(n) )

    ! Collisional pressure
    ! ideal + collisional virial / volume
    p_coll = variable_type ( nam = 'P', val = rho*temp_kinet + vir_sum/t/vol )

    ! Collect together for averaging
    ! Fortran 2003 should automatically allocate this first time
    variables = [ coll_rate, p_coll ]
  END SUBROUTINE calculate

END PROGRAM md_nve_hs

