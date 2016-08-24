! md_nve_hs.f90
! Molecular dynamics, NVE ensemble, hard spheres
PROGRAM md_nve_hs

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
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
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in units sigma = 1, mass = 1

  ! Most important variables
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: vir         ! total collisional virial
  REAL :: kin         ! kinetic energy
  REAL :: temperature ! temperature
  REAL :: t           ! time

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  REAL :: coll_rate, pressure ! quantities to be averaged
  
  INTEGER            :: i, j, k, ncoll, coll, blk, nblock, ioerr
  REAL               :: tij, vir_sum
  REAL, DIMENSION(3) :: total_momentum

  NAMELIST /nml/ nblock, ncoll

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_nve_hs'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE, hard spheres'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Results in units sigma = 1, mass = 1'
  CALL time_stamp ( output_unit )

  ! Set sensible defaults for testing
  nblock = 10
  ncoll  = 10000
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_nve_hs'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Number of blocks',               nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Number of collisions per block', ncoll

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Box (in sigma units)', box
  density = REAL (n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Reduced density', density

  CALL allocate_arrays

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v
  total_momentum = SUM(v,dim=2)
  total_momentum = total_momentum / REAL(n)
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Net momentum/particle', total_momentum
  v   = v - SPREAD(total_momentum,dim=1,ncopies=3)
  kin = 0.5 * SUM ( v**2 )
  temperature = 2.0 * kin / REAL ( 3*(n-1) )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temperature (sigma units)', temperature

  r(:,:) = r(:,:) / box              ! Convert positions to box=1 units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  IF ( overlap ( box ) ) THEN
     WRITE ( unit=error_unit, fmt='(a)' ) 'Particle overlap in initial configuration'
     STOP 'Error in md_nve_hs'
  END IF

  coltime(:) = HUGE(1.0)
  partner(:) = n

  DO i = 1, n
     CALL update ( i, gt, box ) ! initial search for collision partners >i
  END DO

  CALL run_begin ( [ CHARACTER(len=15) :: 'Collision Rate', 'Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin
     vir_sum = 0.0
     t       = 0.0

     DO coll = 1, ncoll ! Begin loop over collisions

        i   = MINLOC ( coltime, dim=1 ) ! locate minimum collision time
        j   = partner(i)                ! collision partner
        tij = coltime(i)                ! time to collision

        t          = t + tij                     ! advance time by tij
        coltime(:) = coltime(:) - tij            ! reduce times to next collision by tij
        r(:,:)     = r(:,:) + tij * v(:,:) / box ! advance all positions by tij (box=1 units)
        r(:,:)     = r(:,:) - ANINT ( r(:,:) )   ! apply periodic boundaries

        CALL collide ( i, j, box, vir ) ! compute collision dynamics

        vir_sum = vir_sum + vir

        DO k = 1, n
           IF ( ( k == i ) .OR. ( partner(k) == i ) .OR. ( k == j ) .OR. ( partner(k) == j ) ) THEN
              CALL update ( k, gt, box ) ! search for partners >k
           ENDIF
        END DO

        CALL update ( i, lt, box ) ! search for partners <i
        CALL update ( j, lt, box ) ! search for partners <j

     END DO ! End loop over collisions

      ! Collisional time averages in sigma units
     coll_rate = 2.0*REAL (ncoll) / t / REAL(n)             ! collision rate per particle
     pressure  = density*temperature + vir_sum / t / box**3 ! ideal + collisional virial / volume
     CALL blk_add ( [coll_rate, pressure] ) ! time averages
     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  WRITE ( unit=output_unit, fmt='(a,t40,2i5)' ) 'Final colliding pair', i, j

  IF ( overlap ( box ) ) STOP 'Particle overlap in final configuration'

  r(:,:) = r(:,:) * box ! Convert positions back from box=1 units
  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r, v )
  CALL time_stamp ( output_unit )

  CALL deallocate_arrays

END PROGRAM md_nve_hs

