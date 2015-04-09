! md_nve_hs.f90 (uses md_nve_hs_module.f90, utility_module.f90)
! Molecular dynamics, NVE ensemble, hard spheres
PROGRAM md_nve_hs
  USE utility_module, ONLY : read_cnf_atoms, write_cnf_atoms, &
       &                     run_begin, run_end, blk_begin, blk_end, blk_add
  USE md_nve_hs_module, ONLY : update, overlap, collide, n, r, v, coltime, partner, lt, gt
  IMPLICIT NONE

  ! Takes in a hard-sphere configuration (positions and velocities)
  ! Checks for overlaps    
  ! Conducts molecular dynamics simulation
  ! Uses no special neighbour lists
  ! ... so is restricted to small number of atoms
  ! Assumes that collisions can be predicted by looking at 
  ! nearest neighbour particles in periodic boundaries
  ! ... so is unsuitable for low densities
  ! Box is taken to be of unit length during the dynamics.
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in units sigma = 1, mass = 1

  ! Most important variables
  REAL :: sigma       ! atomic diameter (in units where box=1)
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: vir         ! total collisional virial
  real :: kin         ! kinetic energy
  real :: temperature ! temperature
  real :: t           ! time

  CHARACTER(len=11), PARAMETER :: cnf_prefix = 'md_hard.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number

  REAL :: coll_rate, pressure ! quantities to be averaged
  
  INTEGER            :: i, j, k, ncoll, coll, blk, nblock
  REAL               :: tij
  REAL               :: vir_sum, sigma_sq
  REAL, DIMENSION(3) :: total_momentum

  NAMELIST /run_parameters/ nblock, ncoll

  WRITE(*,'(''md_hard'')')
  WRITE(*,'(''Molecular dynamics of hard spheres'')')
  WRITE(*,'(''Results in units sigma = 1, mass = 1'')')

  ! Set sensible defaults for testing
  nblock = 10
  ncoll  = 10000
  READ(*,nml=run_parameters)
  WRITE(*,'(''Number of blocks'',              t40,i15)') nblock
  WRITE(*,'(''Number of collisions per block'',t40,i15)') ncoll

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box )
  WRITE(*,'(''Number of particles'', t40,i15  )') n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma = 1.0
  density = REAL (n) * ( sigma / box ) ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density

  ALLOCATE ( r(3,n), v(3,n), coltime(n), partner(n) )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v )
  total_momentum = SUM(v,dim=2)
  total_momentum = total_momentum / REAL(n)
  WRITE(*,'(''Net momentum'',t40,3f15.5)') total_momentum
  v   = v - SPREAD(total_momentum,dim=1,ncopies=3)
  kin = 0.5 * SUM ( v**2 )
  temperature = 2.0 * kin / REAL ( 3*(n-1) )
  WRITE(*,'(''Temperature (sigma units)'',t40,f15.5)') temperature

  ! Convert to box units (time units are unaffected)
  r(:,:) = r(:,:) / box
  v(:,:) = v(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries
  sigma = sigma / box
  sigma_sq = sigma ** 2
  WRITE(*,'(''Sigma (in box units)'',t40,f15.5)'  ) sigma

  IF ( overlap ( sigma_sq ) ) THEN
     STOP 'particle overlap in initial configuration'
  END IF

  coltime(:) = HUGE(1.0)
  partner(:) = n

  DO i = 1, n
     CALL update ( i, gt, sigma_sq ) ! initial search for collision partners >i
  END DO

  CALL run_begin ( ['Coll Rate ','Pressure  '] ) ! must all be character*10 constants

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin
     vir_sum = 0.0
     t = 0.0

     DO coll = 1, ncoll ! Begin loop over collisions

        i   = MINLOC ( coltime, dim=1 ) ! locate minimum collision time
        j   = partner(i)                ! collision partner
        tij = coltime(i)                ! time to collision

        t          = t          + tij           ! advance time by tij
        coltime(:) = coltime(:) - tij           ! reduce times to next collision by tij
        r(:,:)     = r(:,:)     + tij * v(:,:)  ! advance all particles by tij
        r(:,:)     = r(:,:) - ANINT ( r(:,:) )  ! apply periodic boundaries

        CALL collide ( i, j, sigma_sq, vir ) ! compute collision dynamics

        vir_sum = vir_sum + vir

        DO k = 1, n
           IF ( ( k == i ) .OR. ( partner(k) == i ) .OR. ( k == j ) .OR. ( partner(k) == j ) ) THEN
              CALL update ( k, gt, sigma_sq ) ! search for partners >k
           ENDIF
        END DO

        CALL update ( i, lt, sigma_sq ) ! search for partners <i
        CALL update ( j, lt, sigma_sq ) ! search for partners <j

     END DO ! End loop over collisions

      ! Collisional time averages in sigma units
     coll_rate = 2.0*REAL (ncoll) / t / real(n)             ! collision rate per particle
     pressure  = density*temperature + vir_sum / t / box**3 ! ideal + collisional
     CALL blk_add ( [coll_rate, pressure] ) ! time averages
     CALL blk_end ( blk )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v*box ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end

  WRITE(*,'(''Final colliding pair'',t40,2i5)') i, j

  IF ( overlap ( sigma_sq ) ) STOP 'Particle overlap in final configuration'

  ! Convert from box units
  r(:,:) = r(:,:) * box
  v(:,:) = v(:,:) * box
  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r, v )

  DEALLOCATE ( r, v, coltime, partner )

END PROGRAM md_nve_hs

