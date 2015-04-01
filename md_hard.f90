! md_hard.f90 (also uses md_hard_module.f90 and io_module.f90)
! Molecular dynamics of hard spheres
PROGRAM md_hard
  USE md_hard_module
  use io_module
  IMPLICIT NONE

  ! Takes in a hard-sphere configuration (positions and velocities)
  ! Checks for overlaps    
  ! Conducts molecular dynamics simulation for specified number of collisions
  ! Uses no special neighbour lists
  ! ... so is restricted to small number of atoms
  ! Assumes that collisions can be predicted by looking at 
  ! nearest neighbour particles in periodic boundaries
  ! ... so is unsuitable for low densities
  ! Box is taken to be of unit length.                       
  ! However, input configuration, output configuration and all
  ! results are given in units sigma = 1, kT = 1, mass = 1

  REAL                        :: sigma   ! atomic diameter (in units where box=1)
  REAL                        :: box     ! box length (in units where sigma=1)
  REAL                        :: density ! reduced density n*sigma**3/box**3
  CHARACTER(len=7), PARAMETER :: prefix = 'md_hard'
  character(len=7), parameter :: cnfinp = '.cnfinp', cnfout = '.cnfout'

  INTEGER   ::  i, j, k, ncoll, coll
  REAL      ::  tij, t, rate, kinetic_energy
  real, dimension(3) :: total_momentum
  REAL      ::  virial, pvnkt1, virial_avg, temperature, tbc, sigma_sq

  NAMELIST /run_parameters/ ncoll

  WRITE(*,'(''md_hard'')')
  WRITE(*,'(''Molecular dynamics of hard spheres'')')
  WRITE(*,'(''Results in units kT = sigma = 1'')')

  READ(*,nml=run_parameters)
  WRITE(*,'(''Collisions required'',t40,i15)'  ) ncoll

  call read_cnf_atoms ( prefix//cnfinp, n, box )
  WRITE(*,'(''Number of particles'',t40,i15)') n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma = 1.0
  WRITE(*,'(''Sigma (in sigma units)'',t40,f15.5)') sigma
  density = real (n) * ( sigma / box ) ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density

  ALLOCATE ( r(3,n), v(3,n), coltime(n), partner(n) )

  call read_cnf_atoms ( prefix//cnfinp, n, box, r, v )
  total_momentum = sum(v,dim=2)
  total_momentum = total_momentum / real(n)
  write(*,'(''Net momentum'',t40,3f15.5)') total_momentum
  v = v - spread(total_momentum,dim=1,ncopies=3)
  kinetic_energy = 0.5 * SUM ( v**2 )
  temperature = 2.0 * kinetic_energy / REAL ( 3*(n-1) )
  WRITE(*,'(''Initial temperature (sigma units)'',t40,f15.5)') temperature

  ! Convert to box units
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
     CALL update ( i, i+1, n, sigma_sq ) ! initial search for collision partners >i
  END DO

  virial_avg = 0.0 ! zero virial accumulator

  WRITE(*,'(''Start of dynamics'')')

  t = 0.0

  DO coll = 1, ncoll ! Start of main loop

     i   = MINLOC ( coltime, dim=1 ) ! locate minimum collision time
     j   = partner(i)                ! collision partner
     tij = coltime(i)                ! time to collision

     t          = t          + tij           ! advance time by tij
     coltime(:) = coltime(:) - tij           ! reduce times to next collision by tij
     r(:,:)     = r(:,:)     + tij * v(:,:)  ! advance all particles by tij
     r(:,:)     = r(:,:) - ANINT ( r(:,:) )  ! apply periodic boundaries

     CALL collide ( i, j, sigma_sq, virial ) ! compute collision dynamics

     virial_avg = virial_avg + virial

     DO k = 1, n
        IF ( ( k == i ) .OR. ( partner(k) == i ) .OR. ( k == j ) .OR. ( partner(k) == j ) ) THEN
           CALL update ( k, k+1, n, sigma_sq ) ! search for partners >k
        ENDIF
     END DO

     CALL update ( i, 1, i-1, sigma_sq ) ! search for partners <i
     CALL update ( j, 1, j-1, sigma_sq ) ! search for partners <j

  END DO ! End of main loop

  WRITE(*,'(''End of dynamics'')')
  WRITE(*,'(''Final colliding pair '',2i5)') i, j

  IF ( overlap ( sigma_sq ) ) THEN
     WRITE(*,'(''Particle overlap in final configuration'')')
  ENDIF

  kinetic_energy = 0.5 * SUM ( v**2 ) ! in box units
  temperature = 2.0 * kinetic_energy / REAL ( 3*(n-1) ) ! in box units
  pvnkt1 = virial_avg / REAL ( n ) / 3.0 / t / temperature ! dimensionless
  t = t * SQRT ( temperature ) / sigma ! in sigma units
  rate = REAL ( ncoll ) / t ! in sigma units
  tbc  = REAL ( n ) / rate / 2.0 ! in sigma units
  temperature = temperature * box**2 ! in sigma units

  WRITE(*,'(''Final temperature (sigma units)'',t40,f15.5)') temperature
  WRITE(*,'(''PV/NkT - 1 is'',t40,f15.8)') pvnkt1
  WRITE(*,'(''Final time (sigma units)'',t40,f15.8)') t
  WRITE(*,'(''Collision rate (sigma units)'',t40,f15.8)') rate
  WRITE(*,'(''Mean collision time (sigma units)'',t40,f15.8)') tbc

! Convert from box units
  r(:,:) = r(:,:) * box
  v(:,:) = v(:,:) * box
  call write_cnf_atoms ( prefix//cnfout, n, box, r, v )

  DEALLOCATE ( r, v, coltime, partner )

END PROGRAM md_hard

