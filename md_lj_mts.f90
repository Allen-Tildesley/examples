! md_lj_mts.f90
! Molecular dynamics, NVE ensemble, Lennard-Jones atoms
PROGRAM md_lj_mts
  USE utility_module,   ONLY : read_cnf_atoms, write_cnf_atoms, time_stamp, &
       &                       run_begin, run_end, blk_begin, blk_end, blk_add
  USE md_lj_mts_module, ONLY : allocate_arrays, deallocate_arrays, force, r, v, f, n
  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using velocity Verlet algorithm
  ! Uses no special neighbour lists, for clarity
  
  ! This program just illustrates the idea of splitting the non-bonded interactions
  ! using a criterion based on distance, for use in a MTS scheme
  ! This would hardly ever be efficient for a simple LJ potential alone

  ! This program uses LJ units sigma = 1, epsilon = 1, mass = 1 throughout

  ! Most important variables
  REAL :: box         ! box length (in units where sigma=1)
  REAL :: density     ! reduced density n*sigma**3/box**3
  REAL :: dt          ! time step (smallest)
  REAL :: kin         ! total kinetic energy
  REAL :: pressure    ! pressure (LJ sigma=1 units, to be averaged)
  REAL :: temperature ! temperature (LJ sigma=1 units, to be averaged)
  REAL :: energy      ! total energy per atom (LJ sigma=1 units, to be averaged)
  REAL :: lambda      ! healing length for switch function

  INTEGER, PARAMETER        :: k_max = 3   ! number of shells
  REAL,    DIMENSION(k_max) :: r_cut       ! cutoff distance for each shell
  REAL,    DIMENSION(k_max) :: pot         ! total potential energy for each shell
  REAL,    DIMENSION(k_max) :: vir         ! total virial for each shell
  INTEGER, DIMENSION(k_max) :: n_mts       ! successive ratios of number of steps for each shell

  INTEGER :: blk, stp1, stp2, stp3, nstep, nblock, k
  REAL    :: pairs

  CHARACTER(len=13), PARAMETER :: cnf_prefix = 'md_nve_lj.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number

  REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

  NAMELIST /params/ nblock, nstep, r_cut, lambda, dt, n_mts

  WRITE(*,'(''md_nve_lj'')')
  WRITE(*,'(''Molecular dynamics, constant-NVE, Lennard-Jones'')')
  WRITE(*,'(''Results in units epsilon = sigma = 1'')')
  CALL time_stamp

  ! Set sensible default run parameters for testing
  nblock = 10
  nstep  = 1000
  r_cut  = [ 2.4, 3.5, 4.0 ]
  n_mts  = [ 1, 4, 2 ]
  dt     = 0.002
  lambda = 0.1

  READ(*,nml=params)
  WRITE(*,'(''Number of blocks'',          t40,  i15)'   ) nblock
  WRITE(*,'(''Number of steps per block'', t40,  i15)'   ) nstep
  WRITE(*,'(''Potential cutoff distances'',t40,*(f15.5))') r_cut

  DO k = 1, k_max
     IF ( k == 1 ) THEN
        pairs = r_cut(k)**3
     ELSE
        pairs = r_cut(k)**3 - r_cut(k-1)**3
        IF ( r_cut(k)-r_cut(k-1) < lambda ) STOP 'r_cut values must differ by at least lambda'
     END IF
     pairs = REAL(n*(n-1)/2) * (4.0/3.0)*pi * pairs / box**3
     WRITE(*,'(a,i1,t40,i15)') 'Estimated pairs in shell ', k, NINT ( pairs )
  END DO

  WRITE(*,'(''Multiple step ratios'',t40,*(i15))'  ) n_mts(:)
  IF ( n_mts(1) /= 1 ) STOP 'n_mts(1) must be 1'
  IF ( ANY ( n_mts <= 0 ) ) STOP 'n_mts values must be positive'
  DO k = 1, k_max
     WRITE(*,'(''Time step for shell '', i1, t40,  f15.5)' ) k, PRODUCT(n_mts(1:k))*dt
  END DO

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box )
  WRITE(*,'(''Number of particles'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  density = REAL(n) / box ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density
  IF ( r_cut(k_max) > box/2.0  ) STOP 'r_cut(k_max) too large '

  CALL allocate_arrays ( k_max )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v )

  r(:,:) = r(:,:) - ANINT ( r(:,:) / box ) * box ! Periodic boundaries

  ! Calculate forces and pot, vir contributions for each shell
  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, pot(k), vir(k) )
  END DO
  kin         = 0.5*SUM(v**2)
  energy      = ( SUM(pot) + kin ) / REAL ( n )
  temperature = 2.0 * kin / REAL ( 3*(n-1) )
  pressure    = density * temperature + SUM(vir) / box**3
  WRITE(*,'(''Initial total energy (sigma units)'',  t40,f15.5)') energy
  WRITE(*,'(''Initial temperature (sigma units)'',   t40,f15.5)') temperature
  WRITE(*,'(''Initial pressure (sigma units)'',      t40,f15.5)') pressure

  CALL run_begin ( [ CHARACTER(len=15) :: 'Energy', 'Temperature', 'Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     ! The following set of nested loops is specific to k_max=3
     
     DO stp3 = 1, nstep ! Begin loop over steps

        ! Outer shell 3: a single step of size n_mts(3) * n_mts(2) * dt
        v(:,:) = v(:,:) + 0.5 * n_mts(3) * n_mts(2) * dt * f(:,:,3) ! Kick half-step (outer shell)

        DO stp2 = 1, n_mts(3) ! Middle shell 2: n_mts(3) steps of size n_mts(2) * dt

           v(:,:) = v(:,:) + 0.5 * n_mts(2) * dt * f(:,:,2) ! Kick half-step (middle shell)

           DO stp1 = 1, n_mts(2) ! Inner shell 1: n_mts(3) * n_mts(2) steps of size dt

              v(:,:) = v(:,:) + 0.5 * dt * f(:,:,1)                ! Kick half-step (inner shell)
              r(:,:) = r(:,:) + dt * v(:,:)                        ! Drift step
              r(:,:) = r(:,:) - ANINT ( r(:,:)/box ) * box         ! Periodic boundaries
              CALL force ( box, r_cut, lambda, 1, pot(1), vir(1) ) ! Force evaluation (inner shell)
              v(:,:) = v(:,:) + 0.5 * dt * f(:,:,1)                ! Kick half-step (inner shell)

           END DO ! End inner shell 1

           CALL force ( box, r_cut, lambda, 2, pot(2), vir(2) ) ! Force evaluation (middle shell)
           v(:,:) = v(:,:) + 0.5 * n_mts(2) * dt * f(:,:,2)     ! Kick half-step (middle shell)

        END DO ! End middle shell 2

        CALL force ( box, r_cut, lambda, 3, pot(3), vir(3) )        ! Force evaluation (outer shell)
        v(:,:) = v(:,:) + 0.5 * n_mts(3) * n_mts(2) * dt * f(:,:,3) ! Kick half-step (outer shell)
        ! End outer shell 3

        kin         = 0.5*SUM(v**2)
        energy      = ( SUM(pot) + kin ) / REAL ( n )
        temperature = 2.0 * kin / REAL ( 3*(n-1) )
        pressure    = density * temperature + SUM(vir) / box**3

        ! Calculate all variables for this step
        CALL blk_add ( [energy,temperature,pressure] )

     END DO ! End loop over steps

     CALL blk_end ( blk )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk           ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end

  DO k = 1, k_max
     CALL force ( box, r_cut, lambda, k, pot(k), vir(k) )
  END DO
  kin         = 0.5*SUM(v**2)
  energy      = ( SUM(pot) + kin ) / REAL ( n )
  temperature = 2.0 * kin / REAL ( 3*(n-1) )
  pressure    = density * temperature + SUM(vir) / box**3
  WRITE(*,'(''Final total energy (sigma units)'',  t40,f15.5)') energy
  WRITE(*,'(''Final temperature (sigma units)'',   t40,f15.5)') temperature
  WRITE(*,'(''Final pressure (sigma units)'',      t40,f15.5)') pressure
  CALL time_stamp

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r, v )

  CALL deallocate_arrays

END PROGRAM md_lj_mts

