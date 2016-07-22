! mc_nvt_poly_lj.f90
! Monte Carlo, NVT ensemble, polyatomic molecule, LJ atoms
PROGRAM mc_nvt_poly_lj
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit
  USE utility_module,    ONLY : metropolis, read_cnf_mols, write_cnf_mols, random_rotate_quaternion, &
       &                        run_begin, run_end, blk_begin, blk_end, blk_add
  USE mc_poly_lj_module, ONLY : allocate_arrays, deallocate_arrays, energy_1, energy, q_to_d, &
       &                        n, na, r, e, d, ne
  IMPLICIT NONE

  ! Takes in a configuration of polyatomic molecules (positions and quaternions)
  ! Cubic periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature
  ! Uses no special neighbour lists

  ! Box is taken to be of unit length during the Monte Carlo
  ! However, input configuration, output configuration,
  ! most calculations, and all results 
  ! are given in LJ units sigma = 1, epsilon = 1

  ! Most important variables
  REAL :: sigma       ! atomic diameter (in units where box=1)
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
  INTEGER            :: blk, stp, i, nstep, nblock, moves
  REAL               :: pot_old, pot_new, vir_old, vir_new, delta

  REAL, DIMENSION(3)       :: ri   ! position of atom i
  REAL, DIMENSION(0:3)     :: ei   ! orientation of atom i
  INTEGER, PARAMETER       :: natom = 3
  REAL, DIMENSION(3,natom) :: di, db ! bond vectors
  REAL, DIMENSION(3)       :: zeta ! random numbers

  CHARACTER(len=18), PARAMETER :: cnf_prefix = 'mc_nvt_poly_lj.cnf'
  CHARACTER(len=3),  PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)             :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /params/ nblock, nstep, temperature, r_cut, dr_max, de_max

  WRITE(*,'(''mc_nvt_poly_lj'')')
  WRITE(*,'(''Monte Carlo, constant-NVT, polyatomic, Lennard-Jones'')')
  WRITE(*,'(''Results in units epsilon = sigma = 1'')')

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible defaults for testing
  nblock      = 10
  nstep       = 1000
  temperature = 0.7
  r_cut       = 2.5
  dr_max      = 0.15
  de_max      = 0.1
  READ(*,nml=params)
  WRITE(*,'(''Number of blocks'',         t40,i15)'  ) nblock
  WRITE(*,'(''Number of steps per block'',t40,i15)'  ) nstep
  WRITE(*,'(''Temperature'',              t40,f15.5)') temperature
  WRITE(*,'(''Potential cutoff distance'',t40,f15.5)') r_cut
  WRITE(*,'(''Maximum r displacement'',   t40,f15.5)') dr_max
  WRITE(*,'(''Maximum e displacement'',   t40,f15.5)') de_max

  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box )
  WRITE(*,'(''Number of molecules'', t40,i15)'  ) n
  WRITE(*,'(''Box (in sigma units)'',t40,f15.5)') box
  sigma = 1.0
  density = REAL(n) * ( sigma / box ) ** 3
  WRITE(*,'(''Reduced density'',t40,f15.5)') density

  ! Body-fixed bond vectors in sigma units
  na = natom
  WRITE(*,'(''Number of atoms per molecule'',t40,i15)') na
  db(:,1) = [ 0.0, 0.0,  1.0/SQRT(3.0) ]
  db(:,2) = [-0.5, 0.0, -0.5/SQRT(3.0) ]
  db(:,3) = [ 0.5, 0.0, -0.5/SQRT(3.0) ]
  rm_cut = r_cut + 2.0 * SQRT ( MAXVAL ( SUM(db**2,dim=1) ) )
  DO i = 1, 3
     WRITE(*,'(''Body-fixed bond vector '',i1,t40,3f15.5)') i, db(:,i)
  END DO
  WRITE(*,'(''Molecule-molecule cutoff distance'',t40,f15.5)') rm_cut

  ! Convert run and potential parameters to box units
  sigma  = sigma / box
  r_cut  = r_cut / box
  rm_cut = rm_cut / box
  dr_max = dr_max / box
  db     = db / box
  WRITE(*,'(''sigma (in box units)'', t40,f15.5)') sigma
  WRITE(*,'(''rm_cut (in box units)'',t40,f15.5)') rm_cut
  IF ( rm_cut > 0.5 ) STOP 'rm_cut too large '

  CALL allocate_arrays

  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e )

  ! Convert to box units
  r(:,:) = r(:,:) / box
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Convert quaternions to space-fixed bond vectors
  DO i = 1, n
     d(:,:,i) = q_to_d ( e(:,i), db )
  END DO

  CALL energy ( sigma, r_cut, rm_cut, overlap, pot, vir )
  IF ( overlap ) STOP 'Overlap in initial configuration'
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(*,'(''Initial potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Initial pressure (sigma units)'',        t40,f15.5)') pressure

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
           CALL  energy_1 ( ri, di, i, ne, sigma, r_cut, rm_cut, overlap, pot_old, vir_old )
           IF ( overlap ) STOP 'Overlap in current configuration'
           ri(:) = ri(:) + zeta * dr_max   ! trial move to new position
           ri(:) = ri(:) - ANINT ( ri(:) ) ! periodic boundary correction
           ei = random_rotate_quaternion ( de_max, ei ) ! trial rotation
           di = q_to_d ( ei, db ) ! new space-fixed bond vectors
           CALL  energy_1 ( ri, di, i, ne, sigma, r_cut, rm_cut, overlap, pot_new, vir_new )

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
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL energy ( sigma, r_cut, rm_cut, overlap, pot, vir )
  IF ( overlap ) STOP 'Overlap in final configuration'
  potential = pot / REAL ( n )
  pressure  = density * temperature + vir / box**3
  WRITE(*,'(''Final check'')')
  WRITE(*,'(''Final potential energy (sigma units)'',t40,f15.5)') potential
  WRITE(*,'(''Final pressure (sigma units)'',        t40,f15.5)') pressure

  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e )

  CALL deallocate_arrays

END PROGRAM mc_nvt_poly_lj

