! md_npt_lj.f90
! Molecular dynamics, NPT ensemble
PROGRAM md_npt_lj
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE maths_module,     ONLY : random_normal, random_normals
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE md_module,        ONLY : model_description, allocate_arrays, deallocate_arrays, &
       &                       force, r, v, f, n, energy_lrc

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using a measure-preserving algorithm for NPT
  ! Nose-Hoover chains are used, following Martyna et al, Molec Phys, 87, 1117 (1996)
  ! and Tuckerman et al J Phys A, 39, 5629 (2006)
  ! To keep this example reasonably simple, we do not subdivide the timesteps with a
  ! Suzuki-Yoshida decomposition, as described in those papers
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  ! Most important variables
  REAL :: box         ! box length
  REAL :: dt          ! time step
  REAL :: r_cut       ! potential cutoff distance
  REAL :: temperature ! specified temperature
  REAL :: pressure    ! specified pressure
  REAL :: g           ! number of degrees of freedom
  REAL :: conserved   ! conserved energy-like quantity
  REAL :: tau         ! thermostat time scale
  REAL :: tau_baro    ! barostat time scale
  REAL :: w_eps       ! barostat inertia
  REAL :: eps         ! strain
  REAL :: p_eps       ! strain momentum 
  REAL :: box0        ! initial box length 
  REAL :: pot         ! total potential energy
  REAL :: pot_sh      ! total shifted potential energy
  REAL :: kin         ! total kinetic energy
  REAL :: vir         ! total virial
  REAL :: lap         ! total Laplacian
  REAL :: density     ! density (to be averaged)
  REAL :: pres_virial ! virial pressure (to be averaged)
  REAL :: temp_kinet  ! kinetic temperature (to be averaged)
  REAL :: temp_config ! configurational temperature (to be averaged)
  REAL :: energy      ! total energy per atom (to be averaged)
  REAL :: energy_sh   ! total shifted energy per atom (to be averaged)

  INTEGER, PARAMETER    :: m = 3      ! number of Nose-Hoover chain variables
  REAL,    DIMENSION(m) :: q          ! thermal inertias
  REAL,    DIMENSION(m) :: eta        ! thermal coordinates (needed only to calculate conserved quantity)
  REAL,    DIMENSION(m) :: p_eta      ! thermal momenta
  REAL,    DIMENSION(m) :: q_baro     ! barostat thermal inertias
  REAL,    DIMENSION(m) :: eta_baro   ! barostat thermal coordinates
  REAL,    DIMENSION(m) :: p_eta_baro ! barostat thermal momenta

  INTEGER :: blk, stp, nstep, nblock, ioerr
  REAL    :: pot_lrc, vir_lrc

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, r_cut, dt, temperature, pressure, tau, tau_baro

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_npt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NPT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL model_description ( output_unit )
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  r_cut       = 2.5
  dt          = 0.002
  temperature = 0.7 ! specified temperature
  pressure    = 0.5 ! specified pressure
  tau         = 2.0 ! desired thermostat timescale
  tau_baro    = 2.0 ! desired barostat timescale

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_npt_lj'
  END IF
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Specified temperature',     temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Specified pressure',        pressure
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Thermostat timescale',      tau
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Barostat timescale',        tau_baro

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density', density

  CALL allocate_arrays ( box, r_cut )

  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v

  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial values of thermal variables
  g    = REAL ( 3*(n-1) )
  q    = temperature * tau**2   
  q(1) = g * temperature * tau**2
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.5))' ) 'Thermal inertias Q', q
  eta(:) = 0.0
  CALL random_normals ( 0.0, SQRT(temperature), p_eta(:)  )
  p_eta(:) = p_eta(:) * SQRT(q(:))
  q_baro = temperature * tau_baro**2   
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.5))' ) "Barostat thermal inertias Q'", q_baro
  eta_baro(:) = 0.0
  CALL random_normals ( 0.0, SQRT(temperature), p_eta_baro(:)  )
  p_eta_baro(:) = p_eta_baro(:) * SQRT(q_baro(:))
  w_eps = g * temperature * tau_baro**2
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) "Barostat inertia W", w_eps
  box0  = box ! reference box length for strain
  eps   = 1.0 ! initial strain
  p_eps = random_normal ( 0.0, SQRT(temperature*w_eps) ) ! strain momentum

  CALL force ( box, r_cut, pot, pot_sh, vir, lap )
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'Conserved', 'Energy', 'Shifted Energy', &
       &                                  'Temp-kinet', 'Temp-config', 'Density', 'Virial Pressure' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL u4 ( dt/4.0, m, 1 )
        CALL u3 ( dt/2.0 )
        CALL u4 ( dt/4.0, 1, m )

        CALL u2p ( dt/2.0 )
        CALL u2 ( dt/2.0 )
        CALL u1 ( dt )
        CALL force ( box, r_cut, pot, pot_sh, vir, lap )
        CALL u2 ( dt/2.0 )
        CALL u2p ( dt/2.0 )

        CALL u4 ( dt/4.0, m, 1 )
        CALL u3 ( dt/2.0 )
        CALL u4 ( dt/4.0, 1, m )

        CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
        CALL calculate ()

        ! Calculate all variables for this step
        CALL blk_add ( [energy,energy_sh,temp_kinet,pres_virial] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL force ( box, r_cut, pot, pot_sh, vir, lap )
  CALL energy_lrc ( n, box, r_cut, pot_lrc, vir_lrc )
  CALL calculate ( 'Final values' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )

  CALL deallocate_arrays

CONTAINS

  SUBROUTINE u1 ( t ) ! U1 and U1' combined: position and strain drift propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! time over which to propagate (typically dt)

    REAL            :: x, c
    REAL, PARAMETER :: c1 = -1.0/2.0, c2 = 1.0/6.0, c3 = -1.0/24.0

    ! U1 part
    x = t * p_eps / w_eps ! time step * time derivative of strain
    IF ( x < 0.001 ) THEN ! guard against small values
       c = 1.0 + x * ( c1 + x * ( c2 + x * c3 ) ) ! Taylor series to order 3
    ELSE
       c = (1.0-EXP(-x))/x
    END IF ! end guard against small values

    r(:,:) = r(:,:) + c * t * v(:,:) / box ! positions in box=1 units
    r(:,:) = r(:,:) - ANINT ( r(:,:) )     ! periodic boundaries

    ! U1' part
    ! because we divide by box above, it is important to update eps afterwards
    ! if we did not use box-scaled coordinates this would not matter

    eps     = eps + x          ! update strain
    box     = box0 * EXP(eps)  ! update box length
    density = REAL(n) / box**2 ! update density

  END SUBROUTINE u1

  SUBROUTINE u2 ( t ) ! U2: velocity kick step propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! time over which to propagate (typically dt/2)

    REAL            :: x, c, alpha
    REAL, PARAMETER :: c1 = -1.0/2.0, c2 = 1.0/6.0, c3 = -1.0/24.0

    alpha = 1.0 + 3.0 / g
    x = t * alpha * p_eps / w_eps

    IF ( x < 0.001 ) THEN ! guard against small values
       c = 1.0 + x * ( c1 + x * ( c2 + x * c3 ) ) ! Taylor series to order 3
    ELSE
       c = (1.0-EXP(-x))/x
    END IF ! end guard against small values

    v(:,:) = v(:,:)*EXP(-x) + c * t * f(:,:)

  END SUBROUTINE u2

  SUBROUTINE u2p ( t ) ! U2': strain momentum propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! time over which to propagate (typically dt/2)

    REAL :: alpha, pv

    alpha = 1.0 + 3.0 / g
    pv    = alpha * SUM(v**2) / 3.0 + vir
    p_eps = p_eps + ( pv - pressure*box**3 ) * t

  END SUBROUTINE u2p

  SUBROUTINE u3 ( t ) ! U3 and U3' combined: thermostat propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! time over which to propagate (typically dt/2)

    ! U3 part
    v(:,:) = v(:,:) * EXP ( -t * p_eta(1) / q(1) )
    eta(:) = eta(:) + t * p_eta(:) / q(:)

    ! U3' part
    eta_baro(:) = eta_baro(:) + t * p_eta_baro(:) / q_baro(:)
    p_eps       = p_eps * EXP ( -t * p_eta_baro(1) / q_baro(1) )

  END SUBROUTINE u3

  SUBROUTINE u4 ( t, j_start, j_stop ) ! U4 and U4' combined: thermostat propagator
    IMPLICIT NONE
    REAL,    INTENT(in) :: t               ! time over which to propagate (typically dt/4)
    INTEGER, INTENT(in) :: j_start, j_stop ! order in which to tackle variables

    INTEGER         :: j, j_stride
    REAL            :: gj, x, c
    REAL, PARAMETER :: c1 = -1.0/2.0, c2 = 1.0/6.0, c3 = -1.0/24.0

    IF ( j_start > j_stop ) THEN
       j_stride = -1
    ELSE
       j_stride = 1
    END IF

    ! U4 part

    DO j = j_start, j_stop, j_stride ! loop over each momentum in turn

       IF ( j == 1 ) THEN ! the driver Gj for p_eta(1) is different
          gj = SUM(v**2) - g*temperature
       ELSE
          gj = ( p_eta(j-1)**2 / q(j-1) ) - temperature
       END IF

       IF ( j == m ) THEN ! the equation for p_eta(M) is different
          p_eta(j)  = p_eta(j) + t * gj
       ELSE
          x = t * p_eta(j+1)/q(j+1)
          IF ( x < 0.001 ) THEN ! guard against small values
             c = 1.0 + x * ( c1 + x * ( c2 + x * c3 ) ) ! Taylor series to order 3
          ELSE
             c = (1.0-EXP(-x))/x
          END IF ! end guard against small values
          p_eta(j) = p_eta(j)*EXP(-x) + t * gj * c
       END IF

    END DO ! end loop over each momentum in turn

    ! U4' part

    DO j = j_start, j_stop, j_stride ! loop over each momentum in turn

       IF ( j == 1 ) THEN ! the driver Gj for p_eta_baro(1) is different
          gj = p_eps**2/w_eps - temperature
       ELSE
          gj = ( p_eta_baro(j-1)**2 / q_baro(j-1) ) - temperature
       END IF

       IF ( j == m ) THEN ! the equation for p_eta_baro(M) is different
          p_eta_baro(j)  = p_eta_baro(j) + t * gj
       ELSE
          x = t * p_eta_baro(j+1)/q_baro(j+1)
          IF ( x < 0.001 ) THEN ! guard against small values
             c = 1.0 + x * ( c1 + x * ( c2 + x * c3 ) ) ! Taylor series to order 3
          ELSE
             c = (1.0-EXP(-x))/x
          END IF ! end guard against small values
          p_eta_baro(j) = p_eta_baro(j)*EXP(-x) + t * gj * c
       END IF

    END DO ! end loop over each momentum in turn

  END SUBROUTINE u4

  SUBROUTINE calculate ( string ) 
    IMPLICIT NONE
    CHARACTER (len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates variables of interest and (optionally) writes them out

    pot         = pot + pot_lrc
    vir         = vir + vir_lrc
    kin         = 0.5*SUM(v**2)
    energy      = ( pot + kin ) / REAL ( n )
    energy_sh   = ( pot_sh + kin ) / REAL ( n )
    conserved   = pot_sh + kin + SUM(0.5*p_eta**2/q) + SUM(0.5*p_eta_baro**2/q_baro) + 0.5*p_eps**2/w_eps
    conserved   = conserved + temperature * ( g*eta(1) + SUM(eta(2:m)) + SUM(eta_baro) )
    conserved   = conserved / REAL(n)
    temp_kinet  = 2.0 * kin / g
    temp_config = SUM(f**2) / lap
    density     = REAL(n) / box**3
    pres_virial = density * temperature + vir / box**3

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)' ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Conserved quantity', conserved
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Total energy',       energy
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Shifted energy',     energy_sh
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temp-kinet',         temp_kinet
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Temp-config',        temp_config
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Density',            density
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Virial pressure',    pres_virial
    END IF

  END SUBROUTINE calculate

END PROGRAM md_npt_lj

