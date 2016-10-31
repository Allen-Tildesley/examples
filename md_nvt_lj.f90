! md_nvt_lj.f90
! Molecular dynamics, NVT ensemble
PROGRAM md_nvt_lj
  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE maths_module,     ONLY : random_normals
  USE averages_module,  ONLY : time_stamp, run_begin, run_end, blk_begin, blk_end, blk_add
  USE md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                       force, r, v, f, n, potential_type

  IMPLICIT NONE

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using a measure-preserving algorithm for the Nose-Hoover equations
  ! Nose-Hoover chains are used, following Martyna et al, Molec Phys, 87, 1117 (1996)
  ! and Tuckerman et al J Phys A, 39, 5629 (2006)
  ! To keep this example reasonably simple, we do not subdivide the timesteps with a
  ! Suzuki-Yoshida decomposition, as described in those papers
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in and we assume mass=1 throughout
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: density     ! Density
  REAL :: dt          ! Time step
  REAL :: r_cut       ! Potential cutoff distance
  REAL :: temperature ! Specified temperature
  REAL :: g           ! Number of degrees of freedom
  REAL :: tau         ! Thermostat time scale

  INTEGER, PARAMETER    :: m = 3 ! number of Nose-Hoover chain variables
  REAL,    DIMENSION(m) :: q     ! thermal inertias
  REAL,    DIMENSION(m) :: eta   ! thermal coordinates (needed only to calculate conserved quantity)
  REAL,    DIMENSION(m) :: p_eta ! thermal momenta

  ! Quantities to be averaged
  REAL :: en_s    ! Internal energy (cut-and-shifted ) per atom
  REAL :: p_s     ! Pressure (cut-and-shifted)
  REAL :: en_f    ! Internal energy (full, including LRC) per atom
  REAL :: p_f     ! Pressure (full, including LRC)
  REAL :: tk      ! Kinetic temperature
  REAL :: tc      ! Configurational temperature
  REAL :: cn      ! Conserved energy-like quantity per atom

  ! Composite interaction = pot & cut & vir & lap & ovr variables
  TYPE(potential_type) :: total

  INTEGER :: blk, stp, nstep, nblock, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag = 'inp', out_tag = 'out'
  CHARACTER(len=3)            :: sav_tag = 'sav' ! may be overwritten with block number

  NAMELIST /nml/ nblock, nstep, r_cut, dt, temperature, tau

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_nvt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL introduction ( output_unit )
  CALL time_stamp ( output_unit )

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 1000
  r_cut       = 2.5
  dt          = 0.002
  temperature = 0.7 ! specified temperature
  tau         = 2.0 ! desired thermostat timescale

  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Nose-Hoover chains with m = ', m

  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_nvt_lj'
  END IF

  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Specified temperature',     temperature

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

  ! Initial forces, potential, etc plus overlap check
  CALL force ( box, r_cut, total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in md_nvt_lj'
  END IF
  CALL calculate ( 'Initial values' )

  CALL run_begin ( [ CHARACTER(len=15) :: 'E/N (cut&shift)', 'P (cut&shift)', &
       &            'E/N (full)', 'P (full)',  'T (kin)', 'T (con)', 'Cons/N' ] )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL u4_propagator ( dt/4.0, m, 1 )
        CALL u3_propagator ( dt/2.0 )
        CALL u4_propagator ( dt/4.0, 1, m )

        CALL u2_propagator ( dt/2.0 )
        CALL u1_propagator ( dt )

        CALL force ( box, r_cut, total )
        IF ( total%ovr ) THEN
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Error in md_nvt_lj'
        END IF

        CALL u2_propagator ( dt/2.0 )

        CALL u4_propagator ( dt/4.0, m, 1 )
        CALL u3_propagator ( dt/2.0 )
        CALL u4_propagator ( dt/4.0, 1, m )

        CALL calculate ( )

        ! Calculate all variables for this step
        CALL blk_add ( [en_s,p_s,en_f,p_f,tk,tc,cn] )

     END DO ! End loop over steps

     CALL blk_end ( blk, output_unit )
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! save configuration

  END DO ! End loop over blocks

  CALL run_end ( output_unit )

  CALL force ( box, r_cut, total )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in md_nvt_lj'
  END IF
  CALL calculate ( 'Final values' )
  CALL time_stamp ( output_unit )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )

  CALL deallocate_arrays
  CALL conclusion ( output_unit )

CONTAINS

  SUBROUTINE u1_propagator ( t ) ! U1: velocity Verlet drift step propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt)

    r(:,:) = r(:,:) + t * v(:,:) / box  ! Positions in box=1 units
    r(:,:) = r(:,:) - ANINT ( r(:,:) )  ! Periodic boundaries

  END SUBROUTINE u1_propagator

  SUBROUTINE u2_propagator ( t ) ! U2: velocity Verlet kick step propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    v(:,:) = v(:,:) + t * f(:,:)

  END SUBROUTINE u2_propagator

  SUBROUTINE u3_propagator ( t ) ! U3: thermostat propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    v(:,:) = v(:,:) * EXP ( -t * p_eta(1) / q(1) )
    eta(:) = eta(:) + t * p_eta(:) / q(:)

  END SUBROUTINE u3_propagator

  SUBROUTINE u4_propagator ( t, j_start, j_stop ) ! U4: thermostat propagator
    IMPLICIT NONE
    REAL,    INTENT(in) :: t               ! Time over which to propagate (typically dt/4)
    INTEGER, INTENT(in) :: j_start, j_stop ! Order in which to tackle variables

    INTEGER         :: j, j_stride
    REAL            :: gj, x, c
    REAL, PARAMETER :: c1 = -1.0/2.0, c2 = 1.0/6.0, c3 = -1.0/24.0

    IF ( j_start > j_stop ) THEN
       j_stride = -1
    ELSE
       j_stride = 1
    END IF

    DO j = j_start, j_stop, j_stride ! Loop over each momentum in turn

       IF ( j == 1 ) THEN ! The driver Gj for p_eta(1) is different
          gj = SUM(v**2) - g*temperature
       ELSE
          gj = ( p_eta(j-1)**2 / q(j-1) ) - temperature
       END IF

       IF ( j == m ) THEN ! The equation for p_eta(M) is different

          p_eta(j)  = p_eta(j) + t * gj

       ELSE

          x = t * p_eta(j+1)/q(j+1)

          IF ( x < 0.001 ) THEN ! Guard against small values
             c = 1.0 + x * ( c1 + x * ( c2 + x * c3 ) ) ! Taylor series to order 3
          ELSE
             c = (1.0-EXP(-x))/x
          END IF ! End guard against small values

          p_eta(j) = p_eta(j)*EXP(-x) + t * gj * c

       END IF

    END DO ! End loop over each momentum in turn

  END SUBROUTINE u4_propagator

  SUBROUTINE calculate ( string )
    USE md_module, ONLY : potential_lrc, pressure_lrc
    IMPLICIT NONE
    CHARACTER (len=*), INTENT(in), OPTIONAL :: string

    ! This routine calculates variables of interest and (optionally) writes them out

    REAL :: kin

    kin = 0.5*SUM(v**2)
    tk  = 2.0 * kin / g

    en_s = ( total%pot + kin ) / REAL ( n )        ! total%pot is the total cut-and-shifted PE
    en_f = ( total%cut + kin ) / REAL ( n )        ! total%cut is the total cut (but not shifted) PE
    en_f = en_f + potential_lrc ( density, r_cut ) ! Add LRC
    
    p_s  = density * temperature + total%vir / box**3 ! total%vir is the total virial
    p_f  = p_s + pressure_lrc ( density, r_cut )      ! Add LRC

    tc = SUM(f**2) / total%lap ! total%lap is the total Laplacian

    cn = total%pot + kin + SUM(0.5*p_eta**2/q)
    cn = cn + temperature * ( g*eta(1) + SUM(eta(2:m)) )
    cn = cn / REAL(n)

    IF ( PRESENT ( string ) ) THEN
       WRITE ( unit=output_unit, fmt='(a)'           ) string
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (cut&shift)', en_s
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (cut&shift)',   p_s
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'E/N (full)',      en_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'P (full)',        p_f
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (kin)',         tk
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'T (con)',         tc
       WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Cons/N',          cn
    END IF

  END SUBROUTINE calculate

END PROGRAM md_nvt_lj

