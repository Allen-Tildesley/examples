! md_npt_lj.f90
! Molecular dynamics, NPT ensemble
PROGRAM md_npt_lj

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             !
  ! to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     !
  ! published by Oxford University Press ("the publishers").                                       !
  !                                                                                                !
  ! LICENCE                                                                                        !
  ! Creative Commons CC0 Public Domain Dedication.                                                 !
  ! To the extent possible under law, the authors have dedicated all copyright and related         !
  ! and neighboring rights to this software to the PUBLIC domain worldwide.                        !
  ! This software is distributed without any warranty.                                             !
  ! You should have received a copy of the CC0 Public Domain Dedication along with this software.  !
  ! If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               !
  !                                                                                                !
  ! DISCLAIMER                                                                                     !
  ! The authors and publishers make no warranties about the software, and disclaim liability       !
  ! for all uses of the software, to the fullest extent permitted by applicable law.               !
  ! The authors and publishers do not recommend use of this software for any purpose.              !
  ! It is made freely available, solely to clarify points made in the text. When using or citing   !
  ! the software, you should not imply endorsement by the authors or publishers.                   !
  !------------------------------------------------------------------------------------------------!

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

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               maths_module,     ONLY : random_normal, random_normals
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     force, r, v, f, n, potential_type

  IMPLICIT NONE

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: dt          ! Time step
  REAL :: r_cut       ! Potential cutoff distance
  REAL :: temperature ! Specified temperature
  REAL :: pressure    ! Specified pressure
  REAL :: g           ! Number of degrees of freedom
  REAL :: tau         ! Thermostat time scale
  REAL :: tau_baro    ! Barostat time scale
  REAL :: w_eps       ! Barostat inertia
  REAL :: eps         ! Strain
  REAL :: p_eps       ! Strain momentum 
  REAL :: box0        ! Initial box length 

  INTEGER, PARAMETER    :: m = 3      ! Number of Nose-Hoover chain variables
  REAL,    DIMENSION(m) :: q          ! Thermal inertias
  REAL,    DIMENSION(m) :: eta        ! Thermal coordinates (needed only to calculate conserved quantity)
  REAL,    DIMENSION(m) :: p_eta      ! Thermal momenta
  REAL,    DIMENSION(m) :: q_baro     ! Barostat thermal inertias
  REAL,    DIMENSION(m) :: eta_baro   ! Barostat thermal coordinates
  REAL,    DIMENSION(m) :: p_eta_baro ! Barostat thermal momenta

  ! Composite interaction = pot & cut & vir & lap & ovr variables
  TYPE(potential_type) :: total

  INTEGER            :: blk, stp, nstep, nblock, ioerr
  REAL, DIMENSION(3) :: vcm

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, r_cut, dt, temperature, pressure, tau, tau_baro

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_npt_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NPT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL introduction

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 20000
  r_cut       = 2.5
  dt          = 0.005
  temperature = 1.0  ! Specified temperature
  pressure    = 0.99 ! Specified pressure
  tau         = 2.0  ! Desired thermostat timescale
  tau_baro    = 2.0  ! Desired barostat timescale

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_npt_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Potential cutoff distance', r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Specified temperature',     temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Specified pressure',        pressure
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Thermostat timescale',      tau
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Barostat timescale',        tau_baro
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Nose-Hoover chain length',  m

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n) / box**3
  CALL allocate_arrays ( box, r_cut )
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v
  r(:,:) = r(:,:) / box                                     ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) )                        ! Periodic boundaries
  vcm(:) = SUM ( v(:,:), dim=2 ) / REAL(n)                  ! Centre-of mass velocity
  v(:,:) = v(:,:) - SPREAD ( vcm(:), dim = 2, ncopies = n ) ! Set COM velocity to zero

  ! Initial values of thermal variables
  CALL RANDOM_SEED()
  g    = REAL ( 3*(n-1) )
  q    = temperature * tau**2   
  q(1) = g * temperature * tau**2
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.6))' ) 'Thermal inertias Q', q
  eta(:) = 0.0
  CALL random_normals ( 0.0, SQRT(temperature), p_eta(:)  )
  p_eta(:) = p_eta(:) * SQRT(q(:))
  q_baro = temperature * tau_baro**2   
  WRITE ( unit=output_unit, fmt='(a,t40,*(f15.6))' ) "Barostat thermal inertias Q'", q_baro
  eta_baro(:) = 0.0
  CALL random_normals ( 0.0, SQRT(temperature), p_eta_baro(:)  )
  p_eta_baro(:) = p_eta_baro(:) * SQRT(q_baro(:))
  w_eps = g * temperature * tau_baro**2
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) "Barostat inertia W", w_eps
  box0  = box ! Reference box length for strain
  eps   = 0.0 ! Initial strain; generally eps = log ( box/box0 )
  p_eps = random_normal ( 0.0, SQRT(temperature*w_eps) ) ! strain momentum

  ! Initial forces, potential, etc plus overlap check
  CALL force ( box, r_cut, total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in md_npt_lj'
  END IF

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        CALL u4_propagator ( dt/4.0, m, 1 ) ! Inwards order
        CALL u3_propagator ( dt/2.0 )
        CALL u4_propagator ( dt/4.0, 1, m ) ! Outwards order

        CALL u2p_propagator ( dt/2.0 )
        CALL u2_propagator  ( dt/2.0 )

        CALL u1_propagator  ( dt )

        CALL force ( box, r_cut, total )
        IF ( total%ovr ) THEN
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Error in md_npt_lj'
        END IF

        CALL u2_propagator  ( dt/2.0 )
        CALL u2p_propagator ( dt/2.0 )

        CALL u4_propagator ( dt/4.0, m, 1 ) ! Inwards order
        CALL u3_propagator ( dt/2.0 )
        CALL u4_propagator ( dt/4.0, 1, m ) ! Outwards order

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                           ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() )

  CALL force ( box, r_cut, total )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in md_npt_lj'
  END IF

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE u1_propagator ( t ) ! U1 and U1' combined: position and strain drift propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt)

    REAL            :: x, c
    REAL, PARAMETER :: c1 = -1.0/2.0, c2 = 1.0/6.0, c3 = -1.0/24.0

    ! U1 part
    ! The propagator for r(:,:) looks different from the formula given on p142 of the text.
    ! However it is easily derived from that formula, bearing in mind that coordinates in
    ! this program are divided by the box length, which is itself updated in this routine.

    x = t * p_eps / w_eps ! Time step * time derivative of strain

    IF ( x < 0.001 ) THEN ! Guard against small values
       c = 1.0 + x * ( c1 + x * ( c2 + x * c3 ) ) ! Taylor series to order 3
    ELSE
       c = (1.0-EXP(-x))/x
    END IF ! End guard against small values

    r(:,:) = r(:,:) + c * t * v(:,:) / box ! Positions in box=1 units
    r(:,:) = r(:,:) - ANINT ( r(:,:) )     ! Periodic boundaries

    ! U1' part
    ! Because we divide by box above, it is important to update eps afterwards
    ! If we did not use box-scaled coordinates this would not matter

    eps = eps + x         ! Update strain
    box = box0 * EXP(eps) ! Update box length

  END SUBROUTINE u1_propagator

  SUBROUTINE u2_propagator ( t ) ! U2: velocity kick step propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    REAL            :: x, c, alpha
    REAL, PARAMETER :: c1 = -1.0/2.0, c2 = 1.0/6.0, c3 = -1.0/24.0

    alpha = 1.0 + 3.0 / g
    x = t * alpha * p_eps / w_eps

    IF ( x < 0.001 ) THEN ! Guard against small values
       c = 1.0 + x * ( c1 + x * ( c2 + x * c3 ) ) ! Taylor series to order 3
    ELSE
       c = (1.0-EXP(-x))/x
    END IF ! End guard against small values

    v(:,:) = v(:,:)*EXP(-x) + c * t * f(:,:)

  END SUBROUTINE u2_propagator

  SUBROUTINE u2p_propagator ( t ) ! U2': strain momentum propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    REAL :: alpha, pv

    alpha = 1.0 + 3.0 / g
    pv    = alpha * SUM(v**2) / 3.0 + total%vir ! total%vir is the total virial
    p_eps = p_eps + 3.0 * ( pv - pressure*box**3 ) * t

  END SUBROUTINE u2p_propagator

  SUBROUTINE u3_propagator ( t ) ! U3 and U3' combined: thermostat propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    ! U3 part

    v(:,:) = v(:,:) * EXP ( -t * p_eta(1) / q(1) )
    eta(:) = eta(:) + t * p_eta(:) / q(:)

    ! U3' part

    eta_baro(:) = eta_baro(:) + t * p_eta_baro(:) / q_baro(:)
    p_eps       = p_eps * EXP ( -t * p_eta_baro(1) / q_baro(1) )

  END SUBROUTINE u3_propagator

  SUBROUTINE u4_propagator ( t, j_start, j_stop ) ! U4 and U4' combined: thermostat propagator
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

    ! U4 part

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

    ! U4' part

    DO j = j_start, j_stop, j_stride ! Loop over each momentum in turn

       IF ( j == 1 ) THEN ! The driver Gj for p_eta_baro(1) is different
          gj = p_eps**2/w_eps - temperature
       ELSE
          gj = ( p_eta_baro(j-1)**2 / q_baro(j-1) ) - temperature
       END IF

       IF ( j == m ) THEN ! The equation for p_eta_baro(M) is different

          p_eta_baro(j)  = p_eta_baro(j) + t * gj

       ELSE

          x = t * p_eta_baro(j+1)/q_baro(j+1)

          IF ( x < 0.001 ) THEN ! Guard against small values
             c = 1.0 + x * ( c1 + x * ( c2 + x * c3 ) ) ! Taylor series to order 3
          ELSE
             c = (1.0-EXP(-x))/x
          END IF ! End guard against small values

          p_eta_baro(j) = p_eta_baro(j)*EXP(-x) + t * gj * c

       END IF

    END DO ! End loop over each momentum in turn

  END SUBROUTINE u4_propagator

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE lrc_module,      ONLY : potential_lrc, pressure_lrc
    USE averages_module, ONLY : variable_type, msd
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(10) :: variables ! The 10 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: e_s, p_s, e_f, p_f, t_k, t_c, density
    TYPE(variable_type) :: vol_msd, conserved, conserved_msd
    REAL                :: vol, rho, kin, fsq, ext, enp, eng

    ! Preliminary calculations
    vol = box**3        ! Volume
    rho = REAL(n) / vol ! Density
    kin = 0.5*SUM(v**2) ! Kinetic energy
    fsq = SUM(f**2)     ! Total squared force
    ext = SUM(0.5*p_eta**2/q) + SUM(0.5*p_eta_baro**2/q_baro) + 0.5*p_eps**2/w_eps &
         & + temperature * ( g*eta(1) + SUM(eta(2:m)) + SUM(eta_baro) ) ! Extra terms for conserved variable
    eng = kin + total%pot ! Total energy (cut-and-shifted)

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Internal energy (cut-and-shifted ) per atom
    ! Total KE plus total cut-and-shifted PE divided by N
    e_s = variable_type ( nam = 'E/N cut&shifted', val = eng/REAL(n) )

    ! Internal energy (full, including LRC) per atom
    ! LRC plus total KE plus total cut (but not shifted) PE divided by N
    e_f = variable_type ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + (kin+total%cut)/REAL(n) )

    ! Pressure (cut-and-shifted)
    ! Ideal gas contribution plus total virial divided by V
    p_s = variable_type ( nam = 'P cut&shifted', val = rho*temperature + total%vir/vol )

    ! Pressure (full, including LRC)
    ! LRC plus ideal gas contribution plus total virial divided by V
    p_f = variable_type ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*temperature + total%vir/vol )

    ! Kinetic temperature
    t_k = variable_type ( nam = 'T kinetic', val = 2.0*kin/g )

    ! Configurational temperature
    ! Total squared force divided by total Laplacian
    t_c = variable_type ( nam = 'T config', val = fsq/total%lap )

    ! Density
    density = variable_type ( nam = 'Density', val = rho )

    ! Conserved energy-like quantity per atom
    ! Energy plus PV term plus extra terms defined above
    enp = eng + pressure*vol
    conserved = variable_type ( nam = 'Conserved/N', val = (enp+ext)/REAL(n) )

    ! Volume MSD
    vol_msd = variable_type ( nam = 'Volume MSD', val = vol, method = msd, instant = .FALSE. )

    ! Mean-squared deviation of conserved energy-like quantity
    ! Energy plus PV term plus extra terms defined above
    enp = eng + pressure*vol
    conserved_msd = variable_type ( nam = 'Conserved MSD', val = (enp+ext)/REAL(n), &
         &                          method = msd, e_format = .TRUE., instant = .FALSE. )

    ! Collect together for averaging
    variables = [ e_s, p_s, e_f, p_f, t_k, t_c, density, conserved, vol_msd, conserved_msd ]

  END FUNCTION calc_variables

END PROGRAM md_npt_lj

