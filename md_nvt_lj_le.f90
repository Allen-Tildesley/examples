! md_nvt_lj_le.f90
! MD, NVT ensemble, Lees-Edwards boundaries
PROGRAM md_nvt_lj_le

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
  ! Cubic periodic boundary conditions, with Lees-Edwards shear
  ! Conducts molecular dynamics, SLLOD algorithm, with isokinetic thermostat
  ! Refs: Pan et al J Chem Phys 122 094114 (2005)

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in and we assume mass=1 throughout
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     force, r, v, f, n, potential_type

  IMPLICIT NONE

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: density     ! Density
  REAL :: dt          ! Time step
  REAL :: strain_rate ! Strain_rate (velocity gradient) dv_x/dr_y
  REAL :: strain      ! Strain (integrated velocity gradient) dr_x/dr_y

  ! Composite interaction = pot & vir & lap & ovr variables
  TYPE(potential_type) :: total

  INTEGER            :: blk, stp, nstep, nblock, ioerr
  REAL, DIMENSION(3) :: vcm
  REAL, PARAMETER    :: tol = 1.0e-6

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, dt, strain_rate

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_nvt_lj_le'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVT ensemble, Lees-Edwards'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL introduction

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 100000
  dt          = 0.005
  strain_rate = 0.04

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_nvt_lj_le'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step',                 dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Strain rate',               strain_rate

  ! Insist that strain be zero (i.e. an integer) at end of each block
  strain = strain_rate * dt * REAL(nstep)
  strain = strain - ANINT ( strain )
  IF ( ABS ( strain ) > tol ) THEN
     WRITE ( unit=error_unit, fmt='(a,es15.6)') 'Strain must be zero at end of block', strain
     STOP 'Error in md_nvt_lj_le'
  END IF

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density', density
  CALL allocate_arrays ( box )
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v
  strain = 0.0                                              ! For simplicity assume that this is true
  r(:,:) = r(:,:) / box                                     ! Convert positions to box units
  r(1,:) = r(1,:) - ANINT ( r(2,:) ) * strain               ! Extra correction (box=1 units)
  r(:,:) = r(:,:) - ANINT ( r(:,:) )                        ! Periodic boundaries (box=1 units)
  vcm(:) = SUM ( v(:,:), dim=2 ) / REAL(n)                  ! Centre-of mass velocity
  v(:,:) = v(:,:) - SPREAD ( vcm(:), dim = 2, ncopies = n ) ! Set COM velocity to zero

  ! Initial forces, potential, etc plus overlap check
  CALL force ( box, strain, total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in md_nvt_lj_le'
  END IF

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! Isokinetic SLLOD algorithm (Pan et al)

        CALL a_propagator  ( dt/2.0 )
        CALL b1_propagator ( dt/2.0 )

        CALL force ( box, strain, total )
        IF ( total%ovr ) THEN
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Error in md_nvt_lj_le'
        END IF

        CALL b2_propagator ( dt )
        CALL b1_propagator ( dt/2.0 )
        CALL a_propagator  ( dt/2.0 )

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                           ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() )

  CALL force ( box, strain, total )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in md_nvt_lj_le'
  END IF

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE a_propagator ( t ) ! A propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    REAL :: x

    x = t * strain_rate ! Change in strain (dimensionless)

    r(1,:) = r(1,:) + x * r(2,:)        ! Extra strain term
    r(:,:) = r(:,:) + t * v(:,:) / box  ! Drift half-step (positions in box=1 units)
    strain = strain + x                 ! Advance strain and hence boundaries
    strain = strain - ANINT ( strain )  ! Keep strain within (-0.5,0.5)

    r(1,:) = r(1,:) - ANINT ( r(2,:) ) * strain ! Extra PBC correction (box=1 units)
    r(:,:) = r(:,:) - ANINT ( r(:,:) )          ! Periodic boundaries (box=1 units)

  END SUBROUTINE a_propagator

  SUBROUTINE b1_propagator ( t ) ! B1 propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    REAL :: x, c1, c2, g

    x = t * strain_rate ! Change in strain (dimensionless)

    c1 = x * SUM ( v(1,:)*v(2,:) ) / SUM ( v(:,:)**2 )
    c2 = ( x**2 ) * SUM ( v(2,:)**2 ) / SUM ( v(:,:)**2 )
    g  = 1.0 / SQRT ( 1.0 - 2.0*c1 + c2 )

    v(1,:) = v(1,:) - x*v(2,:)
    v(:,:) = g * v(:,:)

  END SUBROUTINE b1_propagator

  SUBROUTINE b2_propagator ( t ) ! B2 propagator
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt)

    REAL :: alpha, beta, h, e, dt_factor, prefactor

    alpha = SUM ( f(:,:)*v(:,:) ) / SUM ( v(:,:)**2 )
    beta  = SQRT ( SUM ( f(:,:)**2 ) / SUM ( v(:,:)**2 ) )
    h     = ( alpha + beta ) / ( alpha - beta )
    e     = EXP ( -beta * t )

    dt_factor = ( 1 + h - e - h / e ) / ( ( 1 - h ) * beta )
    prefactor = ( 1 - h ) / ( e - h / e )

    v(:,:) = prefactor * ( v(:,:) + dt_factor * f(:,:) )

  END SUBROUTINE b2_propagator

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE averages_module, ONLY : variable_type, msd
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(6) :: variables ! The 6 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    ! In this example we simulate using a specified potential (e.g. WCA LJ)
    ! which goes to zero smoothly at the cutoff to highlight energy conservation
    ! so long-range corrections do not arise

    TYPE(variable_type) :: e_s, p_s, t_k, t_c, eta, conserved_msd
    REAL                :: vol, rho, kin, fsq, tmp, kyx, eng
    REAL, PARAMETER     :: tol = 1.0e-6

    ! Preliminary calculations
    vol = box**3                   ! Volume
    rho = REAL(n) / vol            ! Density
    kin = 0.5*SUM(v**2)            ! NB v(:,:) are taken to be peculiar velocities
    fsq = SUM(f**2)                ! Total squared force
    tmp = 2.0 * kin / REAL(3*n-3)  ! Remove three degrees of freedom for momentum conservation
    kyx = SUM(v(1,:)*v(2,:)) / vol ! Kinetic part of off-diagonal pressure tensor
    eng = kin + total%pot          ! Total energy

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Internal energy per atom
    ! Total KE plus total PE divided by N
    e_s = variable_type ( nam = 'E/N', val = eng/REAL(n) )

    ! Pressure
    ! Ideal gas contribution plus total virial divided by V 
    p_s = variable_type ( nam = 'P', val = rho*tmp + total%vir/vol )   

    ! Kinetic temperature
    t_k = variable_type ( nam = 'T kinetic', val = tmp )

    ! Configurational temperature
    ! Total squared force divided by total Laplacian
    t_c = variable_type ( nam = 'T config', val = fsq/total%lap )

    ! Shear viscosity
    IF ( ABS(strain_rate) < tol ) THEN ! Guard against simulation with zero strain rate
       eta = variable_type ( nam = 'Shear viscosity', val = 0.0 )
    ELSE
       eta = variable_type ( nam = 'Shear viscosity', val = -( kyx+total%pyx/vol ) / strain_rate )
    END IF

    ! MSD of conserved kinetic energy
    conserved_msd = variable_type ( nam = 'Conserved MSD', val = kin/REAL(n), &
         &                          method = msd, e_format = .TRUE., instant = .FALSE. )

    ! Collect together for averaging
    variables = [ e_s, p_s, t_k, t_c, eta, conserved_msd ]

  END FUNCTION calc_variables

END PROGRAM md_nvt_lj_le

