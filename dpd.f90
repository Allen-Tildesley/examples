! dpd.f90
! Dissipative particle dynamics
PROGRAM dpd

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
  ! Conducts dissipative particle dynamics using Shardlow or Lowe-Andersen algorithm
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in and we assume mass=1 throughout
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! The range parameter (cutoff distance) is taken as unity

  ! The model is defined in dpd_module
  ! The typical DPD model described by Groot and Warren, J Chem Phys 107, 4423 (1997)
  ! has temperature kT=1, density rho=3, noise level sigma=3, gamma=sigma**2/(2*kT)=4.5
  ! and force strength parameter a=25 (more generally 75*kT/rho).
  ! We recommend a somewhat smaller timestep than their 0.04.
  ! They also give an approximate expression for the pressure, written out at the end for comparison

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               maths_module,     ONLY : lowercase
  USE               dpd_module,       ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     force, lowe, shardlow, p_approx, r, v, f, n, potential_type

  IMPLICIT NONE

  ! Most important variables
  REAL :: box         ! Box length
  REAL :: rho         ! Density
  REAL :: a           ! Force strength parameter
  REAL :: dt          ! Time step
  REAL :: gamma       ! Thermalization rate (inverse time)
  REAL :: temperature ! Temperature (specified)

  ! Composite interaction = pot & vir & lap variables
  TYPE(potential_type) :: total

  INTEGER :: blk, stp, nstep, nblock, ioerr

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number
  CHARACTER(len=10)           :: method

  ! Define a procedure pointer with an interface like that of lowe
  PROCEDURE(lowe), POINTER :: thermalize => NULL()

  NAMELIST /nml/ nblock, nstep, dt, temperature, a, gamma, method

  WRITE ( unit=output_unit, fmt='(a)' ) 'dpd'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Dissipative particle dynamics, constant-NVT ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 and cutoff=1 throughout'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 10000
  dt          = 0.02
  temperature = 1.0
  a           = 75.0 ! actually a*rho/kT: to be multiplied by kT/rho later
  gamma       = 4.5
  method      = 'Lowe'

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in dpd'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',              nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',     nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step',                     dt
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Specified temperature',         temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Force strength a*rho/kT',       a
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Friction / thermal rate gamma', gamma

  IF ( INDEX( lowercase(method), 'shardlow' ) /= 0 ) THEN
     thermalize => shardlow
     WRITE ( unit=output_unit, fmt='(a)' ) 'Shardlow integration method'
     WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'DPD sigma parameter', SQRT ( 2.0 * gamma * temperature )
  ELSE IF ( INDEX( lowercase(method), 'lowe' ) /= 0 ) THEN
     thermalize => lowe
     WRITE ( unit=output_unit, fmt='(a)' ) 'Lowe thermalization method'
     IF ( gamma*dt > 1.0 ) THEN
        WRITE ( unit=error_unit, fmt='(a,f15.6)') 'gamma*dt too large', gamma*dt
        STOP 'Error in dpd'
     END IF
  ELSE
     WRITE ( unit=error_unit, fmt='(a,a)' ) 'Unrecognized thermalization method ', method
     STOP 'Error in dpd'
  END IF

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  rho = REAL(n) / box**3
  a   = a * temperature / rho ! Scale force strength accordingly
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',          rho
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Force strength a', a
  CALL allocate_arrays ( box )
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Initial energy etc
  CALL force ( box, a, total )

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        ! Shardlow or Lowe-Andersen step
        CALL thermalize ( box, temperature, gamma*dt )

        ! Velocity Verlet step
        CALL kick_propagator ( dt/2.0 ) ! Kick half-step
        CALL drift_propagator ( dt )    ! Drift step 
        CALL force ( box, a, total )    ! Force evaluation
        CALL kick_propagator ( dt/2.0 ) ! Kick half-step

        ! Calculate and accumulate quantities for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                           ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! Save configuration

  END DO ! End loop over blocks      

  CALL run_end ( calc_variables() ) ! Output run averages

  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Approx P = ', p_approx ( a, rho, temperature )

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE kick_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    v(:,:) = v(:,:) + t * f(:,:)

  END SUBROUTINE kick_propagator

  SUBROUTINE drift_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt)

    r(:,:) = r(:,:) + t * v(:,:) / box ! Positions in box=1 units
    r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  END SUBROUTINE drift_propagator

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE averages_module, ONLY : variable_type
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(4) :: variables ! The 4 variables listed below

    ! This function returns all variables of interest in an array, for use in the main program

    ! The DPD potential is short ranged, zero at, and beyond, r_cut
    ! so issues of shifted potentials and long-range corrections do not arise

    TYPE(variable_type) :: p_f, e_f, t_k, t_c
    REAL                :: kin, fsq, vol

    ! Preliminary calculations
    kin = 0.5*SUM(v**2) ! Total kinetic energy
    fsq = SUM(f**2)     ! Total squared force
    vol = box**3        ! Volume

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Kinetic temperature
    ! Momentum is conserved, hence 3N-3 degrees of freedom
    t_k = variable_type ( nam = 'T kinetic', val = 2.0*kin/REAL(3*n-3) )

    ! Internal energy per atom
    ! Total KE plus total PE divided N
    e_f = variable_type ( nam = 'E/N', val = (kin+total%pot)/REAL(n) )

    ! Pressure
    ! Ideal gas part plus total virial divided by V
    p_f = variable_type ( nam = 'P', val = rho*temperature + total%vir/vol )

    ! Configurational temperature
    ! Total squared force divided by total Laplacian
    t_c = variable_type ( nam = 'T config', val = fsq/total%lap )

    ! Collect together for averaging
    variables = [ e_f, t_k, t_c, p_f ]

  END FUNCTION calc_variables

END PROGRAM

