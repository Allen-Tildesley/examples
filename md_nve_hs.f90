! md_nve_hs.f90
! Molecular dynamics, NVE ensemble, hard spheres
PROGRAM md_nve_hs

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
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in units sigma = 1, mass = 1

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     update, overlap, collide, n, r, v, coltime, partner, gt, lt

  IMPLICIT NONE

  ! Most important variables
  REAL :: box        ! Box length (in units where sigma=1)
  REAL :: vir        ! Total collisional virial
  REAL :: kin        ! Kinetic energy
  REAL :: temp_kinet ! Temperature (conserved)
  REAL :: dt         ! Time step

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  INTEGER            :: i, j, k, ncoll, col_sum, blk, stp, nblock, nstep, ioerr
  REAL               :: tij, t_now, vir_sum
  REAL, DIMENSION(3) :: vcm

  NAMELIST /nml/ nblock, nstep, dt

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_nve_hs'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE, hard spheres'
  CALL introduction

  ! Set sensible default run parameters for testing
  nblock = 10
  nstep  = 2000
  dt     = 0.05

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_nve_hs'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Timestep',                  dt

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box ) ! First call just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',  n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Box (in sigma units)', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',              REAL(n) / box**3
  CALL allocate_arrays
  CALL read_cnf_atoms ( cnf_prefix//inp_tag, n, box, r, v ) ! Second call gets r and v
  r(:,:) = r(:,:) / box                                     ! Convert positions to box=1 units
  r(:,:) = r(:,:) - ANINT ( r(:,:) )                        ! Periodic boundaries
  vcm(:) = SUM ( v(:,:), dim=2 ) / REAL(n)                  ! Centre-of mass velocity
  v(:,:) = v(:,:) - SPREAD ( vcm(:), dim = 2, ncopies = n ) ! Set COM velocity to zero
  kin        = 0.5 * SUM ( v**2 )
  temp_kinet = 2.0 * kin / REAL ( 3*(n-1) )
  v          = v / SQRT ( temp_kinet ) ! We fix the temperature to be 1.0
  kin        = 0.5 * SUM ( v**2 )
  temp_kinet = 2.0 * kin / REAL ( 3*(n-1) )
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature', temp_kinet

  ! Initial overlap check
  IF ( overlap ( box ) ) THEN
     WRITE ( unit=error_unit, fmt='(a)' ) 'Particle overlap in initial configuration'
     STOP 'Error in md_nve_hs'
  END IF

  ! Initial search for collision partners >i
  coltime(:) = HUGE(1.0)
  partner(:) = n
  DO i = 1, n
     CALL update ( i, box, gt ) 
  END DO

  ! Initialize arrays for averaging and write column headings
  col_sum = 0
  vir_sum = 0.0
  CALL run_begin ( calc_variables() )

  ncoll = 0

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        vir_sum = 0.0 ! Zero collisional virial accumulator for this step
        col_sum = 0   ! Zero collision counter for this step
        t_now   = 0.0 ! Keep track of time within this step

        DO ! Loop over collisions within this step

           i   = MINLOC ( coltime, dim=1 ) ! Locate minimum collision time
           j   = partner(i)                ! Collision partner
           tij = coltime(i)                ! Time to collision

           IF ( t_now + tij > dt ) THEN
              CALL advance ( dt - t_now ) ! Advance to end of time step
              EXIT                        ! Exit loop over collisions
           END IF

           CALL advance ( tij )            ! Advance to time of next collision
           CALL collide ( i, j, box, vir ) ! Compute collision dynamics
           col_sum = col_sum + 1           ! Increment collision counter
           vir_sum = vir_sum + vir         ! Increment collisional virial accumulator

           ! Update collision lists
           DO k = 1, n
              IF ( k==i .OR. k==j .OR. partner(k) == i .OR.  partner(k) == j ) THEN
                 CALL update ( k, box, gt ) ! Search for partners >k
              END IF
           END DO
           CALL update ( i, box, lt ) ! Search for partners <i
           CALL update ( j, box, lt ) ! Search for partners <j

        END DO ! End loop over collisions within this step

        ncoll = ncoll + col_sum

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                           ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag, n, box, r*box, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  WRITE ( unit=output_unit, fmt='(a,t40,i15)' ) 'Total collisions', ncoll

  IF ( overlap ( box ) ) STOP 'Particle overlap in final configuration'

  CALL write_cnf_atoms ( cnf_prefix//out_tag, n, box, r*box, v )

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE advance ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time interval over which to advance configuration

    ! Guard against going back in time
    IF ( t < 0.0 ) THEN ! should never happen
       WRITE ( unit=error_unit, fmt='(a,f15.6)' ) 'Negative time step', t
       STOP 'Error in md_nve_hs/advance'
    END IF

    t_now      = t_now + t                 ! Advance current time by t
    coltime(:) = coltime(:) - t            ! Reduce times to next collision by t
    r(:,:)     = r(:,:) + t * v(:,:) / box ! Advance all positions by t (box=1 units)
    r(:,:)     = r(:,:) - ANINT ( r(:,:) ) ! Apply periodic boundaries

  END SUBROUTINE advance

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE averages_module, ONLY : variable_type
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(2) :: variables ! The 2 variables listed below

    ! This routine calculates all variables of interest
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: coll_rate, p_coll
    REAL                :: vol, rho

    ! Preliminary calculations
    vol = box**3        ! Volume
    rho = REAL(n) / vol ! Density

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Collision rate per particle
    ! We average over the time step
    coll_rate = variable_type ( nam = 'Collision rate', val = 2.0*REAL(col_sum)/dt/REAL(n), instant = .FALSE. )

    ! Collisional pressure
    ! ideal + collisional virial / volume averaged over the time step
    p_coll = variable_type ( nam = 'P', val = rho*temp_kinet + vir_sum/dt/vol, instant = .FALSE. )

    ! Collect together for averaging
    variables = [ coll_rate, p_coll ]

  END FUNCTION calc_variables

END PROGRAM md_nve_hs
