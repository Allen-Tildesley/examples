! md_nvt_poly_lj.f90
! Molecular dynamics, NVT ensemble (NVE option), polyatomic molecules
PROGRAM md_nvt_poly_lj

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

  ! Takes in a configuration of atoms (positions, quaternions, velocities and angular momenta)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics with optional velocity thermalization
  ! Uses no special neighbour lists
  ! The rotational algorithm is described in the text, section 3.3. See A Dullweber, B Leimkuhler, 
  ! R McLachlan, J Chem Phys 107, 5840 (1997) and TF Miller, M Eleftheriou, P Pattnaik, 
  ! A Ndirango, D Newns, GJ Martyna, J Chem Phys 116, 8649 (2002).

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in and we assume mass=1 throughout
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_mols, write_cnf_mols
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               maths_module,     ONLY : q_to_a
  USE               md_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     force, r, e, d, v, ell, f, tau, n, na, db, inertia, potential_type

  IMPLICIT NONE

  ! Most important variables
  REAL    :: box         ! Box length
  REAL    :: dt          ! Time step
  REAL    :: temperature ! Specified temperature
  LOGICAL :: nvt         ! Flag to indicate NVT ensemble

  ! Composite interaction = pot & vir & ovr variables
  TYPE(potential_type) :: total

  INTEGER              :: blk, stp, nstep, nblock, t_interval, ioerr, i, a
  REAL                 :: norm
  REAL, DIMENSION(3)   :: vcm
  REAL, DIMENSION(3,3) :: ai

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf.'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, dt, temperature, t_interval

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_nvt_poly_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVT/NVE ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular mass=1 throughout'

  CALL introduction
  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 20000
  dt          = 0.003
  temperature = 1.0
  t_interval  = 0

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in md_poly_lj'
  END IF

  ! Write out run parameters
  nvt = t_interval > 0 .AND. t_interval < nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',          nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block', nstep
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Time step',                 dt
  IF ( nvt ) THEN
     WRITE ( unit=output_unit, fmt='(a)'           ) 'NVT ensemble'
     WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Thermalization interval', t_interval
     WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature',             temperature
  ELSE
     WRITE ( unit=output_unit, fmt='(a)'   ) 'NVE ensemble'
     t_interval = nstep+1
  END IF

  ! Read in initial configuration and allocate necessary arrays
  CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               REAL(n) / box**3
  CALL allocate_arrays ( box )
  IF ( nvt ) THEN
     CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e ) ! Second call gets r, e
     CALL ran_velocities
  ELSE
     CALL read_cnf_mols ( cnf_prefix//inp_tag, n, box, r, e, v, ell ) ! Second call gets r, e, v and ell
     vcm(:) = SUM ( v(:,:), dim=2 ) / REAL(n)                         ! Centre-of mass velocity
     v(:,:) = v(:,:) - SPREAD ( vcm(:), dim = 2, ncopies = n )        ! Set COM velocity to zero
  END IF
  r(:,:) = r(:,:) / box              ! Convert positions to box units
  r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  ! Calculate all bond vectors
  DO i = 1, n ! Loop over molecules
     norm   = SQRT(SUM(e(:,i)**2)) 
     e(:,i) = e(:,i) / norm ! Ensure normalized quaternions
     ai = q_to_a ( e(:,i) ) ! Rotation matrix for i
     DO a = 1, na ! Loop over all atoms
        d(:,a,i) = MATMUL ( db(:,a), ai ) ! NB: equivalent to ai_T*db, ai_T=transpose of ai
     END DO ! End loop over all atoms
  END DO ! End loop over molecules

  ! Initial forces, potential, etc plus overlap check
  CALL force ( box, total )
  IF ( total%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     STOP 'Error in md_poly_lj'
  END IF

  ! Initialize arrays for averaging and write column headings
  CALL run_begin ( calc_variables() )

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        IF ( nvt .AND. MOD(stp,t_interval) == 0 ) THEN
           CALL ran_velocities
        END IF

        CALL kick_propagator ( 0.5 * dt ) ! Half-kick step

        CALL drift_translate ( dt ) ! Drift step for positions

        ! Succession of drift steps for rotation about body-fixed axes
        ! Depending on the values of the moments of inertia, a different nested
        ! sequence of axes may produce better or worse energy conservation
        CALL drift_rotate ( 1, 0.5*dt )
        CALL drift_rotate ( 2, 0.5*dt )
        CALL drift_rotate ( 3,     dt )
        CALL drift_rotate ( 2, 0.5*dt )
        CALL drift_rotate ( 1, 0.5*dt )

        ! Calculate all bond vectors
        DO i = 1, n ! Loop over molecules
           norm   = SQRT(SUM(e(:,i)**2)) 
           e(:,i) = e(:,i) / norm ! Guard against cumulative roundoff
           ai = q_to_a ( e(:,i) ) ! Rotation matrix for i
           DO a = 1, na ! Loop over all atoms
              d(:,a,i) = MATMUL ( db(:,a), ai ) ! NB: equivalent to ai_T*db, ai_T=transpose of ai
           END DO ! End loop over all atoms
        END DO ! End loop over molecules

        CALL force ( box, total ) ! Force evaluation
        IF ( total%ovr ) THEN
           WRITE ( unit=error_unit, fmt='(a)') 'Overlap in configuration'
           STOP 'Error in md_poly_lj'
        END IF

        CALL kick_propagator ( 0.5 * dt ) ! Half-kick step

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                                  ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk                      ! Number configuration by block
     CALL write_cnf_mols ( cnf_prefix//sav_tag, n, box, r*box, e, v, ell ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  CALL force ( box, total )
  IF ( total%ovr ) THEN ! should never happen
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in final configuration'
     STOP 'Error in md_poly_lj'
  END IF

  CALL write_cnf_mols ( cnf_prefix//out_tag, n, box, r*box, e, v, ell ) ! Write out final configuration

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  SUBROUTINE kick_propagator ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt/2)

    ! Advances velocities and body-fixed angular momenta
    v(:,:)   = v(:,:)   + t * f(:,:)   ! Linear momenta are equivalent to linear velocities
    ell(:,:) = ell(:,:) + t * tau(:,:) ! Space-fixed angular momenta ell and torques tau

  END SUBROUTINE kick_propagator

  SUBROUTINE drift_translate ( t )
    IMPLICIT NONE
    REAL, INTENT(in) :: t ! Time over which to propagate (typically dt)

    ! Advances positions
    
    r(:,:) = r(:,:) + t * v(:,:) / box ! Drift step (positions in box=1 units)
    r(:,:) = r(:,:) - ANINT ( r(:,:) ) ! Periodic boundaries

  END SUBROUTINE drift_translate

  SUBROUTINE drift_rotate ( xyz, t )
    USE maths_module, ONLY : q_to_a, rotate_quaternion
    IMPLICIT NONE
    INTEGER, INTENT(in) :: xyz ! Body-fixed axis about which to rotate
    REAL,    INTENT(in) :: t   ! Time over which to propagate (typically dt or dt/2)

    ! Advances quaternion orientations about a specified axis
    
    INTEGER              :: i
    REAL                 :: w_mag
    REAL, DIMENSION(3)   :: w_hat
    REAL, DIMENSION(3,3) :: ai

    DO i = 1, n
       ai     = q_to_a ( e(:,i) )                              ! Rotation matrix
       w_hat  = ai(xyz,:)                                      ! Space-fixed axis about which to rotate
       w_mag  = DOT_PRODUCT ( ell(:,i), w_hat ) / inertia(xyz) ! Angular velocity about this axis
       e(:,i) = rotate_quaternion ( w_mag*t, w_hat, e(:,i) )   ! Rotate by specified angle
    END DO

  END SUBROUTINE drift_rotate

  SUBROUTINE ran_velocities
    USE maths_module, ONLY : random_normals

    INTEGER              :: i
    REAL                 :: factor
    REAL, DIMENSION(3)   :: v_cm, factors
    REAL, DIMENSION(3,3) :: ai

    CALL random_normals ( 0.0, 1.0, v )                     ! Random velocities
    v_cm(:) = SUM ( v(:,:), dim=2 ) / REAL ( n )            ! Compute centre of mass velocity
    v(:,:)  = v(:,:) - SPREAD ( v_cm(:), dim=2, ncopies=n ) ! Set net momentum to zero
    factor  = SUM ( v**2 ) / REAL(3*n-3)                    ! Estimate of kinetic temperature
    factor  = SQRT ( temperature / factor )                 ! Necessary correction factor
    v       = factor * v                                    ! Make correction

    CALL random_normals ( 0.0, 1.0, ell )                  ! Random body-fixed angular momenta
    factors = SUM ( ell**2, dim=2 ) / (REAL(n)*inertia(:)) ! Estimate of kinetic temperatures
    factors = SQRT ( temperature / factors )               ! Necessary correction factors
    ell     = SPREAD(factors,dim=2,ncopies=n) * ell        ! Make corrections

    ! Convert to space-fixed angular momenta
    DO i = 1, n ! Loop over molecules
       ai       = q_to_a ( e(:,i) )        ! Rotation matrix for i
       ell(:,i) = MATMUL ( ell(:,i), ai  ) ! NB: equivalent to ell_s = ai_T*ell_b, ai_T=transpose of ai
    END DO ! End loop over molecules

  END SUBROUTINE ran_velocities

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE averages_module, ONLY : variable_type, msd
    USE maths_module,    ONLY : q_to_a
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(5) :: variables ! The 5 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type)  :: e_sf, p_sf, t_t, t_r, conserved_msd
    REAL                 :: vol, rho, kin_t, kin_r, eng, tmp_t, tmp_r
    INTEGER              :: i
    REAL, DIMENSION(3)   :: ell_i
    REAL, DIMENSION(3,3) :: ai

    ! Preliminary calculations
    vol   = box**3        ! Volume
    rho   = REAL(n) / vol ! Density
    kin_t = 0.5*SUM(v**2) ! Translational kinetic energy

    kin_r = 0.0
    DO i = 1, n ! Loop over molecules
       ai    = q_to_a ( e(:,i) )       ! Rotation matrix for i
       ell_i = MATMUL ( ai, ell(:,i) ) ! Get body-fixed angular momentum
       kin_r = kin_r + SUM((ell_i**2)/inertia)
    END DO ! End loop over molecules
    kin_r = 0.5*kin_r ! Rotational kinetic energy

    tmp_t = 2.0 * kin_t / REAL(3*n-3) ! Remove three degrees of freedom for momentum conservation
    tmp_r = 2.0 * kin_r / REAL(3*n)   ! 3N degrees of rotational freedom
    eng   = kin_t + kin_r + total%pot ! Total energy for simulated system

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Internal energy (shifted-force potential) per atom
    ! Total translational and rotational KE plus total PE divided by N
    IF ( nvt ) THEN
       e_sf = variable_type ( nam = 'E/N shifted force', val = 3.0*temperature+total%pot/REAL(n) )
    ELSE
       e_sf = variable_type ( nam = 'E/N shifted force', val = eng/REAL(n) )
    END IF

    ! Pressure (shifted-force potential)
    ! Ideal gas contribution plus total virial divided by V
    IF ( nvt ) THEN
       p_sf = variable_type ( nam = 'P shifted force', val = rho*temperature + total%vir/vol )
    ELSE
       p_sf = variable_type ( nam = 'P shifted force', val = rho*tmp_t + total%vir/vol )
    END IF

    ! Kinetic translational temperature
    t_t = variable_type ( nam = 'T translational', val = tmp_t )

    ! Kinetic rotational temperature
    t_r = variable_type ( nam = 'T rotational', val = tmp_r )

    ! Mean-squared deviation of conserved energy per atom
    conserved_msd = variable_type ( nam = 'Conserved MSD', val = eng/REAL(n), &
         &                          method = msd, e_format = .TRUE., instant = .FALSE. )

    ! Collect together for averaging
    variables = [ e_sf, p_sf, t_t, t_r, conserved_msd ]

  END FUNCTION calc_variables

END PROGRAM md_nvt_poly_lj

