! mc_gibbs_lj.f90
! Monte Carlo, Gibbs ensemble
PROGRAM mc_gibbs_lj

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

  ! Takes in a pair of configurations of atoms (positions)
  ! Cubic periodic boundary conditions
  ! Conducts Gibbs ensemble Monte Carlo at the given temperature, total volume and total N
  ! To avoid some inconvenient tests, we disallow configurations in which either box is empty
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Positions r are divided by box length after reading in
  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1

  ! Note that long-range corrections are not included in the acceptance/rejection
  ! of creation and destruction moves

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in mc_module

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, write_cnf_atoms
  USE               averages_module,  ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               maths_module,     ONLY : metropolis, random_integer, random_translate_vector
  USE               mc_module,        ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     potential_1, potential, move, swap, n, r, potential_type

  IMPLICIT NONE

  ! Most important variables
  REAL, DIMENSION(2) :: box         ! Box lengths
  REAL               :: dr_max      ! Maximum MC displacement
  REAL               :: dv_max      ! Maximum MC volume change
  REAL               :: temperature ! Specified temperature
  REAL               :: r_cut       ! Potential cutoff distance

  ! Composite interaction = pot & vir & ovr variables
  TYPE(potential_type), DIMENSION(2) :: total, total_new
  TYPE(potential_type)               :: partial_old, partial_new

  ! Histograms of density, energy, and pressure
  INTEGER, PARAMETER  :: nh = 300
  REAL,    PARAMETER  :: rho_min = 0.0, rho_max = 0.9
  REAL,    PARAMETER  :: rho_del = ( rho_max - rho_min ) / REAL(nh)
  REAL,    PARAMETER  :: eng_min = -3.3, eng_max = 1.2
  REAL,    PARAMETER  :: eng_del = ( eng_max - eng_min ) / REAL(nh)
  REAL, DIMENSION(nh) :: rho_hist, eng_hist

  INTEGER            :: blk, stp, i, nstep, nblock, nswap
  INTEGER            :: iswap, m_acc, x12_try, x12_acc, x21_try, x21_acc, ioerr
  REAL               :: delta, dv, zeta, m1_ratio, m2_ratio, x12_ratio, x21_ratio, v_ratio
  REAL, DIMENSION(3) :: ri
  REAL, DIMENSION(2) :: box_new, vol_new, vol_old

  CHARACTER(len=5), DIMENSION(2), PARAMETER :: cnf_prefix = ['cnf1.','cnf2.']
  CHARACTER(len=3),               PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3),               PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)                          :: sav_tag    = 'sav' ! May be overwritten with block number

  NAMELIST /nml/ nblock, nstep, nswap, temperature, r_cut, dr_max, dv_max

  WRITE( unit=output_unit, fmt='(a)' ) 'mc_gibbs_lj'
  WRITE( unit=output_unit, fmt='(a)' ) 'Monte Carlo, Gibbs ensemble'
  CALL introduction

  CALL RANDOM_SEED () ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nblock      = 10
  nstep       = 10000
  nswap       = 20
  temperature = 1.0
  r_cut       = 2.5
  dr_max      = 0.15
  dv_max      = 10.0

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mc_gibbs_lj'
  END IF

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of blocks',              nblock
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of steps per block',     nstep
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Swap attempts per step',        nswap
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Temperature',                   temperature
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Potential cutoff distance',     r_cut
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum displacement',          dr_max
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Maximum volume change',         dv_max

  ! Read in initial configurations and allocate necessary arrays
  CALL read_cnf_atoms ( cnf_prefix(1)//inp_tag, n(1), box(1) ) ! First call is just to get n and box
  CALL read_cnf_atoms ( cnf_prefix(2)//inp_tag, n(2), box(2) ) ! First call is just to get n and box
  WRITE ( unit=output_unit, fmt='(a,t40,2i15)'   ) 'Number of particles',   n(:)
  WRITE ( unit=output_unit, fmt='(a,t40,2f15.6)' ) 'Simulation box length', box(:)
  WRITE ( unit=output_unit, fmt='(a,t40,2f15.6)' ) 'Density',               REAL(n(:)) / box(:)**3
  CALL allocate_arrays ( box(:), r_cut ) ! Allocate r
  CALL read_cnf_atoms ( cnf_prefix(1)//inp_tag, n(1), box(1), r(:,1:n(1)) )           ! Second call is to get r
  CALL read_cnf_atoms ( cnf_prefix(2)//inp_tag, n(2), box(2), r(:,n(1)+1:n(1)+n(2)) ) ! Second call is to get r
  r(:,1:n(1))           = r(:,1:n(1))           / box(1)  ! Convert positions to box units
  r(:,n(1)+1:n(1)+n(2)) = r(:,n(1)+1:n(1)+n(2)) / box(2)  ! Convert positions to box units
  r(:,:)                = r(:,:) - ANINT ( r(:,:) )       ! Periodic boundaries

  ! Initial energy and overlap check
  total(1) = potential ( 1,      n(1),      box(1), r_cut ) 
  total(2) = potential ( n(1)+1, n(1)+n(2), box(2), r_cut ) 
  IF ( total(1)%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration 1'
     STOP 'Error in mc_gibbs_lj'
  END IF
  IF ( total(2)%ovr ) THEN
     WRITE ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration 2'
     STOP 'Error in mc_gibbs_lj'
  END IF

  ! Initialize arrays for averaging and write column headings
  m1_ratio  = 0.0
  m2_ratio  = 0.0
  x12_ratio = 0.0
  x21_ratio = 0.0
  v_ratio   = 0.0
  CALL run_begin ( calc_variables() )

  ! Zero histograms
  rho_hist(:) = 0.0
  eng_hist(:) = 0.0

  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps

        m_acc = 0

        DO i = 1, n(1) ! Loop over atoms in system 1

           partial_old = potential_1 ( 1, n(1), r(:,i), i, box(1), r_cut ) ! Old atom potential, virial etc

           IF ( partial_old%ovr ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_gibbs_lj'
           END IF

           ri(:) = random_translate_vector ( dr_max/box(1), r(:,i) ) ! Trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )                           ! Periodic boundary correction

           partial_new = potential_1 ( 1, n(1), ri, i, box(1), r_cut ) ! New atom potential, virial etc

           IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

              delta = partial_new%pot - partial_old%pot ! Use cut (but not shifted) potential
              delta = delta / temperature               ! Divide by temperature

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 CALL move ( i, ri )                             ! Update position
                 total(1) = total(1) + partial_new - partial_old ! Update total values
                 m_acc    = m_acc + 1                            ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End test for overlapping configuration

        END DO ! End loop over atoms in system 1

        m1_ratio = REAL(m_acc) / REAL(n(1))

        m_acc = 0

        DO i = n(1)+1, n(1)+n(2) ! Loop over atoms in system 2

           partial_old = potential_1 ( n(1)+1, n(1)+n(2), r(:,i), i, box(2), r_cut ) ! Old atom potential, virial etc

           IF ( partial_old%ovr ) THEN ! should never happen
              WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
              STOP 'Error in mc_gibbs_lj'
           END IF

           ri(:) = random_translate_vector ( dr_max/box(2), r(:,i) ) ! Trial move to new position (in box=1 units)
           ri(:) = ri(:) - ANINT ( ri(:) )                           ! Periodic boundary correction

           partial_new = potential_1 ( n(1)+1, n(1)+n(2), ri, i, box(2), r_cut ) ! New atom potential, virial etc

           IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

              delta = partial_new%pot - partial_old%pot ! Use cut (but not shifted) potential
              delta = delta / temperature               ! Divide by temperature

              IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                 CALL move ( i, ri )                             ! Update position
                 total(2) = total(2) + partial_new - partial_old ! Update total values
                 m_acc    = m_acc + 1                            ! Increment move counter
              END IF ! End accept Metropolis test

           END IF ! End test for overlapping configuration

        END DO ! End loop over atoms in system 2

        m2_ratio = REAL(m_acc) / REAL(n(2))

        x12_try = 0
        x12_acc = 0
        x21_try = 0
        x21_acc = 0

        DO iswap = 1, nswap ! Loop over swap attempts

           CALL RANDOM_NUMBER ( ri )   ! Three uniform random numbers in range (0,1)
           ri = ri - 0.5               ! Now in range (-0.5,+0.5) for box=1 units
           CALL RANDOM_NUMBER ( zeta ) ! Uniform random number on (0,1)

           IF ( zeta > 0.5 ) THEN ! Try swapping 1->2

              x12_try     = x12_try + 1

              IF ( n(1) > 1 ) THEN ! Disallow n(1)->0
                 i           = random_integer ( 1, n(1) )                        ! Choose atom at random in system 1
                 partial_old = potential_1 ( 1, n(1), r(:,i), i, box(1), r_cut ) ! Old atom potential, virial, etc
                 IF ( partial_old%ovr ) THEN ! should never happen
                    WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                    STOP 'Error in mc_gibbs_lj'
                 END IF
                 partial_new = potential_1 ( n(1)+1, n(1)+n(2), ri, 0, box(2), r_cut ) ! New atom potential, virial, etc

                 IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

                    delta = ( partial_new%pot - partial_old%pot ) / temperature ! Use cut (not shifted) potential
                    delta = delta - LOG ( box(2)**3 / REAL ( n(2)+1 ) ) ! Creation in 2
                    delta = delta + LOG ( box(1)**3 / REAL ( n(1) ) )   ! Destruction in 1

                    IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                       CALL swap ( i, ri )               ! Carry out swap
                       total(1) = total(1) - partial_old ! Update total values
                       total(2) = total(2) + partial_new ! Update total values
                       x12_acc  = x12_acc + 1            ! Increment 1->2 move counter
                    END IF ! End accept Metropolis test

                 END IF ! End test for non-overlapping configuration
              END IF ! End test to disallow n(1)->0

           ELSE ! Try swapping 2->1

              x21_try     = x21_try + 1

              IF ( n(2) > 1 ) THEN ! Disallow n(2)->0
                 i           = random_integer ( n(1)+1, n(1)+n(2) )                        ! Choose atom at random in system 2
                 partial_old = potential_1 ( n(1)+1, n(1)+n(2), r(:,i), i, box(2), r_cut ) ! Old atom potential, virial, etc
                 IF ( partial_old%ovr ) THEN ! should never happen
                    WRITE ( unit=error_unit, fmt='(a)') 'Overlap in current configuration'
                    STOP 'Error in mc_gibbs_lj'
                 END IF
                 partial_new = potential_1 ( 1, n(1), ri, 0, box(1), r_cut ) ! New atom potential, virial, etc

                 IF ( .NOT. partial_new%ovr ) THEN ! Test for non-overlapping configuration

                    delta = ( partial_new%pot - partial_old%pot ) / temperature ! Use cut (not shifted) potential
                    delta = delta - LOG ( box(1)**3 / REAL ( n(1)+1 ) ) ! Creation in 1
                    delta = delta + LOG ( box(2)**3 / REAL ( n(2) ) )   ! Destruction in 2

                    IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
                       CALL swap ( i, ri )               ! Carry out swap
                       total(2) = total(2) - partial_old ! Update total values
                       total(1) = total(1) + partial_new ! Update total values
                       x21_acc  = x21_acc + 1            ! Increment 2->1 move counter
                    END IF ! End accept Metropolis test

                 END IF ! End test for non-overlapping configuration

              END IF ! End test to disallow n(2)->0

           END IF ! End choice between trying to swap 1->2 and 2->1

        END DO ! End loop over swap attempts

        x12_ratio = 0.0
        IF ( x12_try > 0 ) x12_ratio = REAL(x12_acc) / REAL(x12_try)
        x21_ratio = 0.0
        IF ( x21_try > 0 ) x21_ratio = REAL(x21_acc) / REAL(x21_try)

        ! Volume move

        v_ratio = 0.0
        CALL RANDOM_NUMBER ( zeta )                ! Uniform on (0,1)
        dv           = dv_max * ( 2.0*zeta - 1.0 ) ! Uniform on (-dv_max,+dv_max)
        vol_old(:)   = box(:)**3                   ! Old volumes
        vol_new(:)   = vol_old(:) + [-dv,dv]       ! New volumes
        box_new(:)   = vol_new(:)**(1.0/3.0)       ! New box lengths
        IF ( MINVAL(box_new) < 2.0*r_cut ) THEN
           WRITE ( unit=error_unit, fmt='(a,2f15.6)') 'Box length too small', box_new
           STOP 'Error in mc_gibbs_lj'
        END IF
        total_new(1) = potential ( 1,      n(1),      box_new(1), r_cut ) 
        total_new(2) = potential ( n(1)+1, n(1)+n(2), box_new(2), r_cut ) 

        IF ( .NOT. ANY ( total_new%ovr ) ) THEN ! Test for non-overlapping configurations

           delta = SUM ( total_new%pot - total%pot )             ! Use cut (but not shifted) potential
           delta = delta / temperature                           ! Divide by temperature
           delta = delta - REAL(n(1))*LOG(vol_new(1)/vol_old(1)) ! Volume scaling in system 1
           delta = delta - REAL(n(2))*LOG(vol_new(2)/vol_old(2)) ! Volume scaling in system 2

           IF ( metropolis ( delta ) ) THEN ! Accept Metropolis test
              total(:) = total_new(:) ! Update total values
              box(:)   = box_new(:)   ! Update box lengths
              v_ratio  = 1.0          ! Set move counter
           END IF ! Reject Metropolis test

        END IF ! End test for non-overlapping configurations

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )
        CALL add_hist

     END DO ! End loop over steps

     CALL blk_end ( blk )                                                                        ! Output block averages
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk                                            ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix(1)//sav_tag, n(1), box(1), box(1)*r(:,1:n(1)) )           ! Save configuration
     CALL write_cnf_atoms ( cnf_prefix(2)//sav_tag, n(2), box(2), box(2)*r(:,n(1)+1:n(1)+n(2)) ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  CALL write_cnf_atoms ( cnf_prefix(1)//out_tag, n(1), box(1), box(1)*r(:,1:n(1))           ) ! Write out final configuration
  CALL write_cnf_atoms ( cnf_prefix(2)//out_tag, n(2), box(2), box(2)*r(:,n(1)+1:n(1)+n(2)) ) ! Write out final configuration

  CALL write_hist

  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  FUNCTION calc_variables () RESULT ( variables )
    USE lrc_module,      ONLY : potential_lrc, pressure_lrc, pressure_delta
    USE averages_module, ONLY : variable_type
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(13) :: variables ! The 13 variables listed below

    ! This function returns all variables of interest in an array, for use in the main program

    ! In this example we simulate using the cut (but not shifted) potential
    ! The values of < p_c >,  < e_c > and < density > should be consistent (for this potential)
    ! For simplicity, long-range corrections are not applied here to give estimates of 
    ! < e_f > and < p_f > for the full (uncut) potential, but this is straightforward to do.
    ! The value of the cut-and-shifted potential is not used, in this example

    TYPE(variable_type) :: m1_r, m2_r, x12_r, x21_r, v_r
    TYPE(variable_type) :: n_1, n_2, density_1, density_2, e1_c, e2_c, p1_c, p2_c
    REAL, DIMENSION(2)  :: vol, rho

    ! Preliminary calculations (m_ratio, total etc are known already)
    vol(:) = box(:)**3
    rho(:) = REAL(n(:)) / vol(:)

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Move, swap, volume exchange acceptance ratios
    m1_r  = variable_type ( nam = 'Move ratio (1)',    val = m1_ratio,  instant = .FALSE. )
    m2_r  = variable_type ( nam = 'Move ratio (2)',    val = m2_ratio,  instant = .FALSE. )
    x12_r = variable_type ( nam = 'Swap ratio (1->2)', val = x12_ratio, instant = .FALSE. )
    x21_r = variable_type ( nam = 'Swap ratio (2->1)', val = x21_ratio, instant = .FALSE. )
    v_r   = variable_type ( nam = 'Volume ratio',      val = v_ratio,   instant = .FALSE. )

    ! Number of particles
    n_1 = variable_type ( nam = 'Number (1)', val = REAL(n(1)) )
    n_2 = variable_type ( nam = 'Number (2)', val = REAL(n(2)) )

    ! Density
    density_1 = variable_type ( nam = 'Density (1)', val = rho(1) )
    density_2 = variable_type ( nam = 'Density (2)', val = rho(2) )

    ! Internal energy per atom for simulated, cut, potential
    ! Ideal gas contribution plus cut (but not shifted) PE divided by N
    e1_c = variable_type ( nam = 'E/N cut (1)', val = 1.5*temperature + total(1)%pot/REAL(n(1)) )
    e2_c = variable_type ( nam = 'E/N cut (2)', val = 1.5*temperature + total(2)%pot/REAL(n(2)) )

    ! Pressure for simulated, cut, potential
    ! delta correction plus ideal gas contribution plus total virial divided by V
    p1_c = variable_type ( nam = 'P cut (1)', val = pressure_delta(rho(1),r_cut) + rho(1)*temperature + total(1)%vir/vol(1) )
    p2_c = variable_type ( nam = 'P cut (2)', val = pressure_delta(rho(2),r_cut) + rho(2)*temperature + total(2)%vir/vol(2) )

    ! Collect together for averaging
    variables = [ m1_r, m2_r, x12_r, x21_r, v_r, n_1, n_2, density_1, density_2, e1_c, e2_c, p1_c, p2_c ]

  END FUNCTION calc_variables

  SUBROUTINE add_hist
    REAL    :: rho, eng
    INTEGER :: k

    rho = REAL(n(1)) / box(1)**3
    k = 1 + FLOOR ( ( rho - rho_min ) / rho_del )
    IF ( k >= 1 .AND. k <= nh ) rho_hist(k) = rho_hist(k) + 1.0

    rho = REAL(n(2)) / box(2)**3
    k = 1 + FLOOR ( ( rho - rho_min ) / rho_del )
    IF ( k >= 1 .AND. k <= nh ) rho_hist(k) = rho_hist(k) + 1.0

    eng = 1.5*temperature + total(1)%pot/REAL(n(1))
    k = 1 + FLOOR ( ( eng - eng_min ) / eng_del )
    IF ( k >= 1 .AND. k <= nh ) eng_hist(k) = eng_hist(k) + 1.0

    eng = 1.5*temperature + total(2)%pot/REAL(n(2))
    k = 1 + FLOOR ( ( eng - eng_min ) / eng_del )
    IF ( k >= 1 .AND. k <= nh ) eng_hist(k) = eng_hist(k) + 1.0

  END SUBROUTINE add_hist

  SUBROUTINE write_hist
    INTEGER                     :: hist_unit, ioerr, k
    REAL                        :: norm, rho, eng
    CHARACTER(len=7), PARAMETER :: filename = 'his.out'

    ! Normalization factor for two data points at each step
    norm     = REAL(2*nstep*nblock)
    rho_hist = rho_hist / ( norm * rho_del )
    eng_hist = eng_hist / ( norm * eng_del )

    OPEN ( newunit = hist_unit, file = filename, status='replace', iostat=ioerr )
    IF ( ioerr /= 0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
       STOP 'Error in write_hist'
    END IF
    DO k = 1, nh
       rho = rho_min + (REAL(k)-0.5) * rho_del
       eng = eng_min + (REAL(k)-0.5) * eng_del
       WRITE ( unit=hist_unit, fmt='(4f15.6)' ) rho, rho_hist(k), eng, eng_hist(k)
    END DO
    CLOSE ( unit=hist_unit )

  END SUBROUTINE write_hist

END PROGRAM mc_gibbs_lj
