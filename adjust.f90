! adjust.f90
! Utility program to allow user to change the density or kinetic energy of MC or MD configuration
PROGRAM adjust

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

  ! Takes in a configuration of atoms or molecules
  ! positions, possibly orientations, and optionally, velocities and angular velocities
  ! Cubic periodic boundary conditions
  ! Adjusts the density by an amount delta_rho
  ! and kinetic energy per particle by an amount delta_kin

  ! Input configuration and output configuration are assumed to be in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1
  ! There is nothing here specific to Lennard-Jones
  ! We assume unit mass and adjust only the translational kinetic energy

  USE, INTRINSIC :: iso_fortran_env,  ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module, ONLY : read_cnf_atoms, read_cnf_mols, write_cnf_atoms, write_cnf_mols
  USE               maths_module,     ONLY : lowercase

  IMPLICIT NONE

  ! Most important variables
  REAL :: box       ! Box length
  REAL :: rho       ! Current density
  REAL :: kin       ! Current kinetic energy per particle
  REAL :: delta_rho ! Desired density change
  REAL :: delta_kin ! Desired kinetic energy change
  REAL :: scale     ! Scaling factor

  LOGICAL :: velocities
  INTEGER :: n, ioerr

  REAL, DIMENSION(:,:), ALLOCATABLE :: r ! Positions (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: v ! Velocities (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: e ! Orientations (3,n) or (0:3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: w ! Angular velocities (3,n)

  CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp'

  CHARACTER(len=10)  :: molecules       ! Character string, input, used to specify molecule_option
  INTEGER            :: molecule_option ! User option for atom, linear, nonlinear or chain molecule
  INTEGER, PARAMETER :: atom = 0, linear = 1, nonlinear = 2, chain = 3

  REAL, PARAMETER :: tol = 1.e-9

  NAMELIST /nml/ delta_rho, delta_kin, velocities, molecules

  ! Set default parameters
  velocities = .FALSE. ! By default, assume MC configuration
  molecules  = 'atom'  ! Options are 'atom', 'chain', 'linear', 'nonlinear'
  delta_rho  = 0.0     ! If neither of these two items is changed ...
  delta_kin  = 0.0     ! ... program will just write out information

  ! Read parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in adjust'
  END IF

  ! Use molecules string to deduce molecule_option
  molecules = lowercase(molecules)

  IF ( INDEX ( molecules, 'chain' ) /= 0 ) THEN

     molecule_option = chain
     WRITE ( unit=output_unit, fmt='(a)' ) 'Chain of atoms, no periodic boundaries, no density'

  ELSE IF ( INDEX ( molecules, 'nonlinear' ) /= 0 ) THEN

     molecule_option = nonlinear
     WRITE ( unit=output_unit, fmt='(a)' ) 'Nonlinear molecules, periodic boundaries'

  ELSE IF ( INDEX ( molecules, 'linear') /= 0 ) THEN

     molecule_option = linear
     WRITE ( unit=output_unit, fmt='(a)' ) 'Linear molecules, periodic boundaries'

  ELSE IF ( INDEX ( molecules, 'atom') /= 0 ) THEN

     molecule_option = atom
     WRITE ( unit=output_unit, fmt='(a)' ) 'Atoms, periodic boundaries'

  ELSE

     WRITE ( unit=error_unit, fmt='(a,a)') 'Unrecognized molecules option: ', molecules
     STOP 'Error in adjust'

  END IF

  IF ( (.NOT. velocities)       .AND. ABS(delta_kin) > tol ) STOP 'No kinetic energy change possible'
  IF ( molecule_option == chain .AND. ABS(delta_rho) > tol ) STOP 'No density change possible'

  ! The read_cnf_atoms routine should work in all cases, just to get n and box (or bond)
  CALL read_cnf_atoms ( filename, n, box )
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  IF ( molecule_option == chain ) THEN
     WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Molecular bond length', box ! box plays role of bond
  ELSE
     rho = REAL(n) / box**3
     WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
     WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               rho
  END IF

  SELECT CASE ( molecule_option )

  CASE ( nonlinear)
     ALLOCATE ( r(3,n), v(3,n), w(3,n), e(0:3,n) )

  CASE default
     ALLOCATE ( r(3,n), v(3,n), w(3,n), e(3,n) )

  END SELECT

  SELECT CASE ( molecule_option )

  CASE ( atom, chain )
     IF ( velocities ) THEN
        CALL read_cnf_atoms ( filename, n, box, r, v )
     ELSE
        CALL read_cnf_atoms ( filename, n, box, r )
     END IF

  CASE ( linear, nonlinear )
     IF ( velocities ) THEN
        CALL read_cnf_mols ( filename, n, box, r, e, v, w )
     ELSE
        CALL read_cnf_mols ( filename, n, box, r, e )
     END IF

  END SELECT

  IF ( velocities ) THEN
     kin = 0.5*SUM(v**2) / REAL(n)
     WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Kinetic energy', kin
  END IF

  IF ( ABS(delta_rho) < tol .AND. ABS(delta_kin) < tol ) THEN

     WRITE ( unit=output_unit, fmt='(a)' ) 'No changes requested'

  ELSE

     IF ( ABS(delta_rho) > tol ) THEN

        IF ( rho+delta_rho < 0.0 ) STOP 'New requested density would be negative'

        scale   = ( rho / (rho+delta_rho) )**(1.0/3.0)
        box     = box * scale
        r(:,:)  = r(:,:) * scale
        rho     = REAL(n) / box**3
        WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'New density', rho

     END IF

     IF ( ABS(delta_kin) > tol ) THEN

        IF ( kin+delta_kin < 0.0 ) STOP 'New requested kinetic energy would be negative'

        scale  = SQRT ( (kin+delta_kin) / kin )
        v(:,:) = v(:,:) * scale
        kin    = 0.5*SUM(v**2) / REAL(n) 
        WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'New kinetic energy', kin
     END IF

     SELECT CASE ( molecule_option )

     CASE ( atom, chain )

        IF ( velocities ) THEN
           CALL write_cnf_atoms ( filename, n, box, r, v )
        ELSE
           CALL write_cnf_atoms ( filename, n, box, r )
        END IF

     CASE ( linear, nonlinear )

        IF ( velocities ) THEN
           CALL write_cnf_mols ( filename, n, box, r, e, v, w )
        ELSE
           CALL write_cnf_mols ( filename, n, box, r, e )
        END IF

     END SELECT

  END IF

  DEALLOCATE ( r, v, e, w )

END PROGRAM adjust

