! adjust_density.f90
! Utility program to allow user to change the density of MC or MD configuration
PROGRAM adjust_density

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE config_io_module, ONLY : read_cnf_atoms, read_cnf_mols, write_cnf_atoms, write_cnf_mols
  USE maths_module,     ONLY : lowercase

  IMPLICIT NONE

  ! Takes in a configuration of atoms or molecules (positions, possibly orientations, and optionally, velocities)
  ! Cubic periodic boundary conditions
  ! Adjusts the density by an amount delta

  ! Input configuration and output configuration are assumed to be in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1
  ! There is nothing here specific to Lennard-Jones

  ! Most important variables
  REAL :: box     ! Box length
  REAL :: density ! Current density
  real :: delta   ! Desired density change
  REAL :: scale   ! Scaling factor

  LOGICAL :: velocities
  INTEGER :: n, ioerr

  REAL, DIMENSION(:,:), ALLOCATABLE :: r ! Positions (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: v ! Velocities (3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: e ! Orientations (3,n) or (0:3,n)
  REAL, DIMENSION(:,:), ALLOCATABLE :: w ! Angular velocities (3,n)

  CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp'

  CHARACTER(len=10)  :: molecules       ! Character string, input, used to specify molecule_option
  INTEGER            :: molecule_option ! User option for atoms, linear, nonlinear or chain molecule
  INTEGER, PARAMETER :: atoms = 0, linear = 1, nonlinear = 2, chain = 3

  REAL, PARAMETER :: tol = 1.e-6

  NAMELIST /nml/ delta, velocities, molecules

  ! Set default parameters
  velocities = .FALSE. ! By default, assume MC configuration
  molecules  = 'atoms' ! Options are 'atoms', 'chain', 'linear', 'nonlinear'
  delta      = 0.0     ! If not changed, program will just write out information

  ! Read parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in adjust_density'
  END IF

  ! Use molecules string to deduce molecule_option
  molecules = lowercase(molecules)

  IF ( INDEX ( molecules, 'chain' ) /= 0 ) THEN

     molecule_option = chain
     WRITE ( unit=output_unit, fmt='(a)' ) 'Chain of atoms, no periodic boundaries, no density'
     STOP 'Error in adjust_density'

  ELSE IF ( INDEX ( molecules, 'nonlinear' ) /= 0 ) THEN

     molecule_option = nonlinear
     WRITE ( unit=output_unit, fmt='(a)' ) 'Nonlinear molecules, periodic boundaries'

  ELSE IF ( INDEX ( molecules, 'linear') /= 0 ) THEN

     molecule_option = linear
     WRITE ( unit=output_unit, fmt='(a)' ) 'Linear molecules, periodic boundaries'

  ELSE IF ( INDEX ( molecules, 'atoms') /= 0 ) THEN

     molecule_option = atoms
     WRITE ( unit=output_unit, fmt='(a)' ) 'Atoms, periodic boundaries'

  ELSE

     WRITE ( unit=error_unit, fmt='(a,a)') 'Unrecognized molecules option: ', molecules
     STOP 'Error in adjust_density'

  END IF

  ! The read_cnf_atoms routine should work in all cases, just to get n and box
  CALL read_cnf_atoms ( filename, n, box )
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of particles',   n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Simulation box length', box
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',               density

  IF ( ABS(delta) < tol ) STOP 'No density change requested'

  IF ( delta+density < 0.0 ) STOP 'New requested density would be negative'
  
  SELECT CASE ( molecule_option )

  CASE ( nonlinear)
     ALLOCATE ( r(3,n), v(3,n), w(3,n), e(0:3,n) )

  CASE default
     ALLOCATE ( r(3,n), v(3,n), w(3,n), e(3,n) )

  END SELECT

  SELECT CASE ( molecule_option )

  CASE ( atoms )
     IF ( velocities ) THEN
        CALL read_cnf_atoms ( filename, n, box, r, v )
     ELSE
        CALL read_cnf_atoms ( filename, n, box, r )
     END IF

  CASE default
     IF ( velocities ) THEN
        CALL read_cnf_mols ( filename, n, box, r, e, v, w )
     ELSE
        CALL read_cnf_mols ( filename, n, box, r, e )
     END IF

  END SELECT

  scale   = ( density / (density+delta) )**(1.0/3.0)
  box     = box * scale
  r(:,:)  = r(:,:) * scale
  density = REAL(n) / box**3
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'New density', density

  SELECT CASE ( molecule_option )

  CASE ( atoms )

     IF ( velocities ) THEN
        CALL write_cnf_atoms ( filename, n, box, r, v )
     ELSE
        CALL write_cnf_atoms ( filename, n, box, r )
     END IF

  CASE default

     IF ( velocities ) THEN
        CALL write_cnf_mols ( filename, n, box, r, e, v, w )
     ELSE
        CALL write_cnf_mols ( filename, n, box, r, e )
     END IF

  END SELECT

  DEALLOCATE ( r, v, e, w )

END PROGRAM adjust_density

