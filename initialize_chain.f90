! initialize_chain.f90
! Sets up initial configuration for MD or MC of chain molecule
PROGRAM initialize_chain

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <dominic.tildesley@epfl.ch> ("the authors"),                          !
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

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  USE, INTRINSIC :: iso_fortran_env,         ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module,        ONLY : write_cnf_atoms
  USE               initialize_chain_module, ONLY : allocate_arrays, deallocate_arrays, &
       &                                            initialize_random, initialize_velocities, n, r, v

  IMPLICIT NONE

  ! Most important variables
  REAL    :: temperature ! Specified temperature
  LOGICAL :: velocities  ! User option requiring velocities for MD (otherwise just positions)
  REAL    :: bond        ! Bond length
  LOGICAL :: soft        ! Flag for soft interactions (i.e. ignore overlaps)

  INTEGER            :: ioerr

  REAL, PARAMETER :: tol = 1.0e-6

  CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp' ! Will be used as an input file by later simulations

  NAMELIST /nml/ n, temperature, bond, velocities, soft

  WRITE ( unit=output_unit, fmt='(a)' ) 'initialize_chain'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Sets up initial configuration file for chain simulations'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass m=1 throughout'
  WRITE ( unit=output_unit, fmt='(a)' ) 'NO periodic boundaries'

  CALL RANDOM_SEED()

  ! Set default parameters
  n                   = 13
  temperature         = 1.0
  bond                = 1.0
  velocities          = .FALSE. ! By default, produce positions only, for MC simulations
  soft                = .FALSE. ! By default, check for overlaps when placing molecules

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in initialize'
  END IF

  WRITE ( unit=output_unit, fmt='(a,t40,i15  )' ) 'n = ',        n
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Bond length', bond
  IF ( soft ) WRITE ( unit=output_unit, fmt='(a)' ) 'Soft option selected - no overlap checking'

  CALL allocate_arrays

  ! Initialize positions and optionally velocities

  CALL initialize_random ( bond, soft ) ! Random with chosen bond length
  IF ( velocities ) CALL initialize_velocities ( temperature )

  ! Write out configuration

  WRITE ( unit=output_unit, fmt='(a,a)' ) 'Writing configuration to filename ', filename

  ! We do not use periodic boundaries for this system
  ! Instead, use "box" variable to store bond length
  IF ( velocities ) THEN
     CALL write_cnf_atoms ( filename, n, bond, r, v )
  ELSE
     CALL write_cnf_atoms ( filename, n, bond, r )
  END IF

  CALL deallocate_arrays

END PROGRAM initialize_chain
