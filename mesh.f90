! mesh.f90
! Assignment of charges to a 3-d mesh
PROGRAM mesh

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

  ! This program assigns a set of charges to a cubic mesh using the
  ! triangular-shaped cloud distribution described by Hockney and Eastwood (1988)
  ! The charges are positioned in a box of unit length.
  ! The charge mesh is indexed from 0 to sc-1 in each coordinate direction

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  USE mesh_module, ONLY : mesh_function

  IMPLICIT NONE

  INTEGER                             :: n   ! Number of charges
  INTEGER                             :: sc  ! Dimension of mesh
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rho ! Mesh cell charge density (0:sc-1,0:sc-1,0:sc-1) 
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: r   ! Charge positions (3,n) 
  REAL, ALLOCATABLE, DIMENSION(:)     :: q   ! Charges (n) 

  INTEGER :: n2, n3, ioerr
  REAL    :: h

  NAMELIST /nml/ n, sc

  WRITE ( unit=output_unit, fmt='(a)' ) 'mesh'
  WRITE ( unit=output_unit, fmt='(a)' ) '3-D mesh assignment of charges'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Unit box length, coordinates in range (0,1)'

  ! Example values for illustration
  n  = 4  ! Number of charges
  sc = 10 ! Dimension of mesh

  ! Read parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  READ ( unit=input_unit, nml=nml, iostat=ioerr )
  IF ( ioerr /= 0 ) THEN
     WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
     IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     STOP 'Error in mesh'
  END IF

  h = 1.0 / REAL( sc ) ! Mesh spacing

  ! Write out run parameters
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Number of charges', n
  WRITE ( unit=output_unit, fmt='(a,t40,i15)'   ) 'Dimension of mesh', sc
  WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Mesh spacing',      h

  ALLOCATE ( rho(0:sc-1,0:sc-1,0:sc-1) ) ! C-style indexing is convenient here
  ALLOCATE ( r(3,n), q(n) )

  ! For illustration we choose random charge positions with coordinates in range (0,1)
  ! In a real application, we would convert positions into this range
  CALL RANDOM_SEED() ! same random number sequence every time
  CALL RANDOM_NUMBER ( r )

  ! For illustration we choose +1 and -1 charges, alternately
  q(1::2) = 1.0
  q(2::2) = -1.0

  rho = mesh_function ( r, q, sc )

  ! Output mesh charge density
  DO n3 = 0, sc-1
     WRITE( unit=output_unit, fmt='(a,i5)' ) 'z-layer ', n3
     DO n2 = 0, sc-1
        WRITE( unit=output_unit, fmt='(*(f10.4))') rho(:,n2,n3)
     END DO
  END DO

  ! Finally check integrated charge density
  WRITE( unit=output_unit, fmt='(a,2f10.6)') 'Total charge = ', SUM ( q ), SUM ( rho ) * (h**3)

  DEALLOCATE ( r, q, rho )
  
END PROGRAM mesh
