! dpd_module.f90
! Dissipative particle dynamics module
MODULE dpd_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: n, r, v, f
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: force, make_ij, lowe, shardlow

  INTEGER                              :: n  ! Number of atoms
  REAL,    DIMENSION(:,:), ALLOCATABLE :: r  ! Positions (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: v  ! Velocities (3,n)
  REAL,    DIMENSION(:,:), ALLOCATABLE :: f  ! Forces (3,n)

  INTEGER                              :: nl ! List size
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ij ! List of ij pairs within range (2,nl)
  INTEGER                              :: np ! Number of pairs within range

  REAL, PARAMETER :: r_cut = 1.0  ! Cutoff distance (unit of length)

CONTAINS

  SUBROUTINE introduction ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)'           ) 'DPD soft potential'
    WRITE ( unit=output_unit, fmt='(a,t40,f15.5)' ) 'Diameter, r_cut = ', r_cut    

  END SUBROUTINE introduction

  SUBROUTINE conclusion ( output_unit )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: output_unit ! Unit for standard output

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays ( box )
    IMPLICIT NONE
    REAL, INTENT(in) :: box ! Simulation box length

    REAL, PARAMETER :: pi = 4.0*ATAN(1.0)

    ! Estimate pair list size, with 30% margin for error
    nl = CEILING ( 1.3*(4.0*pi/3.0)*REAL(n*(n-1)/2)*(r_cut/box)**3 )
    ALLOCATE ( r(3,n), v(3,n), f(3,n), ij(2,nl) )

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    IMPLICIT NONE

    DEALLOCATE ( r, v, f, ij )

  END SUBROUTINE deallocate_arrays

  SUBROUTINE make_ij ( box ) 
    USE maths_module, ONLY : random_integer
    IMPLICIT NONE
    REAL, INTENT(in) :: box ! Simulation box length

    ! Compiles randomized list of pairs within range, and stores in the array ij
    ! with np indicating the number of such pairs

    ! It is assumed that positions in the array r are in units where box = 1

    INTEGER            :: p, q, i, j
    REAL               :: rij_sq
    REAL, DIMENSION(3) :: rij

    p = 0
    DO i = 1, n-1 ! Outer loop over atoms
       DO j = i+1, n ! Inner loop over atoms

          rij(:)  = r(:,i) - r(:,j)
          rij(:) = rij(:) - ANINT ( rij(:) )
          rij(:) = rij(:) * box ! Now in r_cut=1 units
          rij_sq = SUM ( rij(:)**2 )

          IF ( rij_sq < 1.0 ) THEN ! Test for centre-centre spherical cutoff
             p = p + 1
             IF ( p > nl ) THEN
                WRITE ( unit=error_unit, fmt='(a,2i15)') 'Pair list error', p, nl
                STOP 'List error in make_ij'
             END IF
             ij(1,p) = i
             ij(2,p) = j
          END IF ! End test for centre-centre spherical cutoff
       END DO ! End inner loop over atoms
    END DO ! End outer loop over atoms

    np = p ! Store number of pairs

    !  Construct random permutation
    DO p = 1, np-1
       q = random_integer ( p, np )            ! Random integer in [p,np] inclusive
       IF ( p /= q ) ij(:,[p,q]) = ij(:,[q,p]) ! Swap directly using slices
    END DO

  END SUBROUTINE make_ij

  SUBROUTINE force ( box, a, pot, vir, lap )
    IMPLICIT NONE
    REAL, INTENT(in)  :: box ! Simulation box length
    REAL, INTENT(in)  :: a   ! Force strength parameter
    REAL, INTENT(out) :: pot ! Total potential energy
    REAL, INTENT(out) :: vir ! Virial
    REAL, INTENT(out) :: lap ! Laplacian

    ! Calculates potential (pot), virial (vir), Laplacian (lap)
    ! Also calculates forces and stores them in the array f

    ! It is assumed that positions in the array r are in units where box = 1
    ! and that the array ij contains a list of all np pairs within range

    ! Forces are calculated in units where r_cut = 1

    INTEGER            :: p, i, j
    REAL               :: rij_sq, rij_mag, wij
    REAL, DIMENSION(3) :: rij, fij, rij_hat

    f   = 0.0
    pot = 0.0
    vir = 0.0
    lap = 0.0

    DO p = 1, np ! Loop over all pairs within range
       i = ij(1,p)
       j = ij(2,p)

       rij(:)  = r(:,i) - r(:,j)           ! separation vector
       rij(:)  = rij(:) - ANINT ( rij(:) ) ! periodic boundary conditions in box=1 units
       rij(:)  = rij(:) * box              ! now in r_cut=1 units
       rij_sq  = SUM ( rij**2 )            ! squared separation
       IF ( rij_sq > 1.0 ) CYCLE           ! this should never happen
       rij_mag = SQRT ( rij_sq )           ! separation distance
       wij     = 1.0 - rij_mag             ! weight function
       rij_hat = rij / rij_mag             ! unit separation vector
       pot     = pot + 0.5 * wij**2        ! potential
       vir     = vir + wij * rij_mag       ! virial
       lap     = lap + (3.0-2.0/rij_mag)   ! Laplacian
       fij     = wij * rij_hat(:)
       f(:,i)  = f(:,i) + fij
       f(:,j)  = f(:,j) - fij

    END DO ! End loop over all pairs within range

    ! Multiply results by numerical factors
    pot = pot * a
    vir = vir * a / 3.0
    lap = lap * a * 2.0
    f   = f   * a

  END SUBROUTINE force

  SUBROUTINE lowe ( box, temperature, gamma_step )
    USE maths_module, ONLY : random_normal
    IMPLICIT NONE
    REAL, INTENT(in) :: box         ! Simulation box length
    REAL, INTENT(in) :: temperature ! Specified temperature
    REAL, INTENT(in) :: gamma_step  ! Pair selection probability = Gamma * timestep

    ! Apply pairwise Lowe-Andersen thermostat to velocities stored in array v

    ! It is assumed that positions in the array r are in units where box = 1
    ! and that the array ij contains a list of all np pairs within range

    INTEGER            :: i, j, p
    REAL               :: rij_sq, rij_mag, zeta, v_par
    REAL, DIMENSION(3) :: rij, vij, rij_hat

    DO p = 1, np ! Loop over all pairs within range
       CALL RANDOM_NUMBER ( zeta )
       IF ( zeta < gamma_step ) THEN ! Select pair with desired probability
          i = ij(1,p)
          j = ij(2,p)

          ! Choose Gaussian random number with desired standard deviation
          zeta    = random_normal ( 0.0, SQRT(2.0*temperature) )

          rij(:)  = r(:,i) - r(:,j)
          rij(:)  = rij(:) - ANINT ( rij(:) )
          rij(:)  = rij(:) * box ! now in r_cut=1 units
          rij_sq  = SUM ( rij(:)**2 )
          rij_mag = SQRT ( rij_sq )
          rij_hat = rij / rij_mag
          vij     = v(:,i) - v(:,j)
          v_par   = DOT_PRODUCT(vij,rij_hat)
          vij     = 0.5 * ( zeta - v_par ) * rij_hat(:) ! change in velocity
          v(:,i) = v(:,i) + vij
          v(:,j) = v(:,j) - vij

       END IF ! End select pair with desired probability
    END DO ! End loop over all pairs within range

  END SUBROUTINE lowe

  SUBROUTINE shardlow ( box, temperature, gamma_step )
    USE maths_module, ONLY : random_normal
    IMPLICIT NONE
    REAL, INTENT(in) :: box         ! Simulation box length
    REAL, INTENT(in) :: temperature ! Specified temperature
    REAL, INTENT(in) :: gamma_step  ! Gamma * timestep

    ! Implements the Shardlow integration algorithm for velocities stored in array v

    ! It is assumed that positions in the array r are in units where box = 1
    ! and that the array ij contains a list of all np pairs within range

    INTEGER            :: i, j, p
    REAL               :: rij_sq, rij_mag, sqrt_gamma_step, prob, sqrt_prob, zeta, v_par, wij
    REAL, DIMENSION(3) :: rij, vij, rij_hat

    sqrt_gamma_step = SQRT(gamma_step)

    DO p = 1, np ! Loop over all pairs within range
       i = ij(1,p)
       j = ij(2,p)

       zeta  = random_normal ( 0.0, SQRT(2.0*temperature) ) ! normal with desired standard deviation

       rij(:)    = r(:,i) - r(:,j)
       rij(:)    = rij(:) - ANINT ( rij(:) )
       rij(:)    = rij(:) * box ! now in r_cut=1 units
       rij_sq    = SUM ( rij(:)**2 )
       rij_mag   = SQRT ( rij_sq )
       rij_hat   = rij / rij_mag
       wij       = 1.0 - rij_mag         ! weight function
       sqrt_prob = sqrt_gamma_step * wij ! sqrt of p-factor
       prob      = sqrt_prob**2          ! p-factor

       ! First half step
       vij(:) = v(:,i) - v(:,j)
       v_par  = DOT_PRODUCT ( vij(:), rij_hat(:) )
       vij(:) = 0.5 * ( sqrt_prob*zeta - prob*v_par ) * rij_hat(:) ! change in velocity
       v(:,i) = v(:,i) + vij(:)
       v(:,j) = v(:,j) - vij(:)

       ! Second half step
       vij(:) = v(:,i) - v(:,j)
       v_par  = DOT_PRODUCT ( vij, rij_hat )
       vij(:) = 0.5 * ( sqrt_prob*zeta - prob*v_par ) * rij_hat(:) / (1.0+prob) ! change in velocity
       v(:,i) = v(:,i) + vij(:)
       v(:,j) = v(:,j) - vij(:)

    END DO ! End loop over all pairs within range

  END SUBROUTINE shardlow

END MODULE dpd_module
