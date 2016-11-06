! eos_lj_kn_module.f90
! Equation of state for Lennard-Jones: Kolafa-Nezbeda prescription
MODULE eos_lj_module

  IMPLICIT NONE
  PRIVATE
  
  ! Module providing the thermodynamic properties of the Lennard-Jones fluid
  ! Expressions based on
  ! J Kolafa, I Nezbeda, Fluid Phase Equilibria 100, 1 (1994)
  ! henceforth referred to as K&N

  ! Functions are provided to give:
  ! A/N the Helmholtz free energy per atom, alj
  ! G/N=mu, the Gibbs free energy per atom or chemical potential, glj,
  ! P, the pressure, plj,
  ! U/N, the internal energy, ulj
  ! for given temperature t and density rho
  ! Lennard-Jones reduced units (epsilon=1, sigma=1) are used throughout

  ! Public routines
  PUBLIC :: alj, plj, ulj, glj

  ! Private data
  
  REAL, PARAMETER :: pi = 4.0*ATAN(1.0), pi6 = pi/6.0

  REAL, PARAMETER :: gamma = 1.92907278 ! Damping parameter

  ! Coefficients ai(j) in K&N Table 3 (called Cij by them),
  ! We tabulate these separately for each value of i which indexes sqrt(T)
  ! Each ai array corresponds to (1/sqrt(T))**i i.e. to powers 0, -1, -2, and -4 of sqrt(T) in the table
  ! (there are no terms with power -3, hence no a3 array)
  ! Within each array, ai(j) corresponds to rho**(2+j), i.e. to powers j = 2, 3, 4, 5, 6 of rho in the table
  ! and we include the rho**2 factor afterwards
  REAL, DIMENSION(0:4), PARAMETER :: a0 = [   2.01546797,  -28.17881636,   28.28313847,  -10.42402873,   0.0        ]
  REAL, DIMENSION(0:4), PARAMETER :: a1 = [ -19.58371655,   75.62340289, -120.70586598,   93.92740328, -27.37737354 ]
  REAL, DIMENSION(0:4), PARAMETER :: a2 = [  29.34470520, -112.35356937,  170.64908980, -123.06669187,  34.42288969 ]
  REAL, DIMENSION(0:4), PARAMETER :: a4 = [ -13.37031968,   65.38059570, -115.09233113,   88.91973082, -25.62099890 ]

  ! Coefficients d(i) for diameter (hybrid Barker-Henderson hBH), or diam for short,
  ! given in K&N Table 2, left column
  ! d(i) corresponds to sqrt(T)**(-2+i) i.e. to powers -2, -1, 0, 1 of sqrt(T) in the table
  ! and we include the 1/T factor afterwards
  ! There is also a logarithmic term with coefficient dln
  REAL,                 PARAMETER :: dln = -0.063920968
  REAL, DIMENSION(0:3), PARAMETER :: d = [ 0.011117524, -0.076383859, 1.080142248, 0.000693129 ]

  ! Coefficients b(i) for Delta B2 (hybrid Barker-Henderson hBH) or bdel for short,
  ! given in K&N Table 2, right column
  ! b(i) corresponds to (1/sqrt(T))**i i.e. to powers 0, -1, -2, -3, -4, -5, -6, -7 of sqrt(T) in the table
  REAL, DIMENSION(0:7), PARAMETER :: b = [ 0.02459877, 0.0,       -7.02181962,  2.90616279, &
       &                                  -4.13749995, 0.87361369, 0.43102052, -0.58544978 ]

CONTAINS

  FUNCTION alj ( t, rho )
    IMPLICIT NONE
    REAL             :: alj ! Returns LJ Helmholtz free energy per atom, for given
    REAL, INTENT(in) :: t   ! temperature, and
    REAL, INTENT(in) :: rho ! density

    ! NB A/N (ideal gas) includes an additive term 3*kT*log(Lambda)
    ! where Lambda = h*sqrt(2*pi/m*kT) and h is Planck's constant
    ! which in turn gives rise to the ideal gas part of U, 1.5*kT
    ! This T-dependent term should be considered whenever the absolute
    ! values of free energies, entropies, and chemical potentials are discussed
    ! Otherwise, use excess A/N, S/N, G/N=mu relative to ideal gas values
    ! We do not include this term here, effectively setting Lambda=1 for any T
    
    REAL :: eta ! Effective packing fraction

    ! Refer to K&N Eq (28) or (30)

    ! Hard-sphere residual contribution to A/N
    eta = pi6 * rho * diam(t)**3
    alj = ahs ( eta ) * t

    ! Damped perturbed virial expansion term
    alj = alj + EXP(-gamma*rho**2) * rho * t * bdel(t) 

    ! Extra residual term (double power series in rho and sqrt(t))
    ! Expansion in rho starts with rho**2 term
    alj = alj + polynomial ( rho, a0 ) * rho**2
    alj = alj + polynomial ( rho, a1 ) * rho**2 / SQRT(t)
    alj = alj + polynomial ( rho, a2 ) * rho**2 / t
    alj = alj + polynomial ( rho, a4 ) * rho**2 / t**2

    alj = alj + t * ( LOG(rho) - 1.0 ) ! Add ideal gas Helmholtz free energy

  END FUNCTION alj

  FUNCTION glj ( t, rho )
        IMPLICIT NONE
    REAL             :: glj ! Returns LJ chemical potential, for given
    REAL, INTENT(in) :: t   ! temperature, and
    REAL, INTENT(in) :: rho ! density

    ! NB G/N (ideal gas), like A/N, includes an additive term 3*kT*log(Lambda)
    ! We do not include this term here (see note above)

    glj = alj ( t, rho ) + plj ( t, rho ) / rho
    
  END FUNCTION glj
  
  FUNCTION plj ( t, rho )
    IMPLICIT NONE
    REAL             :: plj ! Returns LJ pressure, for given
    REAL, INTENT(in) :: t   ! temperature, and
    REAL, INTENT(in) :: rho ! density

    ! NB in this routine the ideal gas part is included
    
    REAL :: eta ! Effective packing fraction

    ! The coefficients are f(j)*ai(j) f(j)=j+2 (for each power j of rho, starting at rho**2)
    REAL, DIMENSION(0:4), PARAMETER :: f = [2.0,3.0,4.0,5.0,6.0] ! Coefficient factors

    ! Refer to K&N Eq (31)

    ! Hard sphere contribution to P/(rho*T) including ideal-gas part
    eta = pi6 * rho * diam(t)**3
    plj = zhs ( eta )

    ! Damped perturbed virial expansion term for P/(rho*T)
    plj = plj + rho * (1.0-2.0*gamma*rho**2) * EXP(-gamma*rho**2) * bdel(t) 

    plj = plj * t ! Now we are dealing with P/rho

    ! Extra residual term (double power series in rho and sqrt(t))
    ! The coefficients are f(j)*ai(j) for each power j of rho, see above
    ! Expansion in rho starts with rho**2 term
    plj = plj + polynomial ( rho, f*a0 ) * rho**2
    plj = plj + polynomial ( rho, f*a1 ) * rho**2 / SQRT(t)
    plj = plj + polynomial ( rho, f*a2 ) * rho**2 / t
    plj = plj + polynomial ( rho, f*a4 ) * rho**2 / t**2

    plj = plj * rho ! Finally we have P (including ideal-gas part)

  END FUNCTION plj

  FUNCTION ulj ( t, rho)
    IMPLICIT NONE
    REAL             :: ulj ! Returns LJ internal energy per atom, for given
    REAL, INTENT(in) :: t   ! temperature, and
    REAL, INTENT(in) :: rho ! density

    ! NB in this routine we include the ideal gas part
    
    REAL :: eta ! Effective packing fraction

    ! Refer to K&N Eq (32)

    ! Hard sphere contribution to U/N
    eta = pi6 * rho * diam(t)**3
    ulj = 3.0 * (zhs(eta)-1.0) * diam_deriv(t) / diam(t) 

    ! Damped perturbed virial expansion
    ulj = ulj + rho * EXP(-gamma*rho**2) * bdel_deriv(t) 

    ! Extra residual term (double power series in rho and sqrt(t))
    ! The coefficients are (1-i/2)*ai(j) for each power i of sqrt(T)
    ! Powers are 0, -1, -2, -4 hence factors below of 1, 1.5, 2.0, 3.0
    ! Expansion in rho starts with rho**2 term
    ulj = ulj + polynomial ( rho, a0 ) * rho**2 
    ulj = ulj + polynomial ( rho, a1 ) * rho**2 * 1.5 / SQRT(t)
    ulj = ulj + polynomial ( rho, a2 ) * rho**2 * 2.0 / t  
    ulj = ulj + polynomial ( rho, a4 ) * rho**2 * 3.0 / t**2

    ulj = ulj + 1.5*t ! Add ideal gas part
    
  END FUNCTION ulj

  FUNCTION zhs ( eta )
    IMPLICIT NONE
    REAL             :: zhs ! Returns HS P/(rho*kT), for
    REAL, INTENT(in) :: eta ! given packing fraction

    ! Numerator coefficients in K&N Eq (4)
    REAL, DIMENSION(0:4), PARAMETER :: c = [ 1.0, 1.0, 1.0, -2.0/3.0, -2.0/3.0 ]

    ! Refer to K&N Eq (4)

    ! Extended Carnahan-Starling formula for hard spheres (includes ideal-gas part)
    zhs = polynomial ( eta, c ) / (1.0-eta)**3 ! P/(rho*kT)

  END FUNCTION zhs

  FUNCTION ahs ( eta )
    IMPLICIT NONE
    REAL             :: ahs ! Returns excess HS Helmholtz free energy divided by T, A(res)/NkT, for
    REAL, INTENT(in) :: eta ! given packing fraction

    ! NB this is the excess A/N for hard spheres, not including ideal gas part
    
    ! Numerator coefficients in Eq (5)
    REAL, DIMENSION(0:2), PARAMETER :: c = [ 34.0, -33.0, 4.0 ] 

    ! Refer to K&N Eq (5)

    ! Extended Carnahan-Starling formula for hard spheres
    ahs = polynomial ( eta, c ) * eta / 6.0
    ahs = ahs / (1.0-eta)**2
    ahs = ahs + (5.0/3.0) * LOG(1.0-eta) ! Logarithm part of A(res)/NkT

  END FUNCTION ahs

  FUNCTION diam(t)
    IMPLICIT NONE
    REAL             :: diam ! Returns diameter (hybrid Barker-Henderson), for
    REAL, INTENT(in) :: t    ! given temperature

    REAL :: s

    ! Refer to K&N Eq (29)

    s = SQRT(t) ! Expand in powers of sqrt(t)

    ! Coefficients given in K&N Table 2, see module data
    ! Series for diam starts with sqrt(T)**(-2) term, i.e. T**(-1) term
    diam = polynomial ( s, d ) / t ! Include T**(-1) here
    diam = diam + dln * LOG(t)     ! Add logarithmic term

  END FUNCTION diam

  FUNCTION diam_deriv ( t )
    IMPLICIT NONE
    REAL             :: diam_deriv ! Returns derivative of diam with respect to (1/T), for
    REAL, INTENT(in) :: t          ! given temperature

    ! Coefficients for diam, d, given in K&N Table 2, see module data
    ! Derived coefficient factors from K&N Eq (33) represented as f(i)*d(i)
    ! Powers of sqrt(T) are i = -2,  -1, 0,  1    with coefficients stored in d(i+2)
    ! so factors are   (-i/2) =  1, 1/2, 0, -1/2, and they are stored in f(i+2)
    REAL, DIMENSION(0:3), PARAMETER :: f = [ 1.0, 0.5, 0.0, -0.5 ]

    REAL :: s

    s = SQRT(t) ! Expand in powers of sqrt(T)

    ! Original series for diam starts with T**(-1) term, but Eq (33) introduces a factor T
    diam_deriv = polynomial ( s, f*d )
    diam_deriv = diam_deriv - dln * t ! Add derivative of log term here

  END FUNCTION diam_deriv

  FUNCTION bdel ( t )
    IMPLICIT NONE
    REAL             :: bdel ! Returns Delta B2 (hybrid Barker-Henderson) virial expansion term, for
    REAL, INTENT(in) :: t    ! given temperature

    REAL :: s

    ! Refer to K&N Eq (29)

    s = 1.0 / SQRT(t) ! Expand in powers of inverse sqrt(T)

    ! Coefficients given in K&N Table 2, see module data
    bdel = polynomial ( s, b ) ! Polynomial with given coefficients

  END FUNCTION bdel

  FUNCTION bdel_deriv ( t )
    IMPLICIT NONE
    REAL             :: bdel_deriv ! Returns derivative of bdel with respect to (1/T), for
    REAL, INTENT(in) :: t          ! given temperature

    ! Coefficients for bdel, b, given in K&N Table 2, see module data
    ! Derived coefficient factors from Eq (33)
    ! Powers of sqrt(T) are i = 0,  -1,  -2,  -3,  -4,  -5,  -6,  -7, with coefficients stored in b(-i)
    ! so factors are   (-i/2) = 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, and they are stored in f(-i)
    ! however b(0) is multiplied by zero, and b(1) is already zero,
    ! so our interest starts with b(2) and we store just the factors f(2) .. f(7)
    REAL, DIMENSION(2:7), PARAMETER :: f = [ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5 ]

    REAL :: s

    s = 1.0 / SQRT(t) ! Expand in powers of inverse sqrt(T)

    ! The b(2) term in the original series corresponds to T**(-1)
    ! Eq (33) introduces a factor T, so the derivative expansion begins with the constant term
    bdel_deriv = polynomial ( s, f(2:7)*b(2:7) )

  END FUNCTION bdel_deriv

  FUNCTION polynomial ( x, c ) RESULT ( p )
    IMPLICIT NONE
    REAL                           :: p ! Returns polynomial in
    REAL,               INTENT(in) :: x ! supplied variable with coefficients in
    REAL, DIMENSION(:), INTENT(in) :: c ! supplied array

    ! Horner's method
    ! The first coefficient in c is the constant term
    ! The second is the coefficient of x
    ! The third is the coefficient of x**2 and so on

    INTEGER :: i

    p = 0.0
    DO i = SIZE(c), 1, -1 ! Loop over decreasing powers
       p = p * x + c(i)
    END DO ! End loop over decreasing powers

  END FUNCTION polynomial

END MODULE eos_lj_module
