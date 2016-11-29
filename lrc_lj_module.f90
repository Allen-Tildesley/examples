! lrc_lj_module.f90
! Long-range and delta corrections for potential energy and pressure
MODULE lrc_module

  ! The purpose of this module is simply to gather in one place the common
  ! functions for long-range and delta corrections for the Lennard-Jones potential.
  ! If a different potential is used in the simulation, a different file
  ! (with the same module name) containing different expressions should be substituted

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: potential_lrc, pressure_lrc, pressure_delta

  ! Private data
  REAL, PARAMETER :: pi = 4.0 * ATAN(1.0)

CONTAINS

  FUNCTION potential_lrc ( density, r_cut )
    IMPLICIT NONE
    REAL                :: potential_lrc ! Returns long-range correction to potential/atom
    REAL,    INTENT(in) :: density       ! Number density N/V
    REAL,    INTENT(in) :: r_cut         ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones potential per atom
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL :: sr3

    sr3 = 1.0 / r_cut**3

    potential_lrc = pi * ( (8.0/9.0)  * sr3**3  - (8.0/3.0)  * sr3 ) * density

  END FUNCTION potential_lrc

  FUNCTION pressure_lrc ( density, r_cut )
    IMPLICIT NONE
    REAL                :: pressure_lrc ! Returns long-range correction to pressure
    REAL,    INTENT(in) :: density      ! Number density N/V
    REAL,    INTENT(in) :: r_cut        ! Cutoff distance

    ! Calculates long-range correction for Lennard-Jones pressure
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL :: sr3

    sr3 = 1.0 / r_cut**3

    pressure_lrc = pi * ( (32.0/9.0) * sr3**3  - (16.0/3.0) * sr3 ) * density**2

  END FUNCTION pressure_lrc

  FUNCTION pressure_delta ( density, r_cut )
    IMPLICIT NONE
    REAL                :: pressure_delta ! Returns delta correction to pressure
    REAL,    INTENT(in) :: density        ! Number density N/V
    REAL,    INTENT(in) :: r_cut          ! Cutoff distance

    ! Calculates correction for Lennard-Jones pressure
    ! due to discontinuity in the potential at r_cut
    ! density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1

    REAL :: sr3

    sr3 = 1.0 / r_cut**3

    pressure_delta = pi * (8.0/3.0) * ( sr3**3  - sr3 ) * density**2

  END FUNCTION pressure_delta

END MODULE lrc_module
