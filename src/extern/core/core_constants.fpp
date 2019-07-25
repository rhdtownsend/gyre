! Module   : core_constants
! Purpose  : physical constants

module core_constants

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Parameter definitions

  ! Mathematical constants

  real(WP), parameter :: PI = ACOS(-1._WP)
  real(WP), parameter :: TWOPI = 2._WP*PI
  real(WP), parameter :: HALFPI = ASIN(1._WP)

  real(WP), parameter :: DEG_TO_RAD = PI/180._WP
  real(WP), parameter :: RAD_TO_DEG = 1._WP/DEG_TO_RAD

  ! Physical constants (cgs)

  real(WP), parameter :: G_GRAVITY = 6.67428E-8_WP                ! Gravitational constant
  real(WP), parameter :: H_PLANCK =  6.62606896E-27_WP            ! Planck's constant
  real(WP), parameter :: K_BOLTZMANN = 1.3806504E-16_WP           ! Boltzmann's constant
  real(WP), parameter :: C_LIGHT = 2.99792458E10_WP               ! Speed of light in vacuuo
  real(WP), parameter :: SIGMA_STEFAN = 5.670400E-5_WP            ! Stefan's constant
  real(WP), parameter :: A_RADIATION = 4._WP*SIGMA_STEFAN/C_LIGHT ! Radiation constant
  real(WP), parameter :: U_ATOMIC = 1.660538782E-24_WP            ! Atomic mass unit
  real(WP), parameter :: ELECTRON_VOLT = 1.602176487E-12_WP       ! Electron volt
  real(WP), parameter :: E_ELECTRON = 4.80320427E-10_WP           ! Electron charge
  real(WP), parameter :: M_ELECTRON = 9.10938215E-28_WP           ! Electron mass
  real(WP), parameter :: M_PROTON = 1.672621898E-24_WP            ! Proton mass
  real(WP), parameter :: SIGMA_THOMSON = 6.6524586E-25_WP         ! Thomson cross section
  real(WP), parameter :: N_AVOGADRO = 6.0221415E23_WP             ! Avogadro's number
  real(WP), parameter :: SEC_YEAR = 24._WP*365.25_WP*3600._WP     ! Seconds in a year
  
  ! Astronomical constants (cgs)

  real(WP), parameter :: M_SUN = 1.9891E33_WP             ! Solar mass
  real(WP), parameter :: R_SUN = 6.96E10_WP               ! Solar radius
  real(WP), parameter :: L_SUN = 3.826E33_WP              ! Solar luminosity (Allen, 1973)
  real(WP), parameter :: YEAR = 365.25_WP*24._WP*3600._WP ! Solar Year
  real(WP), parameter :: PARSEC = 3.0856776E18_WP         ! Parsec

  ! Filename lengths etc

  integer, parameter :: FILENAME_LEN = 256

  ! Acess specifiers

  private

  public :: PI
  public :: TWOPI
  public :: HALFPI
  public :: DEG_TO_RAD
  public :: RAD_TO_DEG

  public :: G_GRAVITY
  public :: H_PLANCK
  public :: K_BOLTZMANN
  public :: C_LIGHT
  public :: SIGMA_STEFAN
  public :: A_RADIATION
  public :: U_ATOMIC
  public :: ELECTRON_VOLT
  public :: E_ELECTRON
  public :: M_ELECTRON
  public :: M_PROTON
  public :: SIGMA_THOMSON
  public :: N_AVOGADRO
  public :: SEC_YEAR
  
  public :: M_SUN
  public :: R_SUN
  public :: L_SUN
  public :: YEAR
  public :: PARSEC

  public :: FILENAME_LEN

end module core_constants
