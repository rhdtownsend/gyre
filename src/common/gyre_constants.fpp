! Module   : gyre_constants
! Purpose  : (mutable) physical constants

module gyre_constants

  ! Uses

  use core_kinds
  use core_parallel
  use core_constants, &
       G_GRAVITY_ => G_GRAVITY, &
       C_LIGHT_ => C_LIGHT, &
       A_RADIATION_ => A_RADIATION, &
       M_SUN_ => M_SUN, &
       R_SUN_ => R_SUN, &
       L_SUN_ => L_SUN

  ! No implicit typing

  implicit none

  ! Module variables

  ! Physical constants (cgs)

  real(WP), save :: G_GRAVITY = G_GRAVITY_     ! Gravitational constant
  real(WP), save :: C_LIGHT = C_LIGHT_         ! Speed of light in vacuuo
  real(WP), save :: A_RADIATION = A_RADIATION_ ! Radiation constant

  ! Astronomical constants (cgs)

  real(WP), save :: M_SUN = M_SUN_ ! Solar mass
  real(WP), save :: R_SUN = R_SUN_ ! Solar radius
  real(WP), save :: L_SUN = L_SUN_ ! Solar luminosity (Allen, 1973)

  ! Paths

  character(FILENAME_LEN), save :: GYRE_DIR = ''

  ! Access specifiers

  private

  public :: PI
  public :: TWOPI
  public :: HALFPI
  public :: DEG_TO_RAD
  public :: RAD_TO_DEG

  public :: G_GRAVITY
  public :: C_LIGHT
  public :: A_RADIATION
  
  public :: M_SUN
  public :: R_SUN
  public :: L_SUN

  public :: FILENAME_LEN
  public :: GYRE_DIR

end module gyre_constants
