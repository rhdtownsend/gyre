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

  $if ($MPI)
  public :: bcast_constants
  $endif

  ! Procedures

contains

  $if ($MPI)

  subroutine bcast_constants (root_rank)

    integer, intent(in) :: root_rank

    ! Broadcast the mutable constants

    call bcast(G_GRAVITY, root_rank)
    call bcast(C_LIGHT, root_rank)
    call bcast(A_RADIATION, root_rank)

    call bcast(M_SUN, root_rank)
    call bcast(R_SUN, root_rank)
    call bcast(L_SUN, root_rank)

    call bcast(GYRE_DIR, root_rank)

    ! Finish

    return

  end subroutine bcast_constants

  $endif

end module gyre_constants
