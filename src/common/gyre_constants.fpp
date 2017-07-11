! Module   : gyre_constants
! Purpose  : physical constants & environment variables
!
! Copyright 2013-2017 Rich Townsend
!
! This file is part of GYRE. GYRE is free software: you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, version 3.
!
! GYRE is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

$include 'core.inc'

module gyre_constants

  ! Uses

  use core_kinds
  use core_system
  use core_constants, &
       G_GRAVITY_ => G_GRAVITY, &
       C_LIGHT_ => C_LIGHT, &
       A_RADIATION_ => A_RADIATION, &
       M_SUN_ => M_SUN, &
       R_SUN_ => R_SUN, &
       L_SUN_ => L_SUN

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  ! Physical constants (cgs)

  real(WP), save, protected :: G_GRAVITY = G_GRAVITY_     ! Gravitational constant
  real(WP), save, protected :: C_LIGHT = C_LIGHT_         ! Speed of light in vacuuo
  real(WP), save, protected :: A_RADIATION = A_RADIATION_ ! Radiation constant

  ! Astronomical constants (cgs)

  real(WP), save, protected :: M_SUN = M_SUN_ ! Solar mass
  real(WP), save, protected :: R_SUN = R_SUN_ ! Solar radius
  real(WP), save, protected :: L_SUN = L_SUN_ ! Solar luminosity

  ! Paths

  character(FILENAME_LEN), save, protected :: GYRE_DIR = ''

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
  public :: M_PROTON
  public :: K_BOLTZMANN
  
  public :: M_SUN
  public :: R_SUN
  public :: L_SUN

  public :: FILENAME_LEN
  public :: GYRE_DIR

  public :: read_constants
  public :: set_constant

  ! Procedures

contains

  subroutine read_constants (unit)

    integer, intent(in) :: unit

    integer                   :: n_cn
    character(:), allocatable :: gyre_dir_
    integer                   :: status

    namelist /constants/ G_GRAVITY, C_LIGHT, A_RADIATION, &
                         M_SUN, R_SUN, L_SUN, GYRE_DIR

    ! Count the number of constants namelists

    rewind(unit)

    n_cn = 0

    count_loop : do
       read(unit, NML=constants, END=100)
       n_cn = n_cn + 1
    end do count_loop

100 continue

    $ASSERT(n_cn == 1,Input file should contain exactly one &constants namelist)

    ! Read constants

    call get_env('GYRE_DIR', gyre_dir_, status)
    if (status == 0) GYRE_DIR = gyre_dir_
      
    rewind(unit)
    read(unit, NML=constants)

    $ASSERT(GYRE_DIR /= '',GYRE_DIR is not set)

    ! Finish

    return

  end subroutine read_constants

  subroutine set_constant (name, value)

    character(*), intent(in) :: name
    real(WP), intent(in)     :: value

    ! Set the constant

    select case (name)
    case ('G_GRAVITY')
       G_GRAVITY = value
    case ('C_LIGHT')
       C_LIGHT = value
    case ('A_RADIATION')
       A_RADIATION = value
    case ('M_SUN')
       M_SUN = value
    case ('R_SUN')
       R_SUN = value
    case ('L_SUN')
       L_SUN = value
    case default
       $ABORT(Invalid name)
    end select

    ! Finish

    return

  end subroutine set_constant

end module gyre_constants
