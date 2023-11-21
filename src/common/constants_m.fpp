! Module  : constants_m
! Purpose : physical constants & environment variables
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

module constants_m

  ! Uses

  use kinds_m
  use system_m

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  ! Physical constants (cgs)

  real(WP), parameter :: DEFAULT_G_GRAVITY = 6.67430e-8_WP                                ! Gravitational constant (CODATA 2018)
  real(WP), parameter :: DEFAULT_C_LIGHT = 2.99792458E10_WP                               ! Speed of light in vacuuo (CODATA 2018)
  real(WP), parameter :: DEFAULT_SIGMA_STEFAN = 5.670374419E-5_WP                         ! Stefan's constant (CODATA 2018)
  real(WP), parameter :: DEFAULT_A_RADIATION = 4._WP*DEFAULT_SIGMA_STEFAN/DEFAULT_C_LIGHT ! Radiation constant

  real(WP), save, protected :: G_GRAVITY = DEFAULT_G_GRAVITY
  real(WP), save, protected :: C_LIGHT = DEFAULT_C_LIGHT
  real(WP), save, protected :: A_RADIATION = DEFAULT_A_RADIATION

  ! Astronomical constants (cgs)

  real(WP), parameter :: DEFAULT_GM_SUN = 1.3271244E26_WP                 ! Gravitational constant * nominal solar mass (IAU 215 resolution B3)
  real(WP), parameter :: DEFAULT_M_SUN = DEFAULT_GM_SUN/DEFAULT_G_GRAVITY ! Nominal Solar mass
  real(WP), parameter :: DEFAULT_R_SUN = 6.957E10_WP                      ! Nominal solar radius (IAU 2015 resolution B3)
  real(WP), parameter :: DEFAULT_L_SUN = 3.828E33_WP                      ! Nominal solar luminosity (IAU 2015 resolution B3)
  
  real(WP), save, protected :: M_SUN = DEFAULT_M_SUN
  real(WP), save, protected :: R_SUN = DEFAULT_R_SUN
  real(WP), save, protected :: L_SUN = DEFAULT_L_SUN

  ! Lengths

  integer, parameter :: FILENAME_LEN = 256
  integer, parameter :: ITEM_LEN = 32

  ! Paths

  character(FILENAME_LEN), save, protected :: GYRE_DIR = ''

  ! Interfaces

  interface set_constant
     module procedure set_constant_r_
     module procedure set_constant_c_
  end interface set_constant

  ! Access specifiers

  private

  public :: G_GRAVITY
  public :: C_LIGHT
  public :: A_RADIATION
  
  public :: M_SUN
  public :: R_SUN
  public :: L_SUN

  public :: ITEM_LEN
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

  !****

  subroutine set_constant_r_ (name, value)

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

  end subroutine set_constant_r_

  !****

  subroutine set_constant_c_ (name, value)

    character(*), intent(in) :: name
    character(*), intent(in) :: value

    ! Set the constant

    select case (name)
    case ('GYRE_DIR')
       GYRE_DIR = value
    case default
       $ABORT(Invalid name)
    end select

    ! Finish

    return

  end subroutine set_constant_c_

end module constants_m
