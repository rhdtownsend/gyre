! Module   : gyre_math_crmath
! Purpose  : core math functions (crmath version)
!
! Copyright 2020 Rich Townsend & The GYRE Team
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

module gyre_math

  ! Uses

  use core_kinds

  use crmath

  ! No implicit typing

  implicit none

  ! Module variables

  real(DP), save, protected :: PI
  real(DP), save, protected :: TWOPI
  real(DP), save, protected :: HALFPI

  real(DP), save, protected :: DEG_TO_RAD
  real(DP), save, protected :: RAD_TO_DEG

  ! Interfaces

  interface pow
     module procedure pow_r_r_
     module procedure pow_r_c_
     module procedure pow_c_r_
     module procedure pow_c_c_
  end interface pow

  ! Access specifiers

  private

  public :: PI
  public :: HALFPI
  public :: TWOPI
  public :: DEG_TO_RAD
  public :: RAD_TO_DEG
  public :: init_math
  public :: log
  public :: log10
  public :: exp
  public :: sqrt
  public :: pow
  public :: abs
  public :: hypot
  public :: cos
  public :: sin
  public :: tan
  public :: cospi
  public :: sinpi
  public :: tanpi
  public :: acos
  public :: asin
  public :: atan
  public :: atan2
  public :: acospi
  public :: asinpi
  public :: atanpi
  public :: cosh
  public :: sinh
  public :: tanh
  public :: acosh
  public :: asinh
  public :: atanh

  ! Procedures

contains

  subroutine init_math ()

    call crmath_init()

    PI = acos(-1._DP)
    TWOPI = 2._DP*PI
    HALFPI = 0.5_DP*PI

    DEG_TO_RAD = PI/180._DP
    RAD_TO_DEG = 1._DP/DEG_TO_RAD

  end subroutine init_math

  !****

  elemental function pow_r_r_ (x, y) result (pow_xy)

    real(WP), intent(in) :: x
    real(WP), intent(in) :: y
    real(WP)             :: pow_xy

    if (x == 0._WP) then

       pow_xy = 0._WP

    else

       pow_xy = exp(log(x)*y)

    end if

  end function pow_r_r_

  !****
  
  elemental function pow_r_c_ (x, y) result (pow_xy)

    real(WP), intent(in)    :: x
    complex(WP), intent(in) :: y
    complex(WP)             :: pow_xy

    if (x == 0._WP) then

       pow_xy = 0._WP

    else

       pow_xy = exp(log(x)*y)

    end if

  end function pow_r_c_

  !****

  elemental function pow_c_r_ (x, y) result (pow_xy)

    complex(WP), intent(in) :: x
    real(WP), intent(in)    :: y
    complex(WP)             :: pow_xy

    if (x == 0._WP) then

       pow_xy = 0._WP

    else

       pow_xy = exp(log(x)*y)

    end if

  end function pow_c_r_

  !****

  elemental function pow_c_c_ (x, y) result (pow_xy)

    complex(WP), intent(in) :: x
    complex(WP), intent(in) :: y
    complex(WP)             :: pow_xy

    if (x == 0._WP) then

       pow_xy = 0._WP

    else

       pow_xy = exp(log(x)*y)

    end if

  end function pow_c_c_

end module gyre_math
