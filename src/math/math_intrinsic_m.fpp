! Module  : math_intrinsic_m
! Purpose : core math functions (intrinsic version)
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

module math_m

  ! Uses

  use kinds_m

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

  interface cospi
     module procedure cospi_r_
     module procedure cospi_c_
  end interface cospi

  interface sinpi
     module procedure sinpi_r_
     module procedure sinpi_c_
  end interface sinpi

  interface tanpi
     module procedure tanpi_r_
     module procedure tanpi_c_
  end interface tanpi

  interface acospi
     module procedure acospi_r_
     module procedure acospi_c_
  end interface acospi

  interface asinpi
     module procedure asinpi_r_
     module procedure asinpi_c_
  end interface asinpi

  interface atanpi
     module procedure atanpi_r_
     module procedure atanpi_c_
  end interface atanpi

  ! Access specifiers

  private

  public :: PI
  public :: HALFPI
  public :: TWOPI
  public :: DEG_TO_RAD
  public :: RAD_TO_DEG
  public :: init_math
  public :: pow
  public :: cospi
  public :: sinpi
  public :: tanpi
  public :: acospi
  public :: asinpi
  public :: atanpi

  ! Procedures

contains

  subroutine init_math ()

    PI = acos(-1._DP)
    TWOPI = 2._DP*PI
    HALFPI = 0.5_DP*PI

    DEG_TO_RAD = PI/180._DP
    RAD_TO_DEG = 1._DP/DEG_TO_RAD

  end subroutine init_math

  !****

  elemental function pow_r_r_ (x, y) result (pow_x)

    real(dp), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp)             :: pow_x

    pow_x = x**y
    
  end function pow_r_r_

  !****

  elemental function pow_r_c_ (x, y) result (pow_x)

    real(dp), intent(in)    :: x
    complex(dp), intent(in) :: y
    complex(dp)             :: pow_x

    pow_x = x**y
    
  end function pow_r_c_

  !****

  elemental function pow_c_r_ (x, y) result (pow_x)

    complex(DP), intent(in) :: x
    real(DP), intent(in)    :: y
    complex(DP)             :: pow_x

    pow_x = x**y
    
  end function pow_c_r_

  !****

  elemental function pow_c_c_ (x, y) result (pow_x)

    complex(DP), intent(in) :: x
    complex(DP), intent(in) :: y
    complex(DP)             :: pow_x

    pow_x = x**y
    
  end function pow_c_c_

  !****

  elemental function cospi_r_ (x) result (cospi_x)

    real(DP), intent(in) :: x
    real(DP)             :: cospi_x

    cospi_x = COS(x*PI)

  end function cospi_r_
  
  !****

  elemental function cospi_c_ (x) result (cospi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: cospi_x

    cospi_x = COS(x*PI)

  end function cospi_c_
  
  !****

  elemental function sinpi_r_ (x) result (sinpi_x)

    real(DP), intent(in) :: x
    real(DP)             :: sinpi_x

    sinpi_x = SIN(x*PI)

  end function sinpi_r_

  !****

  elemental function sinpi_c_ (x) result (sinpi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: sinpi_x

    sinpi_x = SIN(x*PI)

  end function sinpi_c_

  !****

  elemental function tanpi_r_ (x) result (tanpi_x)

    real(DP), intent(in) :: x
    real(DP)             :: tanpi_x

    tanpi_x = TAN(x*PI)

  end function tanpi_r_

  !****

  elemental function tanpi_c_ (x) result (tanpi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: tanpi_x

    tanpi_x = TAN(x*PI)

  end function tanpi_c_

  !****

  elemental function acospi_r_ (x) result (acospi_x)

    real(DP), intent(in) :: x
    real(DP)             :: acospi_x

    acospi_x = ACOS(x)/PI

  end function acospi_r_

  !****

  elemental function acospi_c_ (x) result (acospi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: acospi_x

    acospi_x = ACOS(x)/PI

  end function acospi_c_
  
  !****

  elemental function asinpi_r_ (x) result (asinpi_x)

    real(DP), intent(in) :: x
    real(DP)             :: asinpi_x

    asinpi_x = ASIN(x)/PI

  end function asinpi_r_

  !****

  elemental function asinpi_c_ (x) result (asinpi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: asinpi_x

    asinpi_x = ASIN(x)/PI

  end function asinpi_c_

  !****

  elemental function atanpi_r_ (x) result (atanpi_x)

    real(DP), intent(in) :: x
    real(DP)             :: atanpi_x

    atanpi_x = ATAN(x)/PI

  end function atanpi_r_

  !****

  elemental function atanpi_c_ (x) result (atanpi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: atanpi_x

    atanpi_x = ATAN(x)/PI

  end function atanpi_c_

end module math_m
