! Module   : gyre_math_intrinsic
! Purpose  : core math functions (crlibm version)
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

  use ISO_C_BINDING
  use IEEE_ARITHMETIC

  ! No implicit typing

  implicit none

  ! Parameters

  complex(DP), parameter :: I = CMPLX(0._DP, 1._DP, DP)

  ! Module variables

  real(DP), save, protected :: PI
  real(DP), save, protected :: TWOPI
  real(DP), save, protected :: HALFPI

  real(DP), save, protected :: DEG_TO_RAD
  real(DP), save, protected :: RAD_TO_DEG

  ! Interfaces

  ! C bindings

  interface

     pure subroutine crlibm_init () bind (C)
       use ISO_C_BINDING
     end subroutine crlibm_init

     pure function log_rz (x) result (log_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: log_x
     end function log_rz

     pure function log10_rz (x) result (log10_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: log10_x
     end function log10_rz

     pure function exp_rd (x) result (exp_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: exp_x
     end function exp_rd

     pure function cos_rz (x) result (cos_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: cos_x
     end function cos_rz

     pure function sin_rz (x) result (sin_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: sin_x
     end function sin_rz

     pure function tan_rz (x) result (tan_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: tan_x
     end function tan_rz

     pure function cospi_rz (x) result (cospi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: cospi_x
     end function cospi_rz

     pure function sinpi_rz (x) result (sinpi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: sinpi_x
     end function sinpi_rz

     pure function tanpi_rz (x) result (tanpi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: tanpi_x
     end function tanpi_rz

     pure function acos_rd (x) result (acos_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: acos_x
     end function acos_rd

     pure function asin_rz (x) result (asin_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: asin_x
     end function asin_rz

     pure function atan_rz (x) result (atan_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: atan_x
     end function atan_rz

     pure function acospi_rd (x) result (acospi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: acospi_x
     end function acospi_rd

     pure function asinpi_rz (x) result (asinpi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: asinpi_x
     end function asinpi_rz

     pure function atanpi_rz (x) result (atanpi_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: atanpi_x
     end function atanpi_rz

     pure function cosh_rz (x) result (cosh_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: cosh_x
     end function cosh_rz

     pure function sinh_rz (x) result (sinh_x) bind (C)
       use ISO_C_BINDING
       real(C_DOUBLE), value, intent(in) :: x
       real(C_DOUBLE)                    :: sinh_x
     end function sinh_rz

  end interface

  ! Generic interfaces

  interface log
     module procedure log_r_
     module procedure log_c_
  end interface log

  interface log10
     module procedure log10_r_
  end interface log10

  interface exp
     module procedure exp_r_
     module procedure exp_c_
  end interface exp

  interface pow
     module procedure pow_r_r_
     module procedure pow_r_c_
     module procedure pow_c_r_
     module procedure pow_c_c_
  end interface pow

  interface cos
     module procedure cos_r_
     module procedure cos_c_
  end interface cos

  interface sin
     module procedure sin_r_
     module procedure sin_c_
  end interface sin

  interface tan
     module procedure tan_r_
     module procedure tan_c_
  end interface tan

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

  interface acos
     module procedure acos_r_
     module procedure acos_c_
  end interface acos

  interface asin
     module procedure asin_r_
     module procedure asin_c_
  end interface asin

  interface atan
     module procedure atan_r_
     module procedure atan_c_
  end interface atan

  interface atan2
     module procedure atan2_r_
  end interface atan2

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

  interface cosh
     module procedure cosh_r_
     module procedure cosh_c_
  end interface cosh

  interface sinh
     module procedure sinh_r_
     module procedure sinh_c_
  end interface sinh

  interface tanh
     module procedure tanh_r_
     module procedure tanh_c_
  end interface tanh

  interface acosh
     module procedure acosh_r_
     module procedure acosh_c_
  end interface acosh

  interface asinh
     module procedure asinh_r_
     module procedure asinh_c_
  end interface asinh

  interface atanh
     module procedure atanh_r_
     module procedure atanh_c_
  end interface atanh

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
  public :: pow
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

    call crlibm_init()

    PI = acos(-1._DP)
    TWOPI = 2._DP*PI
    HALFPI = 0.5_DP*PI

    DEG_TO_RAD = PI/180._DP
    RAD_TO_DEG = 1._DP/DEG_TO_RAD

  end subroutine init_math

  !****

  elemental function log_r_ (x) result (log_x)

    real(DP), intent(in) :: x
    real(DP)             :: log_x

    log_x = log_rz(x)

  end function log_r_

  !****

  elemental function log_c_ (x) result (log_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: log_x

    log_x = CMPLX(log(ABS(x)), atan2(x%im, x%re), DP)

  end function log_c_

  !****

  elemental function log10_r_ (x) result (log10_x)

    real(DP), intent(in) :: x
    real(DP)             :: log10_x

    log10_x = log10_rz(x)

  end function log10_r_

  !****

  elemental function exp_r_ (x) result (exp_x)

    real(DP), intent(in) :: x
    real(DP)             :: exp_x

    exp_x = exp_rd(x)

  end function exp_r_

  !****

  elemental function exp_c_ (x) result (exp_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: exp_x

    exp_x = exp(x%re)*CMPLX(COS(x%im), SIN(x%im), DP)

  end function exp_c_

  !****

  elemental function pow_r_r_ (x, y) result (pow_xy)

    real(DP), intent(in) :: x
    real(DP), intent(in) :: y
    real(DP)             :: pow_xy

    if (x == 0._DP) then

       pow_xy = 0._DP

    else

       pow_xy = exp(log(x)*y)

    end if

  end function pow_r_r_

  !****

  elemental function pow_r_c_ (x, y) result (pow_xy)

    real(DP), intent(in)    :: x
    complex(DP), intent(in) :: y
    complex(DP)             :: pow_xy

    if (x == 0._DP) then

       pow_xy = 0._DP

    else

       pow_xy = exp(log(x)*y)

    end if

  end function pow_r_c_

  !****

  elemental function pow_c_r_ (x, y) result (pow_xy)

    complex(DP), intent(in) :: x
    real(DP), intent(in)    :: y
    complex(DP)             :: pow_xy

    if (x == 0._DP) then

       pow_xy = 0._DP

    else

       pow_xy = exp(log(x)*y)

    end if

  end function pow_c_r_

  !****

  elemental function pow_c_c_ (x, y) result (pow_xy)

    complex(DP), intent(in) :: x
    complex(DP), intent(in) :: y
    complex(DP)             :: pow_xy

    if (x == 0._DP) then

       pow_xy = 0._DP

    else

       pow_xy = exp(log(x)*y)

    end if

  end function pow_c_c_

  !****

  elemental function cos_r_ (x) result (cos_x)

    real(DP), intent(in) :: x
    real(DP)             :: cos_x

    cos_x = cos_rz(x)

  end function cos_r_

  !****

  elemental function cos_c_ (x) result (cos_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: cos_x

    cos_x = CMPLX(cos(x%re)*cosh(x%im), -sin(x%re)*sinh(x%im), DP)

  end function cos_c_

  !****

  elemental function sin_r_ (x) result (sin_x)

    real(DP), intent(in) :: x
    real(DP)             :: sin_x

    sin_x = sin_rz(x)

  end function sin_r_

  !****

  elemental function sin_c_ (x) result (sin_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: sin_x

    sin_x = CMPLX(sin(x%re)*cosh(x%im), cos(x%re)*sinh(x%im), DP)

  end function sin_c_
  
  !****

  elemental function tan_r_ (x) result (tan_x)

    real(DP), intent(in) :: x
    real(DP)             :: tan_x

    tan_x = tan_rz(x)

  end function tan_r_

  !****

  elemental function tan_c_ (x) result (tan_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: tan_x

    real(DP) :: rt
    real(DP) :: it

    rt = tan(x%re)
    it = tanh(x%im)

    tan_x = CMPLX(rt, it, DP)/CMPLX(1._DP, -rt*it, DP)

  end function tan_c_

  !****

  elemental function cospi_r_ (x) result (cospi_x)

    real(DP), intent(in) :: x
    real(DP)             :: cospi_x

    cospi_x = cospi_rz(x)

  end function cospi_r_

  !****

  elemental function cospi_c_ (x) result (cospi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: cospi_x

    cospi_x = CMPLX(cospi(x%re)*cosh(PI*x%im), -sinpi(x%re)*sinh(PI*x%im), DP)

  end function cospi_c_

  !****

  elemental function sinpi_r_ (x) result (sinpi_x)

    real(DP), intent(in) :: x
    real(DP)             :: sinpi_x

    sinpi_x = sinpi_rz(x)

  end function sinpi_r_

  !****

  elemental function sinpi_c_ (x) result (sinpi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: sinpi_x

    sinpi_x = CMPLX(sinpi(x%re)*cosh(PI*x%im), cospi(x%re)*sinh(PI*x%im), DP)

  end function sinpi_c_

  !****

  elemental function tanpi_r_ (x) result (tanpi_x)

    real(DP), intent(in) :: x
    real(DP)             :: tanpi_x

    tanpi_x = tanpi_rz(x)

  end function tanpi_r_

  !****

  elemental function tanpi_c_ (x) result (tanpi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: tanpi_x

    real(DP) :: rt
    real(DP) :: it

    rt = tanpi(x%re)
    it = tanh(PI*x%im)

    tanpi_x = CMPLX(rt, it, DP)/CMPLX(1._DP, -rt*it, DP)

  end function tanpi_c_

  !****

  elemental function acos_r_ (x) result (acos_x)

    real(DP), intent(in) :: x
    real(DP)             :: acos_x

    acos_x = acos_rd(x)

  end function acos_r_

  !****

  elemental function acos_c_ (x) result (acos_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: acos_x

    acos_x = -I*log(x + I*SQRT(1._dp - x**2))
    
  end function acos_c_

  !****

  elemental function asin_r_ (x) result (asin_x)

    real(DP), intent(in) :: x
    real(DP)             :: asin_x

    asin_x = asin_rz(x)

  end function asin_r_

  !****

  elemental function asin_c_ (x) result (asin_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: asin_x

    asin_x = -I*log(I*x + SQRT(1._DP - x**2))

  end function asin_c_

  !****

  elemental function atan_r_ (x) result (atan_x)

    real(DP), intent(in) :: x
    real(DP)             :: atan_x

    atan_x = atan_rz(x)

  end function atan_r_

  !****

  elemental function atan_c_ (x) result (atan_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: atan_x

    atan_x = I*log((I + x)/(I - x))/2._DP

  end function atan_c_

  !****

  elemental function atan2_r_ (y, x) result (atan2_yx)

    real(DP), intent(in) :: y
    real(DP), intent(in) :: x
    real(DP)             :: atan2_yx

    if (x > 0._DP) then
       atan2_yx = atan(y/x)
    elseif (y > 0._WP) then
       atan2_yx = HALFPI - atan(x/y)
    elseif (y < 0._WP) then
       atan2_yx = -HALFPI - atan(x/y)
    elseif (x < 0._WP) then
       atan2_yx = atan(y/x) + PI
    endif

  end function atan2_r_

  !****

  elemental function acospi_r_ (x) result (acospi_x)

    real(DP), intent(in) :: x
    real(DP)             :: acospi_x

    acospi_x = acospi_rd(x)

  end function acospi_r_

  !****

  elemental function acospi_c_ (x) result (acospi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: acospi_x

    acospi_x = acos(x)/PI

  end function acospi_c_

  !****

  elemental function asinpi_r_ (x) result (asinpi_x)

    real(DP), intent(in) :: x
    real(DP)             :: asinpi_x

    asinpi_x = asinpi_rz(x)

  end function asinpi_r_

  !****
  
  elemental function asinpi_c_ (x) result (asinpi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: asinpi_x

    asinpi_x = asin(x)/PI

  end function asinpi_c_
  
  !****

  elemental function atanpi_r_ (x) result (atanpi_x)

    real(DP), intent(in) :: x
    real(DP)             :: atanpi_x

    atanpi_x = atanpi_rz(x)

  end function atanpi_r_

  !****
  
  elemental function atanpi_c_ (x) result (atanpi_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: atanpi_x

    atanpi_x = atan(x)/PI

  end function atanpi_c_
  
  !****

  elemental function cosh_r_ (x) result (cosh_x)

    real(DP), intent(in) :: x
    real(DP)             :: cosh_x

    cosh_x = cosh_rz(x)

  end function cosh_r_

  !****

  elemental function cosh_c_ (x) result (cosh_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: cosh_x

    cosh_x = CMPLX(cosh(x%re)*cos(x%im), sinh(x%re)*sin(x%im), DP)

  end function cosh_c_

  !****

  elemental function sinh_r_ (x) result (sinh_x)

    real(DP), intent(in) :: x
    real(DP)             :: sinh_x

    sinh_x = sinh_rz(x)

  end function sinh_r_

  !****

  elemental function sinh_c_ (x) result (sinh_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: sinh_x

    sinh_x = CMPLX(sinh(x%re)*cos(x%im), cosh(x%re)*sin(x%im), DP)

  end function sinh_c_

  !****

  elemental function tanh_r_ (x) result (tanh_x)

    real(DP), intent(in) :: x
    real(DP)             :: tanh_x

    tanh_x = sinh(x)/cosh(x)

  end function tanh_r_

  !****

  elemental function tanh_c_ (x) result (tanh_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: tanh_x

    real(DP) :: rt
    real(DP) :: it

    rt = tanh(x%re)
    it = tan(x%im)

    tanh_x = CMPLX(rt, it, DP)/CMPLX(1._DP, rt*it, DP)

  end function tanh_c_

  !****

  elemental function acosh_r_ (x) result (acosh_x)

    real(DP), intent(in) :: x
    real(DP)             :: acosh_x

    acosh_x = log(x + SQRT((x - 1._DP)*(x + 1._DP)))

  end function acosh_r_

  !****

  elemental function acosh_c_ (x) result (acosh_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: acosh_x

    acosh_x = log(x + SQRT((x - 1._DP)*(x + 1._DP)))

  end function acosh_c_

  !****

  elemental function asinh_r_ (x) result (asinh_x)

    real(DP), intent(in) :: x
    real(DP)             :: asinh_x

    asinh_x = log(x + SQRT(1._DP + x**2))

  end function asinh_r_

  !****

  elemental function asinh_c_ (x) result (asinh_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: asinh_x

    asinh_x = log(x + SQRT(1._DP + x**2))

  end function asinh_c_

  !****

  elemental function atanh_r_ (x) result (atanh_x)

    real(DP), intent(in) :: x
    real(DP)             :: atanh_x

    atanh_x = log((1._DP + x)/(1._DP - x))/2._DP

  end function atanh_r_

  !****

  elemental function atanh_c_ (x) result (atanh_x)

    complex(DP), intent(in) :: x
    complex(DP)             :: atanh_x

    real(DP) :: rt
    real(DP) :: it

    atanh_x = log((1._DP + x)/(1._DP - x))/2._DP

  end function atanh_c_

end module gyre_math
