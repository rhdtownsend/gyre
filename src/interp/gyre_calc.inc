! Module   : gyre_calc
! Purpose  : calculus-related functions (template)
!
! Copyright 2017 Rich Townsend
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
$include 'core_parallel.inc'

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_calc

  ! Uses

  use core_kinds
  use core_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface deriv
     module procedure deriv_
  end interface deriv

  ! Access specifiers

  private

  public :: deriv

  ! Procedures

contains

  function deriv_ (x, f, deriv_type, df_dx_a, df_dx_b) result (df_dx)

    real(WP), intent(in)            :: x(:)
    $TYPE(WP), intent(in)           :: f(:)
    character(*), intent(in)        :: deriv_type
    $TYPE(WP), optional, intent(in) :: df_dx_a
    $TYPE(WP), optional, intent(in) :: df_dx_b
    $TYPE(WP)                       :: df_dx(SIZE(x))

    $CHECK_BOUNDS(SIZE(f),SIZE(x))

    ! Evaluate the derivative df/dx

    select case (deriv_type)
    case ('SPLINE')
       df_dx = deriv_spline_(x, f, df_dx_a, df_dx_b)
    case('FINDIFF')
       df_dx = deriv_findiff_(x, f, df_dx_a, df_dx_b)
    $if ($T eq 'r')
    case('MONO')
       df_dx = deriv_mono_(x, f, df_dx_a, df_dx_b)
    $endif
    case default
       $ABORT(Invalid deriv_type)
    end select

    ! Finish

    return

  end function deriv_

  !****

  function deriv_spline_ (x, f, df_dx_a, df_dx_b) result (df_dx)

    real(WP), intent(in)            :: x(:)
    $TYPE(WP), intent(in)           :: f(:)
    $TYPE(WP), intent(in), optional :: df_dx_a
    $TYPE(WP), intent(in), optional :: df_dx_b
    $TYPE(WP)                       :: df_dx(SIZE(x))

    integer   :: n
    $TYPE(WP) :: h(SIZE(x)-1)
    $TYPE(WP) :: L(SIZE(x)-1)
    $TYPE(WP) :: D(SIZE(x))
    $TYPE(WP) :: U(SIZE(x)-1)
    $TYPE(WP) :: B(SIZE(x),1)
    integer   :: info

    $CHECK_BOUNDS(SIZE(f),SIZE(x))
    
    ! Calcualte the first derivatives for a cubic spline interp
    ! (ensuring the second derivatives are continuous)

    n = SIZE(x)

    h = x(2:) - x(:n-1)

    ! Set up the tridiagonal matrix and RHS

    ! Inner boundary

    D(1) = 1._WP
    U(1) = 0._WP

    if (PRESENT(df_dx_a)) then
       B(1,1) = df_dx_a
    else
       B(1,1) = (f(2) - f(1))/h(1)
    endif

    ! Internal points

    L(1:n-2) = 2._WP/h(1:n-2)
    D(2:n-1) = 4._WP/h(1:n-2) + 4._WP/h(2:n-1)
    U(2:n-1) = 2._WP/h(2:n-1)

    B(2:n-1,1) = -6._WP*f(1:n-2)/h(1:n-2)**2 + 6._WP*f(2:n-1)/h(1:n-2)**2 + &
                  6._WP*f(3:n  )/h(2:n-1)**2 - 6._WP*f(2:n-1)/h(2:n-1)**2

    ! Outer boundary

    L(n-1) = 0._WP
    D(n) = 1._WP

    if (PRESENT(df_dx_b)) then
       B(n,1) = df_dx_b
    else
       B(n,1) = (f(n) - f(n-1))/h(n-1)
    endif

    ! Solve the tridiagonal system

    call XGTSV(n, 1, L, D, U, B, SIZE(B, 1), info)
    $ASSERT(info == 0,Non-zero return from XTGSV)

    df_dx = B(:,1)

    ! Finish

    return

  end function deriv_spline_

  !****

  function deriv_findiff_ (x, f, df_dx_a, df_dx_b) result (df_dx)

    real(WP), intent(in)            :: x(:)
    $TYPE(WP), intent(in)           :: f(:)
    $TYPE(WP), intent(in), optional :: df_dx_a
    $TYPE(WP), intent(in), optional :: df_dx_b
    $TYPE(WP)                       :: df_dx(SIZE(x))

    integer   :: n
    real(WP)  :: h(SIZE(x)-1)
    $TYPE(WP) :: s(SIZE(x)-1)

    $CHECK_BOUNDS(SIZE(f),SIZE(x))

    ! Calculate the first derivatives via centered finite differences

    n = SIZE(x)

    h = x(2:) - x(:n-1)

    s = (f(2:) - f(:n-1))/h

    if (PRESENT(df_dx_a)) then
       df_dx(1) = df_dx_a
    else
       df_dx(1) = s(1)
    endif

    df_dx(2:n-1) = 0.5_WP*(s(1:n-2) + s(2:n-1))

    if (PRESENT(df_dx_b)) then
       df_dx(n) = df_dx_b
    else
       df_dx(n) = s(n-1)
    endif

    ! Finish

    return

  end function deriv_findiff_

  !****

  $if ($T eq 'r')

  function deriv_mono_ (x, f, df_dx_a, df_dx_b) result (df_dx)

    real(WP), intent(in)            :: x(:)
    $TYPE(WP), intent(in)           :: f(:)
    $TYPE(WP), intent(in), optional :: df_dx_a
    $TYPE(WP), intent(in), optional :: df_dx_b
    $TYPE(WP)                       :: df_dx(SIZE(x))

    integer   :: n
    real(WP)  :: h(SIZE(x)-1)
    $TYPE(WP) :: s(SIZE(x)-1)
    $TYPE(WP) :: p(SIZE(x))
    integer   :: i

    $CHECK_BOUNDS(SIZE(f),SIZE(x))

    ! Calculate the first derivatives using the Steffen (1990, A&A,
    ! 239, 443) monontonicity preserving algorithm

    n = SIZE(x)

    h = x(2:) - x(:n-1)

    s = (f(2:) - f(:n-1))/h

    ! Calculate parabolic gradients

    if (PRESENT(df_dx_a)) then
       p(1) = df_dx_a
    else
       p(1) = s(1)
    endif

    p(2:n-1) = (s(1:n-2)*h(2:n-1) + s(2:n-1)*h(1:n-2))/(h(1:n-2) + h(2:n-1))

    if (PRESENT(df_dx_b)) then
       p(n) = df_dx_b
    else
       p(n) = s(n-1)
    endif

    ! Calculate monotonic gradients

    df_dx(1) = p(1)

    do i = 2,n-1
       df_dx(i) = (SIGN(1._WP, s(i-1)) + SIGN(1._WP, s(i)))* &
                  MIN(ABS(s(i-1)), ABS(s(i)), 0.5*ABS(p(i)))
    end do

    df_dx(n) = p(n)

    ! Finish

    return

  end function deriv_mono_

  $endif

end module gyre_${T}_calc

$endsub

