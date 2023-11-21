! Module  : r_root_m
! Purpose : root finding algorithms (real)
!
! Copyright 2013-2021 Rich Townsend & The GYRE Team
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

module r_root_m

  ! Uses

  use kinds_m

  use ext_m
  use math_m
  use num_par_m
  use status_m

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface solve_root
     module procedure solve_root_r_
     module procedure solve_root_rx_
  end interface solve_root

  interface narrow_bracket
     module procedure narrow_bracket_r_
     module procedure narrow_bracket_rx_
  end interface narrow_bracket

  interface expand_bracket
     module procedure expand_bracket_r_
     module procedure expand_bracket_rx_
  end interface expand_bracket

  ! Access specifiers

  private

  public :: solve_root
  public :: narrow_bracket
  public :: expand_bracket

contains

  $define $SOLVE_ROOT $sub

  $local $T $1
  $local $TYPE $2

  subroutine solve_root_${T}_ (eval_func, x_a, x_b, x_tol, nm_p, x_root, status, n_iter, n_iter_max, relative_tol, f_x_a, f_x_b)

    interface
       subroutine eval_func (x, func, status)
         use kinds_m
         use ext_m
         $TYPE, intent(in)    :: x
         $TYPE, intent(out)   :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE, intent(in)              :: x_a
    $TYPE, intent(in)              :: x_b
    $TYPE, intent(in)              :: x_tol
    class(num_par_t), intent(in)   :: nm_p
    $TYPE, intent(out)             :: x_root
    integer, intent(out)           :: status
    integer, optional, intent(out) :: n_iter
    integer, optional, intent(in)  :: n_iter_max
    logical, optional, intent(in)  :: relative_tol
    $TYPE, optional, intent(in)    :: f_x_a
    $TYPE, optional, intent(in)    :: f_x_b

    $TYPE :: a
    $TYPE :: b
    $TYPE :: f_a
    $TYPE :: f_b

    ! Starting from the bracket [x_a,x_b], which must span a sign
    ! change in the function, solve for a root of the function

    a = x_a
    b = x_b

    if (PRESENT(f_x_a)) then
       f_a = f_x_a
    else
       call eval_func(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_x_b)) then
       f_b = f_x_b
    else
       call eval_func(b, f_b, status)
       if (status /= STATUS_OK) return
    endif

    call narrow_bracket_${T}_(eval_func, a, b, x_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_a, f_b)

    x_root = b

    ! Finish

    return

  end subroutine solve_root_${T}_

  $endsub

  $SOLVE_ROOT(r,real(WP))
  $SOLVE_ROOT(rx,type(r_ext_t))

  !****

  $define $NARROW_BRACKET $sub

  $local $T $1
  $local $TYPE $2

  subroutine narrow_bracket_${T}_ (eval_func, x_a, x_b, x_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_x_a, f_x_b)

    interface
       subroutine eval_func (x, func, status)
         use kinds_m
         use ext_m
         $TYPE, intent(in)    :: x
         $TYPE, intent(out)   :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE, intent(inout)           :: x_a
    $TYPE, intent(inout)           :: x_b
    $TYPE, intent(in)              :: x_tol
    class(num_par_t), intent(in)   :: nm_p
    integer, intent(out)           :: status
    integer, optional, intent(out) :: n_iter
    integer, optional, intent(in)  :: n_iter_max
    logical, optional, intent(in)  :: relative_tol
    $TYPE, optional, intent(inout) :: f_x_a
    $TYPE, optional, intent(inout) :: f_x_b

    ! Narrow the bracket [x_a,x_b] on a root of the function

    select case (nm_p%r_root_solver)
    case ('BRENT')
       call narrow_brent_${T}_(eval_func, x_a, x_b, x_tol, status, n_iter, n_iter_max, relative_tol, f_x_a, f_x_b)
    case default
       $ABORT(Invalid r_root_solver)
    end select

    ! Finish

    return

  end subroutine narrow_bracket_${T}_

  $endsub

  $NARROW_BRACKET(r,real(WP))
  $NARROW_BRACKET(rx,type(r_ext_t))

  !****

  $define $NARROW_BRENT $sub

  $local $T $1
  $local $TYPE $2

  subroutine narrow_brent_${T}_ (eval_func, x_a, x_b, x_tol, status, n_iter, n_iter_max, relative_tol, f_x_a, f_x_b)

    interface
       subroutine eval_func (x, func, status)
         use kinds_m
         use ext_m
         $TYPE, intent(in)    :: x
         $TYPE, intent(out)   :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE, intent(inout)           :: x_a
    $TYPE, intent(inout)           :: x_b
    $TYPE, intent(in)              :: x_tol
    integer, intent(out)           :: status
    integer, optional, intent(out) :: n_iter
    integer, optional, intent(in)  :: n_iter_max
    logical, optional, intent(in)  :: relative_tol
    $TYPE, optional, intent(inout) :: f_x_a
    $TYPE, optional, intent(inout) :: f_x_b

    logical :: relative_tol_
    $TYPE   :: a
    $TYPE   :: b
    $TYPE   :: c
    $TYPE   :: d
    $TYPE   :: e
    $TYPE   :: f_a
    $TYPE   :: f_b
    $TYPE   :: f_c
    $TYPE   :: tol
    $TYPE   :: m
    $TYPE   :: p
    $TYPE   :: q
    $TYPE   :: r
    $TYPE   :: s
    integer :: i_iter

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Narrow the bracket [x_a,x_b] on a root of the function
    ! using Brent's method [based on the ALGOL 60 routine 'zero'
    ! published in [Bre1973]

    ! Set up the initial state

    a = x_a
    b = x_b

    if (PRESENT(f_x_a)) then
       f_a = f_x_a
    else
       call eval_func(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_x_b)) then
       f_b = f_x_b
    else
       call eval_func(b, f_b, status)
       if (status /= STATUS_OK) return
    endif

    c = b
    f_c = f_b

    ! Confirm that the bracket contains a root

    $ASSERT((f_a >= 0._WP .AND. f_b <= 0._WP) .OR. (f_a <= 0._WP .AND. f_b >= 0._WP),Root is not bracketed)

    ! Iterate until convergence to the desired tolerance, or the
    ! maximum number of iterations is exceeded

    i_iter = 0

    status = STATUS_OK

    iterate_loop : do

       i_iter = i_iter + 1

       if (PRESENT(n_iter_max)) then
          if (i_iter > n_iter_max) then
             status = STATUS_ITER_MAX
             exit iterate_loop
          endif
       endif

       ! Reorder c so that it has the opposite sign to b

       if (f_b > 0._WP .EQV. f_c > 0._WP) then
          c = a
          f_c = f_a
          d = b - a
          e = d
       endif

       ! Make sure that the function is smallest in magnitude
       ! at b
          
       if (abs(f_c) < abs(f_b)) then
          a = b
          b = c
          c = a
          f_a = f_b
          f_b = f_c
          f_c = f_a
       endif

       if (relative_tol_) then
          tol = (2._WP*EPSILON(0._WP) + x_tol)*abs(b)
       else
          tol = 2._WP*EPSILON(0._WP)*abs(b) + x_tol
       endif

       m = 0.5_WP*(c - b)

       ! Check for convergence

       if (abs(m) <= tol .OR. f_b == 0._WP) then
          exit iterate_loop
       endif

       ! See if bisection is forced

       if (abs(e) <  tol .OR. abs(f_a) < abs(f_b)) then

          d = m
          e = d

       else

          s = f_b/f_a

          if (a == c) then

             ! Linear interpolation

             p = 2._WP*m*s
             q = 1._WP - s

          else

             ! Inverse quadratic interpolation

             q = f_a/f_c
             r = f_b/f_c

             p = s*(2._WP*m*q*(q - r) - (b - a)*(r - 1._WP))
             q = (q - 1._WP)*(r - 1._WP)*(s - 1._WP)

          endif

          if (p > 0._WP) then
             q = -q
          else
             p = -p
          endif

          s = e
          e = d

          if (2._WP*p < 3._WP*m*q - abs(tol*q) .AND. p < abs(0.5_WP*s*q)) then
             d = p/q
          else
             d = m
             e = d
          endif

       endif

       ! Store the old value of b in a

       a = b
       f_a = f_b

       ! Update b

       if (abs(d) > tol) then
          b = b + d
       else
          if(m > 0._WP) then
             b = b + tol
          else
             b = b - tol
          endif
       endif

       call eval_func(b, f_b, status)
       if (status /= STATUS_OK) exit iterate_loop
       
    end do iterate_loop

    ! Store the results

    x_a = a
    x_b = b

    if (PRESENT(n_iter)) then
       n_iter = i_iter
    end if

    if (PRESENT(f_x_a)) f_x_a = f_a
    if (PRESENT(f_x_b)) f_x_b = f_b

    ! Finish

    return

  end subroutine narrow_brent_${T}_

  $endsub

  $NARROW_BRENT(r,real(WP))
  $NARROW_BRENT(rx,type(r_ext_t))

  !****

  $define $EXPAND_BRACKET $sub

  $local $T $1
  $local $TYPE $2

  subroutine expand_bracket_${T}_ (eval_func, x_a, x_b, status, clamp_a, clamp_b, f_x_a, f_x_b)

    interface
       subroutine eval_func (x, func, status)
         use kinds_m
         use ext_m
         $TYPE, intent(in)    :: x
         $TYPE, intent(out)   :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE, intent(inout)          :: x_a
    $TYPE, intent(inout)          :: x_b
    integer, intent(out)          :: status
    logical, optional, intent(in) :: clamp_a
    logical, optional, intent(in) :: clamp_b
    $TYPE, optional, intent(out)  :: f_x_a
    $TYPE, optional, intent(out)  :: f_x_b

    real(WP), parameter :: EXPAND_FACTOR = 1.6_WP

    logical :: clamp_a_
    logical :: clamp_b_
    $TYPE   :: f_a
    $TYPE   :: f_b
    logical :: move_a

    if (PRESENT(clamp_a)) then
       clamp_a_ = clamp_a
    else
       clamp_a_ = .FALSE.
    endif

    if (PRESENT(clamp_b)) then
       clamp_b_ = clamp_b
    else
       clamp_b_ = .FALSE.
    endif

    $ASSERT(.NOT. (clamp_a_ .AND. clamp_b_),Cannot clamp both points)

    $ASSERT(x_a /= x_b,Invalid initial bracket)

    ! Expand the bracket [x_a,x_b] until it contains a root of the
    ! function

    call eval_func(x_a, f_a, status)
    if (status /= STATUS_OK) return

    call eval_func(x_b, f_b, status)
    if (status /= STATUS_OK) return

    status = STATUS_OK

    expand_loop : do

       if ((f_a > 0._WP .AND. f_b < 0._WP) .OR. &
           (f_a < 0._WP .AND. f_b > 0._WP)) exit expand_loop

       if (clamp_a_) then
          move_a = .FALSE.
       elseif (clamp_b_) then
          move_a = .TRUE.
       else
          move_a = abs(f_b) > abs(f_a)
       endif

       if (move_a) then

          x_a = x_a + EXPAND_FACTOR*(x_a - x_b)

          call eval_func(x_a, f_a, status)
          if (status /= STATUS_OK) exit expand_loop

       else

          x_b = x_b + EXPAND_FACTOR*(x_b - x_a)

          call eval_func(x_b, f_b, status)
          if (status /= STATUS_OK) exit expand_loop

       endif

    end do expand_loop

    ! Store the results

    if (PRESENT(f_x_a)) f_x_a = f_a
    if (PRESENT(f_x_b)) f_x_b = f_b

    ! Finish

    return

  end subroutine expand_bracket_${T}_

  $endsub

  $EXPAND_BRACKET(r,real(WP))
  $EXPAND_BRACKET(rx,type(r_ext_t))

end module r_root_m
