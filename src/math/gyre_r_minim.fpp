! Module   : gyre_r_minim
! Purpose  : minimum finding algorithms (real & r_ext_t)
!
! Copyright 2018-2021 Rich Townsend & The GYRE Team
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

module gyre_r_minim

  ! Uses

  use core_kinds

  use gyre_ext
  use gyre_math
  use gyre_num_par
  use gyre_status

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface solve_minim
     module procedure solve_minim_r_
     module procedure solve_minim_rx_
  end interface solve_minim

  ! Access specifiers

  private

  public :: solve_minim

contains

  $define $SOLVE_MINIM $sub

  $local $T $1
  $local $TYPE $2

  subroutine solve_minim_${T}_ (eval_func, x_a, x_b, x_c, x_tol, nm_p, x_minim, status, n_iter, n_iter_max, relative_tol, f_x_a, f_x_b, f_x_c)

    interface
       subroutine eval_func (x, func, status)
         use core_kinds
         use gyre_ext
         $TYPE, intent(in)    :: x
         $TYPE, intent(out)   :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE, intent(in)             :: x_a
    $TYPE, intent(in)             :: x_b
    $TYPE, intent(in)             :: x_c
    $TYPE, intent(in)             :: x_tol
    class(num_par_t), intent(in)  :: nm_p
    $TYPE, intent(out)            :: x_minim
    integer, intent(out)          :: status
    integer, optional, intent(out):: n_iter
    integer, optional, intent(in) :: n_iter_max
    logical, optional, intent(in) :: relative_tol
    $TYPE, optional, intent(in)   :: f_x_a
    $TYPE, optional, intent(in)   :: f_x_b
    $TYPE, optional, intent(in)   :: f_x_c

    $TYPE :: a
    $TYPE :: b
    $TYPE :: c
    $TYPE :: f_a
    $TYPE :: f_b
    $TYPE :: f_c

    ! Starting from the bracket [x_a,x_b,x_c], find a minimum of
    ! the function

    a = x_a
    b = x_b
    c = x_c

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

    if (PRESENT(f_x_c)) then
       f_c = f_x_c
    else
       call eval_func(c, f_c, status)
       if (status /= STATUS_OK) return
    endif

    call narrow_bracket_${T}_(eval_func, a, b, c, x_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_a, f_b, f_c)

    x_minim = b

    ! Finish

    return

  end subroutine solve_minim_${T}_

  $endsub

  $SOLVE_MINIM(r,real(WP))
  $SOLVE_MINIM(rx,type(r_ext_t))

  !****

  $define $NARROW_BRACKET $sub

  $local $T $1
  $local $TYPE $2

  subroutine narrow_bracket_${T}_ (eval_func, x_a, x_b, x_c, x_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_x_a, f_x_b, f_x_c)

    interface
       subroutine eval_func (x, func, status)
         use core_kinds
         use gyre_ext
         $TYPE, intent(in)    :: x
         $TYPE, intent(out)   :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE, intent(inout)           :: x_a
    $TYPE, intent(inout)           :: x_b
    $TYPE, intent(inout)           :: x_c
    $TYPE, intent(in)              :: x_tol
    class(num_par_t), intent(in)   :: nm_p
    integer, intent(out)           :: status
    integer, optional, intent(out) :: n_iter
    integer, optional, intent(in)  :: n_iter_max
    logical, optional, intent(in)  :: relative_tol
    $TYPE, optional, intent(inout) :: f_x_a
    $TYPE, optional, intent(inout) :: f_x_b
    $TYPE, optional, intent(inout) :: f_x_c

    ! Narrow_Bracket the bracket [x_a,x_b,x_c] on a minimum of the function

    call narrow_bisect_${T}_(eval_func, x_a, x_b, x_c, x_tol, status, n_iter, n_iter_max, relative_tol, f_x_a, f_x_b, f_x_c)

    ! Finish

    return

  end subroutine narrow_bracket_${T}_

  $endsub

  $NARROW_BRACKET(r,real(WP))
  $NARROW_BRACKET(rx,type(r_ext_t))

  !****

  $define $NARROW_BISECT $sub

  $local $T $1
  $local $TYPE $2

  subroutine narrow_bisect_${T}_ (eval_func, x_a, x_b, x_c, x_tol, status, n_iter, n_iter_max, relative_tol, f_x_a, f_x_b, f_x_c)

    interface
       subroutine eval_func (x, func, status)
         use core_kinds
         use gyre_ext
         $TYPE, intent(in)    :: x
         $TYPE, intent(out)   :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE, intent(inout)           :: x_a
    $TYPE, intent(inout)           :: x_b
    $TYPE, intent(inout)           :: x_c
    $TYPE, intent(in)              :: x_tol
    integer, intent(out)           :: status
    integer, optional, intent(out) :: n_iter
    integer, optional, intent(in)  :: n_iter_max
    logical, optional, intent(in)  :: relative_tol
    $TYPE, optional, intent(inout) :: f_x_a
    $TYPE, optional, intent(inout) :: f_x_b
    $TYPE, optional, intent(inout) :: f_x_c

    real(WP), parameter :: R = 0.61803399_WP

    logical :: relative_tol_
    $TYPE   :: a
    $TYPE   :: b
    $TYPE   :: c
    $TYPE   :: d
    $TYPE   :: f_a
    $TYPE   :: f_b
    $TYPE   :: f_c
    $TYPE   :: f_d
    integer :: i_iter
    $TYPE   :: tol

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Narrow the bracket [x_a,x_b,x_c] on a minimum of the function
    ! using a bisection method

    ! Set up the initial state

    a = x_a
    b = x_b
    c = x_c

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

    if (PRESENT(f_x_c)) then
       f_c = f_x_c
    else
       call eval_func(c, f_c, status)
       if (status /= STATUS_OK) return
    endif

    ! Confirm that the bracket is properly ordered and contains a
    ! minimum

    $ASSERT(b > a,Invalid bracket)
    $ASSERT(c > b,Invalid bracket)

    $ASSERT((f_a > f_b) .AND. (f_c > f_b),Minimum is not bracketed)

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

       ! Check for convegence

       if (relative_tol_) then
          tol = (sqrt(EPSILON(0._WP)) + 0.5_WP*x_tol)*(abs(a) + abs(c))
       else
          tol = sqrt(EPSILON(0._WP))*(abs(a) + abs(c)) + x_tol
       endif

       if (abs(c-a) <= tol) exit iterate_loop

       ! Create a new point, and update the bracket accordingly

       if (b-a > c-b) then

          ! New point in left interval
          
          d = (1._WP-R)*b + R*a

          call eval_func(d, f_d, status)
          if (status /= STATUS_OK) return

          if (f_d < f_b) then
             c = b
             b = d
             f_c = f_b
             f_b = f_d
          else
             a = d
             f_a = f_d
          endif

       else

          ! New point in right interval

          d = (1._WP-R)*b + R*c

          call eval_func(d, f_d, status)
          if (status /= STATUS_OK) return

          if (f_d < f_b) then
             a = b
             b = d
             f_a = f_b
             f_b = f_d
          else
             c = d
             f_c = f_d
          endif

       end if

    end do iterate_loop
       
    ! Store the results

    x_a = a
    x_b = b
    x_c = c

    if (PRESENT(n_iter)) then
       n_iter = i_iter
    end if

    if (PRESENT(f_x_a)) f_x_a = f_a
    if (PRESENT(f_x_b)) f_x_b = f_b
    if (PRESENT(f_x_c)) f_x_c = f_c

    ! Finish

    return

  end subroutine narrow_bisect_${T}_

  $endsub

  $NARROW_BISECT(r,real(WP))
  $NARROW_BISECT(rx,type(r_ext_t))

end module gyre_r_minim
