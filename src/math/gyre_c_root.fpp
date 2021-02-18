! Module   : gyre_c_root
! Purpose  : root finding algorithms (complex & c_ext_t)
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

module gyre_c_root

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

  interface solve_root
     module procedure solve_root_c_
     module procedure solve_root_cx_
  end interface solve_root

  interface narrow_bracket
     module procedure narrow_bracket_c_
     module procedure narrow_bracket_cx_
  end interface narrow_bracket

  interface expand_bracket
     module procedure expand_bracket_c_ 
     module procedure expand_bracket_cx_ 
  end interface expand_bracket

  ! Access specifiers

  private

  public :: solve_root
  public :: narrow_bracket
  public :: expand_bracket

contains

  $define $SOLVE $sub

  $local $T $1
  $local $TYPE_R $2
  $local $TYPE_C $3

  subroutine solve_root_${T}_ (eval_func, z_a, z_b, z_tol, nm_p, z_root, status, n_iter, n_iter_max, relative_tol, f_z_a, f_z_b)

    interface
       subroutine eval_func (z, func, status)
         use core_kinds
         use gyre_ext
         $TYPE_C, intent(in)  :: z
         $TYPE_C, intent(out) :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE_C, intent(in)            :: z_a
    $TYPE_C, intent(in)            :: z_b
    $TYPE_R, intent(in)            :: z_tol
    class(num_par_t), intent(in)   :: nm_p
    $TYPE_C, intent(out)           :: z_root
    integer, intent(out)           :: status
    integer, optional, intent(out) :: n_iter
    integer, optional, intent(in)  :: n_iter_max
    logical, optional, intent(in)  :: relative_tol
    $TYPE_C, optional, intent(in)  :: f_z_a
    $TYPE_C, optional, intent(in)  :: f_z_b

    $TYPE_C :: a
    $TYPE_C :: b
    $TYPE_C :: f_a
    $TYPE_C :: f_b

    ! Starting from the bracket [z_a,z_b], find a root of the function

    a = z_a
    b = z_b

    if (PRESENT(f_z_a)) then
       f_a = f_z_a
    else
       call eval_func(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_z_b)) then
       f_b = f_z_b
    else
       call eval_func(b, f_b, status)
       if (status /= STATUS_OK) return
    endif

    call narrow_bracket_${T}_(eval_func, a, b, z_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_a, f_b)

    z_root = b

    ! Finish

    return

  end subroutine solve_root_${T}_

  $endsub

  $SOLVE(c,real(WP),complex(WP))
  $SOLVE(cx,type(r_ext_t),type(c_ext_t))

  !****

  $define $NARROW_BRACKET $sub

  $local $T $1
  $local $TYPE_R $2
  $local $TYPE_C $3

  subroutine narrow_bracket_${T}_ (eval_func, z_a, z_b, z_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_z_a, f_z_b)

    interface
       subroutine eval_func (z, func, status)
         use core_kinds
         use gyre_ext
         $TYPE_C, intent(in)  :: z
         $TYPE_C, intent(out) :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE_C, intent(inout)           :: z_a
    $TYPE_C, intent(inout)           :: z_b
    $TYPE_R, intent(in)              :: z_tol
    class(num_par_t), intent(in)     :: nm_p
    integer, intent(out)             :: status
    integer, optional, intent(out)   :: n_iter
    integer, optional, intent(in)    :: n_iter_max
    logical, optional, intent(in)    :: relative_tol
    $TYPE_C, optional, intent(inout) :: f_z_a
    $TYPE_C, optional, intent(inout) :: f_z_b

    ! Narrow the bracket [z_a,z_b] toward a root of the function

    select case (nm_p%c_root_solver)
    case ('SECANT')
       call narrow_secant_${T}_(eval_func, z_a, z_b, z_tol, status, n_iter, n_iter_max, relative_tol, f_z_a, f_z_b)
    case ('RIDDERS')
       call narrow_ridders_${T}_(eval_func, z_a, z_b, z_tol, status, n_iter, n_iter_max, relative_tol, f_z_a, f_z_b)
    case default
       $ABORT(Invalid c_root_solver)
    end select

    ! Finish

    return

  end subroutine narrow_bracket_${T}_

  $endsub

  $NARROW_BRACKET(c,real(WP),complex(WP))
  $NARROW_BRACKET(cx,type(r_ext_t),type(c_ext_t))

  !****

  $define $NARROW_SECANT $sub

  $local $T $1
  $local $TYPE_R $2
  $local $TYPE_C $3

  subroutine narrow_secant_${T}_ (eval_func, z_a, z_b, z_tol, status, n_iter, n_iter_max, relative_tol, f_z_a, f_z_b)
 
    interface
       subroutine eval_func (z, func, status)
         use core_kinds
         use gyre_ext
         $TYPE_C, intent(in)  :: z
         $TYPE_C, intent(out) :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE_C, intent(inout)           :: z_a
    $TYPE_C, intent(inout)           :: z_b
    $TYPE_R, intent(in)              :: z_tol
    integer, intent(out)             :: status
    integer, optional, intent(out)   :: n_iter
    integer, optional, intent(in)    :: n_iter_max
    logical, optional, intent(in)    :: relative_tol
    $TYPE_C, optional, intent(inout) :: f_z_a
    $TYPE_C, optional, intent(inout) :: f_z_b

    logical :: relative_tol_
    $TYPE_C :: a
    $TYPE_C :: b
    $TYPE_C :: c
    $TYPE_C :: f_a
    $TYPE_C :: f_b
    $TYPE_C :: f_c
    integer :: i_iter
    $TYPE_C :: f_dz
    $TYPE_C :: rho
    $TYPE_R :: tol

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Narrow the bracket [z_a,z_b] toward a root of the function !
    ! using the secant method

    ! Set up the initial state

    a = z_a
    b = z_b

    if (PRESENT(f_z_a)) then
       f_a = f_z_a
    else
       call eval_func(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_z_b)) then
       f_b = f_z_b
    else
       call eval_func(b, f_b, status)
       if (status /= STATUS_OK) return
    endif

    if (abs(f_a) < abs(f_b)) then

       c = a
       a = b
       b = c

       f_c = f_a
       f_a = f_b
       f_b = f_c

    endif

    ! Iterate until convergence to the desired tolerance, or the
    ! maximum number of iterations is exceeded

    i_iter = 0

    status = STATUS_OK

    iterate_loop : do

       if (f_b == 0._WP) exit iterate_loop

       i_iter = i_iter + 1

       if (PRESENT(n_iter_max)) then
          if (i_iter > n_iter_max) then
             status = STATUS_ITER_MAX
             exit iterate_loop
          endif
       endif

       ! Calculate the correction

       f_dz = f_b*(b - a)

       rho = f_b - f_a

       ! Check for a singular correction

       if (abs(b*rho) < 8._WP*EPSILON(0._WP)*abs(f_dz)) then
          $ABORT(Singular correction in secant)
       endif

       ! Update the root

       a = b
       f_a = f_b

       b = b - f_dz/rho
       call eval_func(b, f_b, status)
       if (status /= STATUS_OK) exit iterate_loop

       ! Check for convergence

       if (relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + z_tol)*abs(b)
       else
          tol = 4._WP*EPSILON(0._WP)*abs(b) + z_tol
       endif

       if (abs(b - a) <= tol) exit iterate_loop

    end do iterate_loop

    ! Store the results

    z_a = a
    z_b = b

    if (PRESENT(n_iter)) n_iter = i_iter

    if (PRESENT(f_z_a)) f_z_a = f_a
    if (PRESENT(f_z_b)) f_z_b = f_b

    ! Finish

  end subroutine narrow_secant_${T}_

  $endsub

  $NARROW_SECANT(c,real(WP),complex(WP))
  $NARROW_SECANT(cx,type(r_ext_t),type(c_ext_t))

  !****

  $define $NARROW_RIDDERS $sub

  $local $T $1
  $local $TYPE_R $2
  $local $TYPE_C $3

  subroutine narrow_ridders_${T}_ (eval_func, z_a, z_b, z_tol, status, n_iter, n_iter_max, relative_tol, f_z_a, f_z_b)

    interface
       subroutine eval_func (z, func, status)
         use core_kinds
         use gyre_ext
         $TYPE_C, intent(in)  :: z
         $TYPE_C, intent(out) :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE_C, intent(inout)           :: z_a
    $TYPE_C, intent(inout)           :: z_b
    $TYPE_R, intent(in)              :: z_tol
    integer, intent(out)             :: status
    integer, optional, intent(out)   :: n_iter
    integer, optional, intent(in)    :: n_iter_max
    logical, optional, intent(in)    :: relative_tol
    $TYPE_C, optional, intent(inout) :: f_z_a
    $TYPE_C, optional, intent(inout) :: f_z_b

    logical :: relative_tol_
    $TYPE_C :: a
    $TYPE_C :: b
    $TYPE_C :: c
    $TYPE_C :: f_a
    $TYPE_C :: f_b
    $TYPE_C :: f_c
    integer :: i_iter
    $TYPE_C :: exp_Q_p
    $TYPE_C :: exp_Q_m
    $TYPE_C :: exp_Q
    $TYPE_C :: f_dz
    $TYPE_C :: rho
    $TYPE_R :: tol

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    $ASSERT(z_a /= z_b,Invalid initial bracket)

    ! Narrow the bracket [z_a,z_b] toward a root of the function using
    ! a complex Ridders' method (with secant updates, rather than
    ! regula falsi)

    ! Set up the initial state

    a = z_a
    b = z_b

    if (PRESENT(f_z_a)) then
       f_a = f_z_a
    else
       call eval_func(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_z_b)) then
       f_b = f_z_b
    else
       call eval_func(b, f_b, status)
       if (status /= STATUS_OK) return
    endif

    if (abs(f_a) < abs(f_b)) then

       c = a
       a = b
       b = c

       f_c = f_a
       f_a = f_b
       f_b = f_c

    endif

    ! Iterate until convergence to the desired tolerance, or the
    ! maximum number of iterations is exceeded

    i_iter = 0

    status = STATUS_OK

    iterate_loop : do

       if (f_b == 0._WP) exit iterate_loop

       i_iter = i_iter + 1

       if (PRESENT(n_iter_max)) then
          if (i_iter > n_iter_max) then
             status = STATUS_ITER_MAX
             exit iterate_loop
          end if
       endif

       ! Calculate the mid-point values

       c =  0.5_WP*(a + b)

       call eval_func(c, f_c, status)
       if (status /= STATUS_OK) exit iterate_loop

       ! Solve for the re-scaling exponential

       exp_Q_p = (f_c + sqrt(f_c*f_c - f_a*f_b))/f_b
       exp_Q_m = (f_c - sqrt(f_c*f_c - f_a*f_b))/f_b

       if (abs(exp_Q_p-1._WP) < abs(exp_Q_m-1._WP)) then
          exp_Q = exp_Q_p
       else
          exp_Q = exp_Q_m
       endif

       ! Apply the secant method to the re-scaled problem
 
       f_dz = f_b*(exp_Q*exp_Q)*(b - a)

       rho = f_b*(exp_Q*exp_Q) - f_a

       ! Check for a singular correction

       if (abs(b*rho) < 8._WP*EPSILON(0._WP)*abs(f_dz)) then
          $ABORT(Singular correction in secant)
       endif

       ! Update the root

       a = b
       f_a = f_b

       b = b - f_dz/rho
       call eval_func(b, f_b, status)
       if (status /= STATUS_OK) exit iterate_loop

       ! Check for convergence

       if (relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + z_tol)*abs(b)
       else
          tol = 4._WP*EPSILON(0._WP)*abs(b) + z_tol
       endif

       if (abs(b - a) <= tol) exit iterate_loop

    end do iterate_loop

    ! Store the results

    z_a = a
    z_b = b

    if (PRESENT(n_iter)) n_iter = i_iter

    if (PRESENT(f_z_a)) f_z_a = f_a
    if (PRESENT(f_z_b)) f_z_b = f_b

    ! Finish

    return

  end subroutine narrow_ridders_${T}_

  $endsub

  $NARROW_RIDDERS(c,real(WP),complex(WP))
  $NARROW_RIDDERS(cx,type(r_ext_t),type(c_ext_t))

  !****

  $define $EXPAND_BRACKET $sub

  $local $T $1
  $local $TYPE_R $2
  $local $TYPE_C $3

  subroutine expand_bracket_${T}_ (eval_func, z_a, z_b, f_z_tol, status, clamp_a, clamp_b, relative_tol, f_z_a, f_z_b)

    interface
       subroutine eval_func (z, func, status)
         use core_kinds
         use gyre_ext
         $TYPE_C, intent(in)  :: z
         $TYPE_C, intent(out) :: func
         integer, intent(out) :: status
       end subroutine eval_func
    end interface
    $TYPE_C, intent(inout)         :: z_a
    $TYPE_C, intent(inout)         :: z_b
    $TYPE_R, intent(in)            :: f_z_tol
    integer, intent(out)           :: status
    $TYPE_C, optional, intent(out) :: f_z_a
    $TYPE_C, optional, intent(out) :: f_z_b
    logical, optional, intent(in)  :: clamp_a
    logical, optional, intent(in)  :: clamp_b
    logical, optional, intent(in)  :: relative_tol

    real(WP), parameter :: EXPAND_FACTOR = 1.6_WP

    logical :: relative_tol_
    logical :: clamp_a_
    logical :: clamp_b_
    $TYPE_C :: f_a
    $TYPE_C :: f_b
    $TYPE_R :: tol
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

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    $ASSERT(z_a /= z_b,Invalid initial bracket)

    ! Expand the bracket [z_a,z_b] until the difference between f(z_a)
    ! and f(z_b) exceeds the tolerance

    call eval_func(z_a, f_a, status)
    if (status /= STATUS_OK) return

    call eval_func(z_b, f_b, status)
    if (status /= STATUS_OK) return

    status = STATUS_OK

    expand_loop : do

       if (relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + f_z_tol)*max(abs(f_a), abs(f_b))
       else
          tol = 4._WP*EPSILON(0._WP)*max(abs(f_a), abs(f_b)) + f_z_tol
       endif

       if (abs(f_a - f_b) > tol) exit expand_loop

       if (clamp_a_) then
          move_a = .FALSE.
       elseif (clamp_b_) then
          move_a = .TRUE.
       else
          move_a = abs(f_b) > abs(f_a)
       endif

       if (move_a) then
          z_a = z_a + EXPAND_FACTOR*(z_a - z_b)
          call eval_func(z_a, f_a, status)
          if (status /= STATUS_OK) exit expand_loop
       else
          z_b = z_b + EXPAND_FACTOR*(z_b - z_a)
          call eval_func(z_b, f_b, status)
          if (status /= STATUS_OK) exit expand_loop
       endif

    end do expand_loop

    ! Store f_a and f_b

    if (PRESENT(f_z_a)) f_z_a = f_a
    if (PRESENT(f_z_b)) f_z_b = f_b

    ! Finish

    return

  end subroutine expand_bracket_${T}_

  $endsub

  $EXPAND_BRACKET(c,real(WP),complex(WP))
  $EXPAND_BRACKET(cx,type(r_ext_t),type(c_ext_t))

end module gyre_c_root
