! Module   : gyre_c_extfunc
! Purpose  : root finding with extended-range arithmetic (complex)
!
! Copyright 2013-2014 Rich Townsend
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

module gyre_c_extfunc

  ! Uses

  use core_kinds

  use gyre_ext_arith

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: c_extfunc_t
   contains
     private
     procedure(eval_), deferred :: eval
     procedure                  :: expand => expand_
     procedure                  :: narrow => narrow_ridders_
     procedure                  :: root => root_
  end type c_extfunc_t

  ! Interfaces

  abstract interface
    function eval_ (this, cx) result (f_cx)
      use gyre_c_ext
      import c_extfunc_t
      class(c_extfunc_t), intent(inout) :: this
      type(c_ext_t), intent(in)         :: cx
      type(c_ext_t)                     :: f_cx
    end function eval_
  end interface

  ! Access specifiers

  private

  public :: c_extfunc_t

  ! Procedures

contains

  subroutine expand_ (this, cx_a, cx_b, f_cx_tol, f_cx_a, f_cx_b, clamp_a, clamp_b, relative_tol)

    class(c_extfunc_t), intent(inout)    :: this
    type(c_ext_t), intent(inout)         :: cx_a
    type(c_ext_t), intent(inout)         :: cx_b
    type(r_ext_t), intent(in)            :: f_cx_tol
    type(c_ext_t), optional, intent(out) :: f_cx_a
    type(c_ext_t), optional, intent(out) :: f_cx_b
    logical, optional, intent(in)        :: clamp_a
    logical, optional, intent(in)        :: clamp_b
    logical, optional, intent(in)        :: relative_tol

    real(WP), parameter :: EXPAND_FACTOR = 1.6_WP

    logical       :: relative_tol_
    logical       :: clamp_a_
    logical       :: clamp_b_
    type(c_ext_t) :: f_a
    type(c_ext_t) :: f_b
    type(r_ext_t) :: tol
    logical       :: move_a

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

    $ASSERT(cx_a /= cx_b,Invalid initial pair)

    ! Expand the pair [cx_a,cx_b] until the difference between f(cx_a)
    ! and f(cx_b) exceeds the tolerance

    f_a = this%eval(cx_a)
    f_b = this%eval(cx_b)

    expand_loop : do

       if (relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + f_cx_tol)*MAX(ABS(f_a), ABS(f_b))
       else
          tol = 4._WP*EPSILON(0._WP)*MAX(ABS(f_a), ABS(f_b)) + f_cx_tol
       endif

       if (ABS(f_a - f_b) > tol) exit expand_loop

       if (clamp_a_) then
          move_a = .FALSE.
       elseif (clamp_b_) then
          move_a = .TRUE.
       else
          move_a = ABS(f_b) > ABS(f_a)
       endif

       if (move_a) then
          cx_a = cx_a + EXPAND_FACTOR*(cx_a - cx_b)
          f_a = this%eval(cx_a)
       else
          cx_b = cx_b + EXPAND_FACTOR*(cx_b - cx_a)
          f_b = this%eval(cx_b)
       endif

    end do expand_loop

    ! Store f_a and f_b

    if (PRESENT(f_cx_a)) f_cx_a = f_a
    if (PRESENT(f_cx_b)) f_cx_b = f_b

    ! Finish

    return

  end subroutine expand_ec_

!****

  subroutine narrow_secant_ (this, cx_a, cx_b, cx_tol, f_cx_a, f_cx_b, n_iter, relative_tol)

    class(c_extfunc_t), intent(inout)             :: this
    type(c_ext_t), intent(inout)           :: cx_a
    type(c_ext_t), intent(inout)           :: cx_b
    type(r_ext_t), intent(in)                 :: cx_tol
    type(c_ext_t), optional, intent(inout) :: f_cx_a
    type(c_ext_t), optional, intent(inout) :: f_cx_b
    integer, optional, intent(inout)             :: n_iter
    logical, optional, intent(in)                :: relative_tol

    integer             :: n_iter_
    logical             :: relative_tol_
    type(c_ext_t) :: a
    type(c_ext_t) :: b
    type(c_ext_t) :: c
    type(c_ext_t) :: f_a
    type(c_ext_t) :: f_b
    type(c_ext_t) :: f_c
    integer             :: i
    type(c_ext_t) :: f_dz
    type(c_ext_t) :: rho
    type(r_ext_t)    :: tol

    if (PRESENT(n_iter)) then
       n_iter_ = n_iter
    else
       n_iter_ = 50
    end if

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    $ASSERT(cx_a /= cx_b,Invalid initial pair)

    ! Use the secant method to narrow the pair [cx_a,cx_b] on the
    ! function root

    ! Set up the initial state

    a = cx_a
    b = cx_b

    if (PRESENT(f_cx_a)) then
       f_a = f_cx_a
    else
       f_a = this%eval(a)
    endif

    if (PRESENT(f_cx_b)) then
       f_b = f_cx_b
    else
       f_b = this%eval(b)
    endif

    if (ABS(f_a) < ABS(f_b)) then

       c = a
       a = b
       b = c

       f_c = f_a
       f_a = f_b
       f_b = f_c

    endif

    ! Iterate until the correction drops below the threshold, or the
    ! maximum number of iterations is exceeded

    iterate_loop : do i = 1,n_iter_

       ! Calculate the correction

       f_dz = f_b*(b - a)

       rho = f_b - f_a

       ! Check for a singular correction

       if (ABS(b*rho) < 8._WP*EPSILON(0._WP)*ABS(f_dz)) then
          $ABORT(Singular correction in secant)
       endif

       ! Update the root

       a = b
       f_a = f_b

       b = b - f_dz/rho
       f_b = this%eval(b)

       ! Check for convergence

       if (relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + cx_tol)*ABS(b)
       else
          tol = 4._WP*EPSILON(0._WP)*ABS(b) + cx_tol
       endif

       if ((ABS(b - a) <= tol .OR. f_b == 0._WP)) exit iterate_loop

    end do iterate_loop

    if (PRESENT(n_iter)) then
       n_iter = i
    else
       $ASSERT(i <= n_iter_,Too many iterations)
    endif

    ! Store the results

    cx_a = a
    cx_b = b

    if (PRESENT(f_cx_a)) f_cx_a = f_a
    if (PRESENT(f_cx_b)) f_cx_b = f_b

    ! Finish

  end subroutine narrow_secant_

!****

  subroutine narrow_ridders_ (this, cx_a, cx_b, rx_tol, f_cx_a, f_cx_b, n_iter, relative_tol)

    class(c_extfunc_t), intent(inout)      :: this
    type(c_ext_t), intent(inout)           :: cx_a
    type(c_ext_t), intent(inout)           :: cx_b
    type(r_ext_t), intent(in)              :: rx_tol
    type(c_ext_t), optional, intent(inout) :: f_cx_a
    type(c_ext_t), optional, intent(inout) :: f_cx_b
    integer, optional, intent(inout)       :: n_iter
    logical, optional, intent(in)          :: relative_tol

    integer       :: n_iter_
    logical       :: relative_tol_
    type(c_ext_t) :: a
    type(c_ext_t) :: b
    type(c_ext_t) :: c
    type(c_ext_t) :: f_a
    type(c_ext_t) :: f_b
    type(c_ext_t) :: f_c
    integer       :: i
    type(c_ext_t) :: exp_Q_p
    type(c_ext_t) :: exp_Q_m
    type(c_ext_t) :: exp_Q
    type(c_ext_t) :: f_dz
    type(c_ext_t) :: rho
    type(r_ext_t) :: tol

    if (PRESENT(n_iter)) then
       n_iter_ = n_iter
    else
       n_iter_ = 50
    end if

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    $ASSERT(cx_a /= cx_b,Invalid initial pair)

    ! Use a complex Ridders' method (with secant updates, rather than
    ! regula falsi) to narrow the pair [cx_a,cx_b] on the function
    ! root

    ! Set up the initial state

    a = cx_a
    b = cx_b

    if (PRESENT(f_cx_a)) then
       f_a = f_cx_a
    else
       f_a = this%eval(a)
    endif

    if (PRESENT(f_cx_b)) then
       f_b = f_cx_b
    else
       f_b = this%eval(b)
    endif

    if (ABS(f_a) < ABS(f_b)) then

       c = a
       a = b
       b = c

       f_c = f_a
       f_a = f_b
       f_b = f_c

    endif

    ! Iterate until the correction drops below the threshold, or the
    ! maximum number of iterations is exceeded

    iterate_loop : do i = 1,n_iter_

       ! Calculate the mid-point values

       c =  0.5_WP*(a + b)
       f_c = this%eval(c)

       ! Solve for the re-scaling exponential

       exp_Q_p = (f_c + SQRT(f_c*f_c - f_a*f_b))/f_b
       exp_Q_m = (f_c - SQRT(f_c*f_c - f_a*f_b))/f_b

       if (ABS(exp_Q_p-1._WP) < ABS(exp_Q_m-1._WP)) then
          exp_Q = exp_Q_p
       else
          exp_Q = exp_Q_m
       endif

       ! Apply the secant method to the re-scaled problem
 
       f_dz = f_b*(exp_Q*exp_Q)*(b - a)

       rho = f_b*(exp_Q*exp_Q) - f_a

       ! Check for a singular correction

       if (ABS(b*rho) < 8._WP*EPSILON(0._WP)*ABS(f_dz)) then
          $ABORT(Singular correction in secant)
       endif

       ! Update the root

       a = b
       f_a = f_b

       b = b - f_dz/rho
       f_b = this%eval(b)

       ! Check for convergence

       if (relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + rx_tol)*ABS(b)
       else
          tol = 4._WP*EPSILON(0._WP)*ABS(b) + rx_tol
       endif

       if ((ABS(b - a) <= tol .OR. f_b == 0._WP)) exit iterate_loop

    end do iterate_loop

    if (PRESENT(n_iter)) then
       n_iter = i
    else
       $ASSERT(i <= n_iter_,Too many iterations)
    endif

    ! Store the results

    cx_a = a
    cx_b = b

    if (PRESENT(f_cx_a)) f_cx_a = f_a
    if (PRESENT(f_cx_b)) f_cx_b = f_b

    ! Finish

  end subroutine narrow_ridders_

!****

  function root_ (this, cx_a, cx_b, rx_tol, f_cx_a, f_cx_b, n_iter, relative_tol) result (cx)

    class(c_extfunc_t), intent(inout)   :: this
    type(c_ext_t), intent(in)           :: cx_a
    type(c_ext_t), intent(in)           :: cx_b
    type(r_ext_t), intent(in)           :: rx_tol
    type(c_ext_t), optional, intent(in) :: f_cx_a
    type(c_ext_t), optional, intent(in) :: f_cx_b
    integer, optional, intent(inout)    :: n_iter
    logical, optional, intent(in)       :: relative_tol
    type(c_ext_t)                       :: cx

    type(c_ext_t) :: a
    type(c_ext_t) :: b
    type(c_ext_t) :: f_a
    type(c_ext_t) :: f_b

    ! Find a function root

    a = cx_a
    b = cx_b

    if (PRESENT(f_cx_a)) then
       f_a = f_cx_a
    else
       f_a = this%eval(a)
    endif

    if (PRESENT(f_cx_b)) then
       f_b = f_cx_b
    else
       f_b = this%eval(b)
    endif

    call this%narrow(a, b, rx_tol, f_a, f_b, n_iter, relative_tol)

    cx = b

    ! Finish

    return

  end function root_

end module gyre_c_extfunc
