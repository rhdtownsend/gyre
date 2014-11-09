! Module   : gyre_r_extfunc
! Purpose  : root finding with extended-range arithmetic (real)
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

module gyre_r_extfunc

  ! Uses

  use core_kinds

  use gyre_ext

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: r_extfunc_t
   contains
     private
     procedure(eval_), deferred :: eval
     procedure                  :: expand => expand_
     procedure                  :: narrow => narrow_
     procedure                  :: root => root_
  end type r_extfunc_t

  ! Interfaces

  abstract interface
    function eval_ (this, ex) result (f_ex)
      use gyre_r_ext
      import r_extfunc_t
      class(r_extfunc_t), intent(inout) :: this
      type(r_ext_t), intent(in)         :: rx
      type(r_ext_t)                     :: f_rx
    end function eval_
  end interface

  ! Access specifiers

  private

  public :: r_extfunc_t

  ! Procedures

contains

  subroutine expand_ (this, rx_a, rx_b, f_rx_a, f_rx_b, clamp_a, clamp_b)

    class(r_extfunc_t), intent(inout)    :: this
    type(r_ext_t), intent(inout)         :: rx_a
    type(r_ext_t), intent(inout)         :: rx_b
    type(r_ext_t), optional, intent(out) :: f_rx_a
    type(r_ext_t), optional, intent(out) :: f_rx_b
    logical, optional, intent(in)        :: clamp_a
    logical, optional, intent(in)        :: clamp_b

    real(WP), parameter :: EXPAND_FACTOR = 1.6_WP

    logical          :: clamp_a_
    logical          :: clamp_b_
    type(r_ext_t) :: f_a
    type(r_ext_t) :: f_b
    logical          :: move_a

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

    $ASSERT(rx_a /= rx_b,Invalid initial bracket)

    ! Expand the bracket [rx_a,rx_b] until it contains a function root

    f_a = this%eval(rx_a)
    f_b = this%eval(rx_b)

    expand_loop : do

       if ((f_a > 0._WP .AND. f_b < 0._WP) .OR. &
           (f_a < 0._WP .AND. f_b > 0._WP)) exit expand_loop

       if (clamp_a_) then
          move_a = .FALSE.
       elseif (clamp_b_) then
          move_a = .TRUE.
       else
          move_a = ABS(f_b) > ABS(f_a)
       endif

       if (move_a) then
          rx_a = rx_a + EXPAND_FACTOR*(rx_a - rx_b)
          f_a = this%eval(rx_a)
       else
          rx_b = rx_b + EXPAND_FACTOR*(rx_b - rx_a)
          f_b = this%eval(rx_b)
       endif

    end do expand_loop

    ! Store f_a and f_b

    if (PRESENT(f_rx_a)) f_rx_a = f_a
    if (PRESENT(f_rx_b)) f_rx_b = f_b

    ! Finish

    return

  end subroutine expand_

!****

  subroutine narrow_ (this, ex_a, rx_b, rx_tol, f_rx_a, f_rx_b, n_iter, relative_tol)

    class(r_extfunc_t), intent(inout)      :: this
    type(r_ext_t), intent(inout)           :: rx_a
    type(r_ext_t), intent(inout)           :: rx_b
    type(r_ext_t), intent(in)              :: rx_tol
    type(r_ext_t), optional, intent(inout) :: f_rx_a
    type(r_ext_t), optional, intent(inout) :: f_rx_b
    integer, optional, intent(inout)       :: n_iter
    logical, optional, intent(in)          :: relative_tol

    integer       :: n_iter_
    logical       :: relative_tol_
    type(r_ext_t) :: a
    type(r_ext_t) :: b
    type(r_ext_t) :: c
    type(r_ext_t) :: d
    type(r_ext_t) :: e
    type(r_ext_t) :: f_a
    type(r_ext_t) :: f_b
    type(r_ext_t) :: f_c
    type(r_ext_t) :: tol
    type(r_ext_t) :: m
    type(r_ext_t) :: p
    type(r_ext_t) :: q
    type(r_ext_t) :: r
    type(r_ext_t) :: s
    integer       :: i

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

    ! Use Brent's method [based on the ALGOL 60 routine 'zero'
    ! published in Brent (1973, "Algorithms for Minimization without
    ! Derivatives", Prentice Hall, Englewood Cliffs] to narrow the
    ! bracket [rx_a,rx_b] on the function root

    ! Set up the initial state

    a = rx_a
    b = rx_b

    if (PRESENT(f_rx_a)) then
       f_a = f_rx_a
    else
       f_a = this%eval(a)
    endif

    if (PRESENT(f_rx_b)) then
       f_b = f_rx_b
    else
       f_b = this%eval(b)
    endif

    c = b
    f_c = f_b

    ! Confirm that the bracket contains a root

    $ASSERT((f_a >= 0._WP .AND. f_b <= 0._WP) .OR. (f_a <= 0._WP .AND. f_b >= 0._WP),Root is not bracketed)

    ! Iterate until the correction drops below the threshold, or the
    ! maximum number of iterations is exceeded

    iterate_loop : do i = 1,n_iter_

       ! Reorder c so that it has the opposite sign to b

       if (f_b > 0._WP .EQV. f_c > 0._WP) then
          c = a
          f_c = f_a
          d = b - a
          e = d
       endif

       ! Make sure that the function is smallest in magnitude
       ! at b
          
       if (ABS(f_c) < ABS(f_b)) then
          a = b
          b = c
          c = a
          f_a = f_b
          f_b = f_c
          f_c = f_a
       endif

       if (relative_tol_) then
          tol = (2._WP*EPSILON(0._WP) + rx_tol)*ABS(b)
       else
          tol = 2._WP*EPSILON(0._WP)*ABS(b) + rx_tol
       endif

       m = 0.5_WP*(c - b)

       ! Check for convergence

       if (ABS(m) <= tol .OR. f_b == 0._WP) then
          exit iterate_loop
       endif

       ! See if bisection is forced

       if (ABS(e) <  tol .OR. ABS(f_a) < ABS(f_b)) then

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

          if (2._WP*p < 3._WP*m*q - ABS(tol*q) .AND. p < ABS(0.5_WP*s*q)) then
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

       if (ABS(d) > tol) then
          b = b + d
       else
          if(m > 0._WP) then
             b = b + tol
          else
             b = b - tol
          endif
       endif

       f_b = this%eval(b)

    end do iterate_loop

    if (PRESENT(n_iter)) then
       n_iter = i
    else
       $ASSERT(i <= n_iter_,Too many iterations)
    endif

    ! Store the results

    rx_a = a
    rx_b = b

    if (PRESENT(f_rx_a)) f_rx_a = f_a
    if (PRESENT(f_rx_b)) f_rx_b = f_b

    ! Finish

    return

  end subroutine narrow_

!****

  function root_ (this, rx_a, rx_b, rx_tol, f_rx_a, f_rx_b, n_iter, relative_tol) result (ex)

    class(r_extfunc_t), intent(inout)   :: this
    type(r_ext_t), intent(in)           :: rx_a
    type(r_ext_t), intent(in)           :: rx_b
    type(r_ext_t), intent(in)           :: rx_tol
    type(r_ext_t), optional, intent(in) :: f_rx_a
    type(r_ext_t), optional, intent(in) :: f_rx_b
    integer, optional, intent(inout)    :: n_iter
    logical, optional, intent(in)       :: relative_tol
    type(r_ext_t)                       :: ex

    type(r_ext_t) :: a
    type(r_ext_t) :: b
    type(r_ext_t) :: f_a
    type(r_ext_t) :: f_b

    ! Find a function root

    a = rx_a
    b = rx_b

    if (PRESENT(f_rx_a)) then
       f_a = f_rx_a
    else
       f_a = this%eval(a)
    endif

    if (PRESENT(f_rx_b)) then
       f_b = f_rx_b
    else
       f_b = this%eval(b)
    endif

    call this%narrow(a, b, rx_tol, f_a, f_b, n_iter, relative_tol)

    ex = b

    ! Finish

    return

  end function root_

end module gyre_ext_func
