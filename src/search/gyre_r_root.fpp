! Module   : gyre_r_root
! Purpose  : root finding algorithms for discriminant functions (real)
!
! Copyright 2013-2015 Rich Townsend
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

module gyre_r_root

  ! Uses

  use core_kinds

  use gyre_ext_func
  use gyre_ext
  use gyre_num_par
  use gyre_status

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface solve
     module procedure solve_
  end interface solve

  interface narrow
     module procedure narrow_
  end interface narrow

  interface expand
     module procedure expand_
  end interface expand

  ! Access specifiers

  private

  public :: solve
  public :: expand
  public :: narrow

contains

  subroutine solve_ (rf, rx_a, rx_b, rx_tol, nm_p, rx_root, status, n_iter, n_iter_max, relative_tol, f_rx_a, f_rx_b)

    class(r_ext_func_t), intent(inout)  :: rf
    type(r_ext_t), intent(in)           :: rx_a
    type(r_ext_t), intent(in)           :: rx_b
    type(r_ext_t), intent(in)           :: rx_tol
    type(r_ext_t), intent(out)          :: rx_root
    class(num_par_t), intent(in)        :: nm_p
    integer, intent(out)                :: status
    integer, optional, intent(out)      :: n_iter
    integer, optional, intent(in)       :: n_iter_max
    logical, optional, intent(in)       :: relative_tol
    type(r_ext_t), optional, intent(in) :: f_rx_a
    type(r_ext_t), optional, intent(in) :: f_rx_b

    type(r_ext_t) :: a
    type(r_ext_t) :: b
    type(r_ext_t) :: f_a
    type(r_ext_t) :: f_b

    ! Starting from the bracket [rx_a,rx_b], find a root of the
    ! function rf

    a = rx_a
    b = rx_b

    if (PRESENT(f_rx_a)) then
       f_a = f_rx_a
    else
       call rf%eval(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_rx_b)) then
       f_b = f_rx_b
    else
       call rf%eval(b, f_b, status)
       if (status /= STATUS_OK) return
    endif

    call narrow_(rf, a, b, rx_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_a, f_b)

    rx_root = b

    ! Finish

    return

  end subroutine solve_

!****

  subroutine narrow_ (rf, rx_a, rx_b, rx_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_rx_a, f_rx_b)

    class(r_ext_func_t), intent(inout)     :: rf
    type(r_ext_t), intent(inout)           :: rx_a
    type(r_ext_t), intent(inout)           :: rx_b
    type(r_ext_t), intent(in)              :: rx_tol
    class(num_par_t), intent(in)           :: nm_p
    integer, intent(out)                   :: status
    integer, optional, intent(out)         :: n_iter
    integer, optional, intent(in)          :: n_iter_max
    logical, optional, intent(in)          :: relative_tol
    type(r_ext_t), optional, intent(inout) :: f_rx_a
    type(r_ext_t), optional, intent(inout) :: f_rx_b

    ! Narrow the bracket [rx_a,rx_b] on a root of the function rf

    select case (nm_p%r_root_solver)
    case ('BRENT')
       call narrow_brent_(rf, rx_a, rx_b, rx_tol, status, n_iter, n_iter_max, relative_tol, f_rx_a, f_rx_b)
    case default
       $ABORT(Invalid r_root_solver)
    end select

    ! Finish

    return

  end subroutine narrow_

!****

  subroutine narrow_brent_ (rf, rx_a, rx_b, rx_tol, status, n_iter, n_iter_max, relative_tol, f_rx_a, f_rx_b)

    class(r_ext_func_t), intent(inout)     :: rf
    type(r_ext_t), intent(inout)           :: rx_a
    type(r_ext_t), intent(inout)           :: rx_b
    type(r_ext_t), intent(in)              :: rx_tol
    integer, intent(out)                   :: status
    integer, optional, intent(out)         :: n_iter
    integer, optional, intent(in)          :: n_iter_max
    logical, optional, intent(in)          :: relative_tol
    type(r_ext_t), optional, intent(inout) :: f_rx_a
    type(r_ext_t), optional, intent(inout) :: f_rx_b

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
    integer       :: i_iter

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Narrow the bracket [rx_a,rx_b] on a root of the function rf
    ! using Brent's method [based on the ALGOL 60 routine 'zero'
    ! published in [Bre1973]

    ! Set up the initial state

    a = rx_a
    b = rx_b

    if (PRESENT(f_rx_a)) then
       f_a = f_rx_a
    else
       call rf%eval(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_rx_b)) then
       f_b = f_rx_b
    else
       call rf%eval(b, f_b, status)
       if (status /= STATUS_OK) return
    endif

    c = b
    f_c = f_b

    ! Confirm that the bracket contains a root

    $ASSERT((f_a >= 0._WP .AND. f_b <= 0._WP) .OR. (f_a <= 0._WP .AND. f_b >= 0._WP),Root is not bracketed)

    ! Iterate until the correction drops below the threshold, or the
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

       call rf%eval(b, f_b, status)
       if (status /= STATUS_OK) exit iterate_loop
       
    end do iterate_loop

    ! Store the results

    rx_a = a
    rx_b = b

    if (PRESENT(n_iter)) then
       n_iter = i_iter
    end if

    if (PRESENT(f_rx_a)) f_rx_a = f_a
    if (PRESENT(f_rx_b)) f_rx_b = f_b

    ! Finish

    return

  end subroutine narrow_brent_

!****

  subroutine expand_ (rf, rx_a, rx_b, status, clamp_a, clamp_b, f_rx_a, f_rx_b)

    class(r_ext_func_t), intent(inout)   :: rf
    type(r_ext_t), intent(inout)         :: rx_a
    type(r_ext_t), intent(inout)         :: rx_b
    integer, intent(out)                 :: status
    logical, optional, intent(in)        :: clamp_a
    logical, optional, intent(in)        :: clamp_b
    type(r_ext_t), optional, intent(out) :: f_rx_a
    type(r_ext_t), optional, intent(out) :: f_rx_b

    real(WP), parameter :: EXPAND_FACTOR = 1.6_WP

    logical       :: clamp_a_
    logical       :: clamp_b_
    type(r_ext_t) :: f_a
    type(r_ext_t) :: f_b
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

    $ASSERT(rx_a /= rx_b,Invalid initial bracket)

    ! Expand the bracket [rx_a,rx_b] until it contains a root of the
    ! function rf

    call rf%eval(rx_a, f_a, status)
    if (status /= STATUS_OK) return

    call rf%eval(rx_b, f_b, status)
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
          move_a = ABS(f_b) > ABS(f_a)
       endif

       if (move_a) then

          rx_a = rx_a + EXPAND_FACTOR*(rx_a - rx_b)

          call rf%eval(rx_a, f_a, status)
          if (status /= STATUS_OK) exit expand_loop

       else

          rx_b = rx_b + EXPAND_FACTOR*(rx_b - rx_a)

          call rf%eval(rx_b, f_b, status)
          if (status /= STATUS_OK) exit expand_loop

       endif

    end do expand_loop

    ! Store the results

    if (PRESENT(f_rx_a)) f_rx_a = f_a
    if (PRESENT(f_rx_b)) f_rx_b = f_b

    ! Finish

    return

  end subroutine expand_

end module gyre_r_root
