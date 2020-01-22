! Module   : gyre_r_min
! Purpose  : minimum finding algorithms (real)
!
! Copyright 2018-2020 Rich Townsend & The GYRE Team
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

module gyre_r_min

  ! Uses

  use core_kinds

  use gyre_ext_func
  use gyre_ext
  use gyre_math
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

  ! Access specifiers

  private

  public :: solve
!  public :: narrow

contains

  subroutine solve_ (rf, rx_a, rx_b, rx_c, rx_tol, nm_p, rx_min, status, n_iter, n_iter_max, relative_tol, f_rx_a, f_rx_b, f_rx_c)

    class(r_ext_func_t), intent(inout)  :: rf
    type(r_ext_t), intent(in)           :: rx_a
    type(r_ext_t), intent(in)           :: rx_b
    type(r_ext_t), intent(in)           :: rx_c
    type(r_ext_t), intent(in)           :: rx_tol
    type(r_ext_t), intent(out)          :: rx_min
    class(num_par_t), intent(in)        :: nm_p
    integer, intent(out)                :: status
    integer, optional, intent(out)      :: n_iter
    integer, optional, intent(in)       :: n_iter_max
    logical, optional, intent(in)       :: relative_tol
    type(r_ext_t), optional, intent(in) :: f_rx_a
    type(r_ext_t), optional, intent(in) :: f_rx_b
    type(r_ext_t), optional, intent(in) :: f_rx_c

    type(r_ext_t) :: a
    type(r_ext_t) :: b
    type(r_ext_t) :: c
    type(r_ext_t) :: f_a
    type(r_ext_t) :: f_b
    type(r_ext_t) :: f_c

    ! Starting from the bracket [rx_a,rx_b,rx_c], find a minimum of
    ! the function rf

    a = rx_a
    b = rx_b
    c = rx_c

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

    if (PRESENT(f_rx_c)) then
       f_c = f_rx_c
    else
       call rf%eval(c, f_c, status)
       if (status /= STATUS_OK) return
    endif

    call narrow_(rf, a, b, c, rx_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_a, f_b, f_c)

    rx_min = b

    ! Finish

    return

  end subroutine solve_

  !****

  subroutine narrow_ (rf, rx_a, rx_b, rx_c, rx_tol, nm_p, status, n_iter, n_iter_max, relative_tol, f_rx_a, f_rx_b, f_rx_c)

    class(r_ext_func_t), intent(inout)     :: rf
    type(r_ext_t), intent(inout)           :: rx_a
    type(r_ext_t), intent(inout)           :: rx_b
    type(r_ext_t), intent(inout)           :: rx_c
    type(r_ext_t), intent(in)              :: rx_tol
    class(num_par_t), intent(in)           :: nm_p
    integer, intent(out)                   :: status
    integer, optional, intent(out)         :: n_iter
    integer, optional, intent(in)          :: n_iter_max
    logical, optional, intent(in)          :: relative_tol
    type(r_ext_t), optional, intent(inout) :: f_rx_a
    type(r_ext_t), optional, intent(inout) :: f_rx_b
    type(r_ext_t), optional, intent(inout) :: f_rx_c

    ! Narrow the bracket [rx_a,rx_b,rx_c] on a minimum of the function rf

    call narrow_bisect_(rf, rx_a, rx_b, rx_c, rx_tol, status, n_iter, n_iter_max, relative_tol, f_rx_a, f_rx_b, f_rx_c)

    ! Finish

    return

  end subroutine narrow_

  !****

  subroutine narrow_bisect_ (rf, rx_a, rx_b, rx_c, rx_tol, status, n_iter, n_iter_max, relative_tol, f_rx_a, f_rx_b, f_rx_c)

    class(r_ext_func_t), intent(inout)     :: rf
    type(r_ext_t), intent(inout)           :: rx_a
    type(r_ext_t), intent(inout)           :: rx_b
    type(r_ext_t), intent(inout)           :: rx_c
    type(r_ext_t), intent(in)              :: rx_tol
    integer, intent(out)                   :: status
    integer, optional, intent(out)         :: n_iter
    integer, optional, intent(in)          :: n_iter_max
    logical, optional, intent(in)          :: relative_tol
    type(r_ext_t), optional, intent(inout) :: f_rx_a
    type(r_ext_t), optional, intent(inout) :: f_rx_b
    type(r_ext_t), optional, intent(inout) :: f_rx_c

    real(WP), parameter :: R = 0.61803399_WP

    logical       :: relative_tol_
    type(r_ext_t) :: a
    type(r_ext_t) :: b
    type(r_ext_t) :: c
    type(r_ext_t) :: d
    type(r_ext_t) :: f_a
    type(r_ext_t) :: f_b
    type(r_ext_t) :: f_c
    type(r_ext_t) :: f_d
    integer       :: i_iter
    type(r_ext_t) :: tol

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Narrow the bracket [rx_a,rx_b,rx_c] on a minimum of the function rf
    ! using a bisection method

    ! Set up the initial state

    a = rx_a
    b = rx_b
    c = rx_c

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

    if (PRESENT(f_rx_c)) then
       f_c = f_rx_c
    else
       call rf%eval(c, f_c, status)
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
          tol = (sqrt(EPSILON(0._WP)) + 0.5_WP*rx_tol)*(abs(a) + abs(c))
       else
          tol = sqrt(EPSILON(0._WP))*(abs(a) + abs(c)) + rx_tol
       endif

       if (abs(c-a) <= tol) exit iterate_loop

       ! Create a new point, and update the bracket accordingly

       if (b-a > c-b) then

          ! New point in left interval
          
          d = (1._WP-R)*b + R*a

          call rf%eval(d, f_d, status)
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

          call rf%eval(d, f_d, status)
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

    rx_a = a
    rx_b = b
    rx_c = c

    if (PRESENT(n_iter)) then
       n_iter = i_iter
    end if

    if (PRESENT(f_rx_a)) f_rx_a = f_a
    if (PRESENT(f_rx_b)) f_rx_b = f_b
    if (PRESENT(f_rx_c)) f_rx_c = f_c

    ! Finish

    return

  end subroutine narrow_bisect_

end module gyre_r_min
