! Module   : gyre_c_root
! Purpose  : root finding algorithms for discriminant functions (complex)
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

module gyre_c_root

  ! Uses

  use core_kinds

  use gyre_cimplex
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
  public :: narrow
  public :: expand

contains

  subroutine solve_ (cf, np, cx_a, cx_b, cx_tol, cx_root, status, n_iter, n_iter_max, relative_tol, f_cx_a, f_cx_b)

    class(c_ext_func_t), intent(inout)  :: cf
    class(num_par_t), intent(in)        :: np
    type(c_ext_t), intent(in)           :: cx_a
    type(c_ext_t), intent(in)           :: cx_b
    type(r_ext_t), intent(in)           :: cx_tol
    type(c_ext_t), intent(out)          :: cx_root
    integer, intent(out)                :: status
    integer, optional, intent(out)      :: n_iter
    integer, optional, intent(in)       :: n_iter_max
    logical, optional, intent(in)       :: relative_tol
    type(c_ext_t), optional, intent(in) :: f_cx_a
    type(c_ext_t), optional, intent(in) :: f_cx_b

    type(c_ext_t) :: a
    type(c_ext_t) :: b
    type(c_ext_t) :: f_a
    type(c_ext_t) :: f_b

    ! Starting from the pair [cx_a,cx_b], find a root of the function
    ! cf

    a = cx_a
    b = cx_b

    if (PRESENT(f_cx_a)) then
       f_a = f_cx_a
    else
       call cf%eval(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_cx_b)) then
       f_b = f_cx_b
    else
       call cf%eval(b, f_b, status)
       if (status /= STATUS_OK) return
    endif

    call narrow(cf, np, a, b, cx_tol, status, n_iter, n_iter_max, relative_tol, f_a, f_b)

    cx_root = b

    ! Finish

    return

  end subroutine solve_

!****

  subroutine narrow_ (cf, np, cx_a, cx_b, cx_tol, status, n_iter, n_iter_max, relative_tol, f_cx_a, f_cx_b)

    class(c_ext_func_t), intent(inout)     :: cf
    class(num_par_t), intent(in)           :: np
    type(c_ext_t), intent(inout)           :: cx_a
    type(c_ext_t), intent(inout)           :: cx_b
    type(r_ext_t), intent(in)              :: cx_tol
    integer, intent(out)                   :: status
    integer, optional, intent(out)         :: n_iter
    integer, optional, intent(in)          :: n_iter_max
    logical, optional, intent(in)          :: relative_tol
    type(c_ext_t), optional, intent(inout) :: f_cx_a
    type(c_ext_t), optional, intent(inout) :: f_cx_b

    ! Narrow the pair [cx_a,cx_b] toward a root of the function cf

    select case (np%c_root_solver)
    case ('SECANT')
       call narrow_secant_(cf, np, cx_a, cx_b, cx_tol, status, n_iter, n_iter_max, relative_tol, f_cx_a, f_cx_b)
    case ('RIDDERS')
       call narrow_ridders_(cf, np, cx_a, cx_b, cx_tol, status, n_iter, n_iter_max, relative_tol, f_cx_a, f_cx_b)
    case ('SIMPLEX')
       call narrow_simplex_(cf, np, cx_a, cx_b, cx_tol, status, n_iter, n_iter_max, relative_tol, f_cx_a, f_cx_b)
    case default
       $ABORT(Invalid c_root_solver)
    end select

    ! Finish

    return

  end subroutine narrow_

!****

  subroutine narrow_secant_ (cf, np, cx_a, cx_b, cx_tol, status, n_iter, n_iter_max, relative_tol, f_cx_a, f_cx_b)
 
    class(c_ext_func_t), intent(inout)     :: cf
    class(num_par_t), intent(in)           :: np
    type(c_ext_t), intent(inout)           :: cx_a
    type(c_ext_t), intent(inout)           :: cx_b
    type(r_ext_t), intent(in)              :: cx_tol
    integer, intent(out)                   :: status
    integer, optional, intent(out)         :: n_iter
    integer, optional, intent(in)          :: n_iter_max
    logical, optional, intent(in)          :: relative_tol
    type(c_ext_t), optional, intent(inout) :: f_cx_a
    type(c_ext_t), optional, intent(inout) :: f_cx_b

    logical       :: relative_tol_
    type(c_ext_t) :: a
    type(c_ext_t) :: b
    type(c_ext_t) :: c
    type(c_ext_t) :: f_a
    type(c_ext_t) :: f_b
    type(c_ext_t) :: f_c
    integer       :: i_iter
    type(c_ext_t) :: f_dz
    type(c_ext_t) :: rho
    type(r_ext_t) :: tol

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Narrow the pair [cx_a,cx_b] toward a root of the function cf
    ! using the secant method

    ! Set up the initial state

    a = cx_a
    b = cx_b

    if (PRESENT(f_cx_a)) then
       f_a = f_cx_a
    else
       call cf%eval(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_cx_b)) then
       f_b = f_cx_b
    else
       call cf%eval(b, f_b, status)
       if (status /= STATUS_OK) return
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

       if (ABS(b*rho) < 8._WP*EPSILON(0._WP)*ABS(f_dz)) then
          $ABORT(Singular correction in secant)
       endif

       ! Update the root

       a = b
       f_a = f_b

       b = b - f_dz/rho
       call cf%eval(b, f_b, status)
       if (status /= STATUS_OK) exit iterate_loop

       ! Check for convergence

       if (relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + cx_tol)*ABS(b)
       else
          tol = 4._WP*EPSILON(0._WP)*ABS(b) + cx_tol
       endif

       if (ABS(b - a) <= tol) exit iterate_loop

    end do iterate_loop

    ! Store the results

    cx_a = a
    cx_b = b

    if (PRESENT(n_iter)) n_iter = i_iter

    if (PRESENT(f_cx_a)) f_cx_a = f_a
    if (PRESENT(f_cx_b)) f_cx_b = f_b

    ! Finish

  end subroutine narrow_secant_

!****

  subroutine narrow_ridders_ (cf, np, cx_a, cx_b, cx_tol, status, n_iter, n_iter_max, relative_tol, f_cx_a, f_cx_b)

    class(c_ext_func_t), intent(inout)     :: cf
    class(num_par_t), intent(in)           :: np
    type(c_ext_t), intent(inout)           :: cx_a
    type(c_ext_t), intent(inout)           :: cx_b
    type(r_ext_t), intent(in)              :: cx_tol
    integer, intent(out)                   :: status
    integer, optional, intent(out)         :: n_iter
    integer, optional, intent(in)          :: n_iter_max
    logical, optional, intent(in)          :: relative_tol
    type(c_ext_t), optional, intent(inout) :: f_cx_a
    type(c_ext_t), optional, intent(inout) :: f_cx_b

    logical       :: relative_tol_
    type(c_ext_t) :: a
    type(c_ext_t) :: b
    type(c_ext_t) :: c
    type(c_ext_t) :: f_a
    type(c_ext_t) :: f_b
    type(c_ext_t) :: f_c
    integer       :: i_iter
    type(c_ext_t) :: exp_Q_p
    type(c_ext_t) :: exp_Q_m
    type(c_ext_t) :: exp_Q
    type(c_ext_t) :: f_dz
    type(c_ext_t) :: rho
    type(r_ext_t) :: tol

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    $ASSERT(cx_a /= cx_b,Invalid initial pair)

    ! Narrow the pair [cx_a,cx_b] toward a root of the function cf
    ! using a complex Ridders' method (with secant updates, rather
    ! than regula falsi)

    ! Set up the initial state

    a = cx_a
    b = cx_b

    if (PRESENT(f_cx_a)) then
       f_a = f_cx_a
    else
       call cf%eval(a, f_a, status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_cx_b)) then
       f_b = f_cx_b
    else
       call cf%eval(b, f_b, status)
       if (status /= STATUS_OK) return
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

       call cf%eval(c, f_c, status)
       if (status /= STATUS_OK) exit iterate_loop

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
       call cf%eval(b, f_b, status)
       if (status /= STATUS_OK) exit iterate_loop

       ! Check for convergence

       if (relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + cx_tol)*ABS(b)
       else
          tol = 4._WP*EPSILON(0._WP)*ABS(b) + cx_tol
       endif

       if (ABS(b - a) <= tol) exit iterate_loop

    end do iterate_loop

    ! Store the results

    cx_a = a
    cx_b = b

    if (PRESENT(n_iter)) n_iter = i_iter

    if (PRESENT(f_cx_a)) f_cx_a = f_a
    if (PRESENT(f_cx_b)) f_cx_b = f_b

    ! Finish

    return

  end subroutine narrow_ridders_

!****

  subroutine narrow_simplex_ (cf, np, cx_a, cx_b, cx_tol, status, n_iter, n_iter_max, relative_tol, f_cx_a, f_cx_b)

    class(c_ext_func_t), intent(inout)     :: cf
    class(num_par_t), intent(in)           :: np
    type(c_ext_t), intent(inout)           :: cx_a
    type(c_ext_t), intent(inout)           :: cx_b
    type(r_ext_t), intent(in)              :: cx_tol
    integer, intent(out)                   :: status
    integer, optional, intent(out)         :: n_iter
    integer, optional, intent(in)          :: n_iter_max
    logical, optional, intent(in)          :: relative_tol
    type(c_ext_t), optional, intent(inout) :: f_cx_a
    type(c_ext_t), optional, intent(inout) :: f_cx_b

    type(c_ext_t)   :: cx(3)
    type(c_ext_t)   :: f_cx(3) 
    type(cimplex_t) :: cm
   
    $ASSERT(cx_a /= cx_b,Invalid initial pair)

    ! Narrow the pair [cx_a,cx_b] toward a root of the function cf
    ! using the simplex algorithm to minimize |cf|

    ! Set up the initial state

    cx(1) = cx_a
    cx(2) = cx_b
    cx(3) = 0.5_WP*(cx(1) + cx(2)) + 0.5_WP*(cx(2) - cx(1))*CMPLX(0._WP, 1._WP, KIND=WP)

    if (PRESENT(f_cx_a)) then
       f_cx(1) = f_cx_a
    else
       call cf%eval(cx(1), f_cx(1), status)
       if (status /= STATUS_OK) return
    endif

    if (PRESENT(f_cx_b)) then
       f_cx(2) = f_cx_b
    else
       call cf%eval(cx(2), f_cx(2), status)
       if (status /= STATUS_OK) return
    endif

    call cf%eval(cx(3), f_cx(3), status)
    if (status /= STATUS_OK) return

    ! Set up the cimplex_t

    cm = cimplex_t(cf, cx, f_cx)

    ! Refine it

    call refine(cm, cx_tol, status, n_iter, n_iter_max, relative_tol)

    ! Store the results

    cx = cm%verts()
    f_cx = cm%values()

    cx_a = cx(2)
    cx_b = cx(1)

    if (PRESENT(f_cx_a)) f_cx_a = f_cx(2)
    if (PRESENT(f_cx_b)) f_cx_b = f_cx(1)

    ! Finish

    return

  end subroutine narrow_simplex_

!****

  subroutine expand_ (cf, cx_a, cx_b, f_cx_tol, status, clamp_a, clamp_b, relative_tol, f_cx_a, f_cx_b)

    class(c_ext_func_t), intent(inout)   :: cf
    type(c_ext_t), intent(inout)         :: cx_a
    type(c_ext_t), intent(inout)         :: cx_b
    type(r_ext_t), intent(in)            :: f_cx_tol
    integer, intent(out)                 :: status
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

    call cf%eval(cx_a, f_a, status)
    if (status /= STATUS_OK) return

    call cf%eval(cx_b, f_b, status)
    if (status /= STATUS_OK) return

    status = STATUS_OK

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
          call cf%eval(cx_a, f_a, status)
          if (status /= STATUS_OK) exit expand_loop
       else
          cx_b = cx_b + EXPAND_FACTOR*(cx_b - cx_a)
          call cf%eval(cx_b, f_b, status)
          if (status /= STATUS_OK) exit expand_loop
       endif

    end do expand_loop

    ! Store f_a and f_b

    if (PRESENT(f_cx_a)) f_cx_a = f_a
    if (PRESENT(f_cx_b)) f_cx_b = f_b

    ! Finish

    return

  end subroutine expand_

end module gyre_c_root
