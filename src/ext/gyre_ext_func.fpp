! Module   : gyre_ext_func
! Purpose  : monovariate functions with extented-range arithmetic
!
! Copyright 2013 Rich Townsend
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

module gyre_ext_func

  ! Uses

  use core_kinds

  use gyre_ext_arith

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ext_func_t
   contains
     procedure                      :: eval_er
     procedure(eval_ec_i), deferred :: eval_ec
     generic                        :: eval => eval_er, eval_ec
     procedure                      :: expand_bracket
     procedure                      :: narrow_bracket
     procedure                      :: narrow_pair
     procedure                      :: root_r
     procedure                      :: root_c
     generic                        :: root => root_r, root_c
  end type ext_func_t 

  ! Interfaces

  abstract interface
    function eval_ec_i (this, ez) result (f_ez)
      use gyre_ext_arith
      import ext_func_t
      class(ext_func_t), intent(inout) :: this
      type(ext_complex_t), intent(in)  :: ez
      type(ext_complex_t)              :: f_ez
    end function eval_ec_i
  end interface

  ! Access specifiers

  private

  public :: ext_func_t

  ! Procedures

contains

  function eval_er (this, ex) result (f_ex)

    class(ext_func_t), intent(inout) :: this
    type(ext_real_t), intent(in)     :: ex
    type(ext_real_t)                 :: f_ex

    ! Evaluate the real function based on the complex function
    ! this%eval_c

    f_ex = ext_real(this%eval(ext_complex(ex)))

    ! Finish

    return

  end function eval_er

!****

  subroutine expand_bracket (this, ex_a, ex_b, f_ex_a, f_ex_b, clamp_a, clamp_b)

    class(ext_func_t), intent(inout)        :: this
    type(ext_real_t), intent(inout)         :: ex_a
    type(ext_real_t), intent(inout)         :: ex_b
    type(ext_real_t), intent(out), optional :: f_ex_a
    type(ext_real_t), intent(out), optional :: f_ex_b
    logical, intent(in), optional           :: clamp_a
    logical, intent(in), optional           :: clamp_b

    real(WP), parameter :: EXPAND_FACTOR = 1.6_WP

    logical          :: clamp_a_
    logical          :: clamp_b_
    type(ext_real_t) :: f_a
    type(ext_real_t) :: f_b
    logical          :: move_a

    if(PRESENT(clamp_a)) then
       clamp_a_ = clamp_a
    else
       clamp_a_ = .FALSE.
    endif

    if(PRESENT(clamp_b)) then
       clamp_b_ = clamp_b
    else
       clamp_b_ = .FALSE.
    endif

    $ASSERT(.NOT. (clamp_a_ .AND. clamp_b_),Cannot clamp both points)

    $ASSERT(ex_a /= ex_b,Invalid initial bracket)

    ! Expand the bracket [x_a,x_b] until it contains a root of the
    ! real function this%eval_r(x)

    f_a = this%eval(ex_a)
    f_b = this%eval(ex_b)

    expand_loop : do

       if((f_a > 0._WP .AND. f_b < 0._WP) .OR. &
          (f_a < 0._WP .AND. f_b > 0._WP)) exit expand_loop

       if(clamp_a_) then
          move_a = .FALSE.
       elseif(clamp_b_) then
          move_a = .TRUE.
       else
          move_a = ABS(f_b) > ABS(f_a)
       endif

       if(move_a) then
          ex_a = ex_a + EXPAND_FACTOR*(ex_a - ex_b)
          f_a = this%eval(ex_a)
       else
          ex_b = ex_b + EXPAND_FACTOR*(ex_b - ex_a)
          f_b = this%eval(ex_b)
       endif

    end do expand_loop

    ! Store f_a and f_b

    if(PRESENT(f_ex_a)) f_ex_a = f_a
    if(PRESENT(f_ex_b)) f_ex_b = f_b

    ! Finish

    return

  end subroutine expand_bracket

!****

  subroutine narrow_bracket (this, ex_a, ex_b, ex_tol, f_ex_a, f_ex_b, n_iter, relative_tol)

    class(ext_func_t), intent(inout)          :: this
    type(ext_real_t), intent(inout)           :: ex_a
    type(ext_real_t), intent(inout)           :: ex_b
    type(ext_real_t), intent(in)              :: ex_tol
    type(ext_real_t), optional, intent(inout) :: f_ex_a
    type(ext_real_t), optional, intent(inout) :: f_ex_b
    integer, optional, intent(inout)          :: n_iter
    logical, optional, intent(in)             :: relative_tol

    integer          :: n_iter_
    logical          :: relative_tol_
    type(ext_real_t) :: a
    type(ext_real_t) :: b
    type(ext_real_t) :: c
    type(ext_real_t) :: d
    type(ext_real_t) :: e
    type(ext_real_t) :: f_a
    type(ext_real_t) :: f_b
    type(ext_real_t) :: f_c
    type(ext_real_t) :: tol
    type(ext_real_t) :: m
    type(ext_real_t) :: p
    type(ext_real_t) :: q
    type(ext_real_t) :: r
    type(ext_real_t) :: s
    integer          :: i

    if(PRESENT(n_iter)) then
       n_iter_ = n_iter
    else
       n_iter_ = 50
    end if

    if(PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Use Brent's method [based on the ALGOL 60 routine 'zero'
    ! published in Brent (1973, "Algorithms for Minimization without
    ! Derivatives", Prentice Hall, Englewood Cliffs] to narrow the
    ! bracket [x_a,x_b] bracket on the real function this%eval_r(x)

    ! Set up the initial state

    a = ex_a
    b = ex_b

    if(PRESENT(f_ex_a)) then
       f_a = f_ex_a
    else
       f_a = this%eval(a)
    endif

    if(PRESENT(f_ex_b)) then
       f_b = f_ex_b
    else
       f_b = this%eval(b)
    endif

    c = b
    f_c = f_b

    ! Check that a root does indeed lie within the bracket

    $ASSERT((f_a >= 0._WP .AND. f_b <= 0._WP) .OR. (f_a <= 0._WP .AND. f_b >= 0._WP),Root is not bracketed)

    ! Iterate until the correction drops below the threshold, or the
    ! maximum number of iterations is exceeded

    iterate_loop : do i = 1,n_iter_

       ! Reorder c so that it has the opposite sign to b

       if(f_b > 0._WP .EQV. f_c > 0._WP) then
          c = a
          f_c = f_a
          d = b - a
          e = d
       endif

       ! Make sure that the function is smallest in magnitude
       ! at b
          
       if(ABS(f_c) < ABS(f_b)) then
          a = b
          b = c
          c = a
          f_a = f_b
          f_b = f_c
          f_c = f_a
       endif

       if(relative_tol_) then
          tol = (2._WP*EPSILON(0._WP) + ex_tol)*ABS(b)
       else
          tol = 2._WP*EPSILON(0._WP)*ABS(b) + ex_tol
       endif

       m = 0.5_WP*(c - b)

       ! Check for convergence

       if(ABS(m) <= tol .OR. f_b == 0._WP) then
          exit iterate_loop
       endif

       ! See if bisection is forced

       if(ABS(e) <  tol .OR. ABS(f_a) < ABS(f_b)) then

          d = m
          e = d

       else

          s = f_b/f_a

          if(a == c) then

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

          if(p > 0._WP) then
             q = -q
          else
             p = -p
          endif

          s = e
          e = d

          if(2._WP*p < 3._WP*m*q - ABS(tol*q) .AND. p < ABS(0.5_WP*s*q)) then
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

       if(ABS(d) > tol) then
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

    if(PRESENT(n_iter)) then
       n_iter = i
    else
       $ASSERT(i <= n_iter_,Too many iterations)
    endif

    ! Store the results

    ex_a = a
    ex_b = b

    if(PRESENT(f_ex_a)) f_ex_a = f_a
    if(PRESENT(f_ex_b)) f_ex_b = f_b

    ! Finish

    return

  end subroutine narrow_bracket

!****

  subroutine narrow_pair (this, ez_a, ez_b, ez_tol, f_ez_a, f_ez_b, n_iter, relative_tol)

    class(ext_func_t), intent(inout)             :: this
    type(ext_complex_t), intent(inout)           :: ez_a
    type(ext_complex_t), intent(inout)           :: ez_b
    type(ext_real_t), intent(in)                 :: ez_tol
    type(ext_complex_t), optional, intent(inout) :: f_ez_a
    type(ext_complex_t), optional, intent(inout) :: f_ez_b
    integer, optional, intent(inout)             :: n_iter
    logical, optional, intent(in)                :: relative_tol

    integer             :: n_iter_
    logical             :: relative_tol_
    type(ext_complex_t) :: a
    type(ext_complex_t) :: b
    type(ext_complex_t) :: c
    type(ext_complex_t) :: f_a
    type(ext_complex_t) :: f_b
    type(ext_complex_t) :: f_c
    integer             :: i
    type(ext_complex_t) :: f_dz
    type(ext_complex_t) :: rho
    type(ext_real_t)    :: tol

    if(PRESENT(n_iter)) then
       n_iter_ = n_iter
    else
       n_iter_ = 50
    end if

    if(PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Use the secant method to narrow the pair [z_a,z_b] on a root of
    ! the complex function this%eval_c(z)

    ! Set up the initial state

    a = ez_a
    b = ez_b

    if(PRESENT(f_ez_a)) then
       f_a = f_ez_a
    else
       f_a = this%eval(a)
    endif

    if(PRESENT(f_ez_b)) then
       f_b = f_ez_b
    else
       f_b = this%eval(b)
    endif

    if(ABS(f_a) < ABS(f_b)) then

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

       if(ABS(b*rho) < 8._WP*EPSILON(0._WP)*ABS(f_dz)) then
          $ABORT(Singular correction in secant)
       endif

       ! Update the root

       a = b
       f_a = f_b

       b = b - f_dz/rho
       f_b = this%eval(b)

       ! Check for convergence

       if(relative_tol_) then
          tol = (4._WP*EPSILON(0._WP) + ez_tol)*ABS(b)
       else
          tol = 4._WP*EPSILON(0._WP)*ABS(b) + ez_tol
       endif

       if((ABS(b - a) <= tol .OR. f_b == 0._WP)) exit iterate_loop

    end do iterate_loop

    if(PRESENT(n_iter)) then
       n_iter = i
    else
       $ASSERT(i <= n_iter_,Too many iterations)
    endif

    ! Store the results

    ez_a = a
    ez_b = b

    if(PRESENT(f_ez_a)) f_ez_a = f_a
    if(PRESENT(f_ez_b)) f_ez_b = f_b

    ! Finish

  end subroutine narrow_pair

!****

  function root_r (this, ex_a, ex_b, ex_tol, f_ex_a, f_ex_b, n_iter, relative_tol) result (ex)

    class(ext_func_t), intent(inout)       :: this
    type(ext_real_t), intent(in)           :: ex_a
    type(ext_real_t), intent(in)           :: ex_b
    type(ext_real_t), intent(in)           :: ex_tol
    type(ext_real_t), optional, intent(in) :: f_ex_a
    type(ext_real_t), optional, intent(in) :: f_ex_b
    integer, optional, intent(inout)       :: n_iter
    logical, optional, intent(in)          :: relative_tol
    type(ext_real_t)                       :: ex

    type(ext_real_t) :: a
    type(ext_real_t) :: b
    type(ext_real_t) :: f_a
    type(ext_real_t) :: f_b

    ! Find a root of the real function this%eval_r(x)

    a = ex_a
    b = ex_b

    if(PRESENT(f_ex_a)) then
       f_a = f_ex_a
    else
       f_a = this%eval(a)
    endif

    if(PRESENT(f_ex_b)) then
       f_b = f_ex_b
    else
       f_b = this%eval(b)
    endif

    call this%narrow_bracket(a, b, ex_tol, f_a, f_b, n_iter, relative_tol)

    ex = b

    ! Finish

    return

  end function root_r

!****

  function root_c (this, ez_a, ez_b, ez_tol, f_ez_a, f_ez_b, n_iter, relative_tol) result (ez)

    class(ext_func_t), intent(inout)          :: this
    type(ext_complex_t), intent(in)           :: ez_a
    type(ext_complex_t), intent(in)           :: ez_b
    type(ext_real_t), intent(in)              :: ez_tol
    type(ext_complex_t), optional, intent(in) :: f_ez_a
    type(ext_complex_t), optional, intent(in) :: f_ez_b
    integer, optional, intent(inout)          :: n_iter
    logical, optional, intent(in)             :: relative_tol
    type(ext_complex_t)                       :: ez

    type(ext_complex_t) :: a
    type(ext_complex_t) :: b
    type(ext_complex_t) :: f_a
    type(ext_complex_t) :: f_b

    ! Find a root of the real function this%eval_r(x)

    a = ez_a
    b = ez_b

    if(PRESENT(f_ez_a)) then
       f_a = f_ez_a
    else
       f_a = this%eval(a)
    endif

    if(PRESENT(f_ez_b)) then
       f_b = f_ez_b
    else
       f_b = this%eval(b)
    endif

    call this%narrow_pair(a, b, ez_tol, f_a, f_b, n_iter, relative_tol)

    ez = b

    ! Finish

    return

  end function root_c

end module gyre_ext_func
