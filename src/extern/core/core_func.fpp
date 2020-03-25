! Module   : core_func
! Purpose  : monovariate functions
 
$include 'core.inc'

module core_func

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: func_t
   contains
     private
     procedure                    :: eval_r_
     procedure(eval_c_), deferred :: eval_c_
     generic, public              :: eval => eval_r_, eval_c_
     procedure, public            :: expand_bracket => expand_bracket_
     procedure, public            :: narrow_bracket => narrow_bracket_
     procedure                    :: root_r_
     procedure                    :: root_c_
     generic, public              :: root => root_r_, root_c_
     procedure, public            :: minimum => minimum_r_
     procedure, public            :: maximum => maximum_r_
  end type func_t

  ! Interfaces

  abstract interface
    function eval_c_ (this, z) result (f_z)
      use core_kinds
      import func_t
      class(func_t), intent(inout) :: this
      complex(WP), intent(in)      :: z
      complex(WP)                  :: f_z
    end function eval_c_
  end interface

  ! Access specifiers

  private

  public :: func_t

  ! Procedures

contains

    function eval_r_ (this, x) result (f_x)

    class(func_t), intent(inout) :: this
    real(WP), intent(in)         :: x
    real(WP)                     :: f_x

    ! Evaluate the real function based on the complex function
    ! this%eval(z)

    f_x = REAL(this%eval(CMPLX(x, KIND=WP)), WP)

    ! Finish

    return

  end function eval_r_

!****

  subroutine expand_bracket_ (this, x_a, x_b, f_x_a, f_x_b, clamp_a, clamp_b)

    class(func_t), intent(inout)    :: this
    real(WP), intent(inout)         :: x_a
    real(WP), intent(inout)         :: x_b
    real(WP), intent(out), optional :: f_x_a
    real(WP), intent(out), optional :: f_x_b
    logical, intent(in), optional   :: clamp_a
    logical, intent(in), optional   :: clamp_b

    real(WP), parameter :: EXPAND_FACTOR = 1.6_WP

    logical  :: clamp_a_
    logical  :: clamp_b_
    real(WP) :: f_a
    real(WP) :: f_b
    logical  :: move_a

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

    $ASSERT(x_a /= x_b,Invalid initial bracket)

    ! Expand the bracket [x_a,x_b] until it contains a root of the real function this%eval(x)

    f_a = this%eval(x_a)
    f_b = this%eval(x_b)

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
          x_a = x_a + EXPAND_FACTOR*(x_a - x_b)
          f_a = this%eval(x_a)
       else
          x_b = x_b + EXPAND_FACTOR*(x_b - x_a)
          f_b = this%eval(x_b)
       endif

    end do expand_loop

    ! Store f_a and f_b

    if(PRESENT(f_x_a)) f_x_a = f_a
    if(PRESENT(f_x_b)) f_x_b = f_b

    ! Finish

    return

  end subroutine expand_bracket_

!****

  subroutine narrow_bracket_ (this, x_a, x_b, x_tol, f_x_a, f_x_b, n_iter, relative_tol)

    class(func_t), intent(inout)      :: this
    real(WP), intent(inout)           :: x_a
    real(WP), intent(inout)           :: x_b
    real(WP), intent(in)              :: x_tol
    real(WP), optional, intent(inout) :: f_x_a
    real(WP), optional, intent(inout) :: f_x_b
    integer, optional, intent(inout)  :: n_iter
    logical, optional, intent(in)     :: relative_tol

    real(WP), parameter :: EPS = EPSILON(0._WP)

    integer  :: n_iter_
    logical  :: relative_tol_
    real(WP) :: a
    real(WP) :: b
    real(WP) :: c
    real(WP) :: d
    real(WP) :: e
    real(WP) :: f_a
    real(WP) :: f_b
    real(WP) :: f_c
    real(WP) :: tol
    real(WP) :: m
    real(WP) :: p
    real(WP) :: q
    real(WP) :: r
    real(WP) :: s
    integer  :: i

    if(PRESENT(n_iter)) then
       n_iter_ = n_iter
    else
       n_iter_ = 75
    end if

    if(PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Use Brent's method [based on the ALGOL 60 routine 'zero'
    ! published in Brent (1973, "Algorithms for Minimization without
    ! Derivatives", Prentice Hall, Englewood Cliffs] to narrow the
    ! bracket [x_a,x_b] bracket on the real function this%eval(x)

    ! Set up the initial state

    a = x_a
    b = x_b

    if(PRESENT(f_x_a)) then
       f_a = f_x_a
    else
       f_a = this%eval(a)
    endif

    if(PRESENT(f_x_b)) then
       f_b = f_x_b
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

       if(f_b > 0. .EQV. f_c > 0.) then
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
          tol = (2.*EPS + x_tol)*ABS(b)
       else
          tol = 2.*EPS*ABS(b) + x_tol
       endif

       m = 0.5*(c - b)

       ! Check for convergence

       if(ABS(m) <= tol .OR. f_b == 0.) exit iterate_loop

       ! See if bisection is forced

       if(ABS(e) <  tol .OR. ABS(f_a) < ABS(f_b)) then

          d = m
          e = d

       else

          s = f_b/f_a

          if(a == c) then

             ! Linear interpolation

             p = 2.*m*s
             q = 1. - s

          else

             ! Inverse quadratic interpolation

             q = f_a/f_c
             r = f_b/f_c

             p = s*(2.*m*q*(q - r) - (b - a)*(r - 1.))
             q = (q - 1.)*(r - 1.)*(s - 1.)

          endif

          if(p > 0.) then
             q = -q
          else
             p = -p
          endif

          s = e
          e = d

          if(2.*p < 3.*m*q - ABS(tol*q) .AND. p < ABS(0.5*s*q)) then
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

       b = b + MERGE(d, MERGE(tol, -tol, m > 0), ABS(d) > tol)

       f_b = this%eval(b)

    end do iterate_loop

    if(PRESENT(n_iter)) then
       n_iter = i
    else
       $ASSERT(i <= n_iter_,Too many iterations)
    endif

    ! Store the results

    x_a = a
    x_b = b

    if(PRESENT(f_x_a)) f_x_a = f_a
    if(PRESENT(f_x_b)) f_x_b = f_b

    ! Finish

    return

  end subroutine narrow_bracket_

!****

  function root_r_ (this, x_a, x_b, x_tol, f_x_a, f_x_b, n_iter, relative_tol) result (x)

    class(func_t), intent(inout)     :: this
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    real(WP), intent(in)             :: x_tol
    real(WP), optional, intent(in)   :: f_x_a
    real(WP), optional, intent(in)   :: f_x_b
    integer, optional, intent(inout) :: n_iter
    logical, optional, intent(in)    :: relative_tol
    real(WP)                         :: x

    real(WP) :: a
    real(WP) :: b
    real(WP) :: f_a
    real(WP) :: f_b

    ! Find a root of the real function this%eval(x)

    a = x_a
    b = x_b

    if(PRESENT(f_x_a)) then
       f_a = f_x_a
    else
       f_a = this%eval(a)
    endif

    if(PRESENT(f_x_b)) then
       f_b = f_x_b
    else
       f_b = this%eval(b)
    endif

    call this%narrow_bracket(a, b, x_tol, f_a, f_b, n_iter, relative_tol)

    x = b

    ! Finish

    return

  end function root_r_

!****

  function root_c_ (this, z_a, z_b, z_tol, f_z_a, f_z_b, n_iter, relative_tol) result (z)

    class(func_t), intent(inout)      :: this
    complex(WP), intent(in)           :: z_a
    complex(WP), intent(in)           :: z_b
    real(WP), intent(in)              :: z_tol
    complex(WP), optional, intent(in) :: f_z_a
    complex(WP), optional, intent(in) :: f_z_b
    integer, optional, intent(inout)  :: n_iter
    logical, optional, intent(in)     :: relative_tol
    complex(WP)                       :: z

    real(WP), parameter :: EPS = EPSILON(0._WP)

    integer     :: n_iter_
    logical     :: relative_tol_
    complex(WP) :: a
    complex(WP) :: b
    complex(WP) :: c
    complex(WP) :: f_a
    complex(WP) :: f_b
    complex(WP) :: f_c
    integer     :: i
    complex(WP) :: f_dz
    complex(WP) :: rho
    real(WP)    :: tol

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

    ! Use the secant method to find a complex root of the function
    ! this%eval(z)

    ! Set up the initial state

    a = z_a
    b = z_b

    if(PRESENT(f_z_a)) then
       f_a = f_z_a
    else
       f_a = this%eval(a)
    endif

    if(PRESENT(f_z_b)) then
       f_b = f_z_b
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

       if(ABS(b*rho) < 8._WP*EPS*ABS(f_dz)) then
          $ABORT(Singular correction in secant)
       endif

       ! Update the root

       a = b
       f_a = f_b

       b = b - f_dz/rho
       f_b = this%eval(b)

       ! Check for convergence

       if(relative_tol_) then
          tol = (4._WP*EPS + z_tol)*ABS(b)
       else
          tol = 4._WP*EPS*ABS(b) + z_tol
       endif

       if((ABS(b - a) <= tol .OR. f_b == 0._WP)) exit iterate_loop

    end do iterate_loop

    if(PRESENT(n_iter)) then
       n_iter = i
    else
       $ASSERT(i <= n_iter_,Too many iterations)
    endif

    ! Store the result

    z = b

    ! Finish

  end function root_c_

!****

  function minimum_r_ (this, x_a, x_b, x_c, x_tol, relative_tol) result (x)

    class(func_t), intent(inout)  :: this
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    real(WP), intent(in)          :: x_c
    real(WP), intent(in)          :: x_tol
    logical, intent(in), optional :: relative_tol
    real(WP)                      :: x

    ! Find a local minimum of the real function this%eval(x)
    
    x = extremum_(this, x_a, x_b, x_c, x_tol, .TRUE., relative_tol)

    ! Finish

    return

  end function minimum_r_

!****

  function maximum_r_ (this, x_a, x_b, x_c, x_tol, relative_tol) result (x)

    class(func_t), intent(inout)  :: this
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    real(WP), intent(in)          :: x_c
    real(WP), intent(in)          :: x_tol
    logical, intent(in), optional :: relative_tol
    real(WP)                      :: x

    ! Find a local maximum of the real function this%eval(x)
    
    x = extremum_(this, x_a, x_b, x_c, x_tol, .FALSE., relative_tol)

    ! Finish

    return

  end function maximum_r_

!****

  function extremum_ (rf, x_a, x_b, x_c, x_tol, minimum, relative_tol) result (x)

    class(func_t), intent(inout)  :: rf
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    real(WP), intent(in)          :: x_c
    real(WP), intent(in)          :: x_tol
    logical, intent(in)           :: minimum
    logical, intent(in), optional :: relative_tol
    real(WP)                      :: x

    real(WP), parameter :: EPS = EPSILON(0._WP)
    real(WP), parameter :: R = 0.61803399_WP
    real(WP), parameter :: C = 1._WP - R

    logical  :: relative_tol_
    real(WP) :: x_0
    real(WP) :: x_1
    real(WP) :: x_2
    real(WP) :: x_3
    real(WP) :: f_1
    real(WP) :: f_2
    real(WP) :: tol

    if(PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Use a golden section search to find a local extremum of the real
    ! function this%eval(x)

    ! Set up the starting locations

    x_0 = x_a
    x_3 = x_c

    if(ABS(x_c-x_b) > ABS(x_b-x_a)) then
       x_1 = x_b
       x_2 = x_b + C*(x_c-x_b)
    else
       x_2 = x_b
       x_1 = x_b - C*(x_b-x_a)
    end if

    ! Iterate until convergence

    f_1 = rf%eval(x_1)*MERGE(1._WP, -1._WP, minimum)
    f_2 = rf%eval(x_2)*MERGE(1._WP, -1._WP, minimum)

    iterate_loop : do

       ! Check for a converged bracket

       if(relative_tol_) then
          tol = (SQRT(EPS) + 0.5_WP*x_tol)*(ABS(x_1) + ABS(x_2))
       else
          tol = SQRT(EPS)*(ABS(x_1) + ABS(x_2)) + x_tol
       endif

       if (ABS(x_3-x_0) <= tol) exit iterate_loop

       ! Update the bracket

       if (f_2 < f_1) then
          call shft3_(x_0, x_1, x_2, R*x_2+C*x_3)
          call shft2_(f_1, f_2, rf%eval(x_2)*MERGE(1._WP, -1._WP, minimum))
       else
          call shft3_(x_3, x_2, x_1, R*x_1+C*x_0)
          call shft2_(f_2, f_1, rf%eval(x_1)*MERGE(1._WP, -1._WP, minimum))
       end if

    end do iterate_loop

    ! Select the lower point

    if (f_1 < f_2) then
       x = x_1
    else
       x = x_2
    endif

    ! Finish

    return

  contains

    subroutine shft2_ (a, b, c)

      real(WP), intent(out)   :: a
      real(WP), intent(inout) :: b
      real(WP), intent(in)    :: c

      ! Shift b to a and c to b

      a=b
      b=c

      ! Finish

      return

    end subroutine shft2_

    subroutine shft3_ (a, b, c, d)

      real(WP), intent(out)   :: a 
      real(WP), intent(inout) :: b 
      real(WP), intent(inout) :: c
      real(WP), intent(in)    :: d

      ! Shift b to a, c to b and d to c

      a=b
      b=c
      c=d

    end subroutine shft3_

  end function extremum_

end module core_func
