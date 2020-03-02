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
     procedure                     :: eval_r
     procedure(eval_c_i), deferred :: eval_c
     generic                       :: eval => eval_r, eval_c
     procedure                     :: root_r
     procedure                     :: root_c
     generic                       :: root => root_r, root_c
     procedure                     :: minimum => minimum_r
     procedure                     :: maximum => maximum_r
  end type func_t

  ! Interfaces

  abstract interface
    function eval_c_i (this, z) result (f_z)
      use core_kinds
      import func_t
      class(func_t), intent(inout) :: this
      complex(WP), intent(in)      :: z
      complex(WP)                  :: f_z
    end function eval_c_i
  end interface

  ! Access specifiers

  private

  public :: func_t

  ! Procedures

contains

    function eval_r (this, x) result (f_x)

    class(func_t), intent(inout) :: this
    real(WP), intent(in)         :: x
    real(WP)                     :: f_x

    ! Evaluate the real function based on the complex function this%eval_c

    f_x = REAL(this%eval_c(CMPLX(x, KIND=WP)), WP)

    ! Finish

    return

  end function eval_r

!****

  function root_r (this, x_a, x_b, x_tol, f_x_a, f_x_b, n_iter, relative_tol) result (x)

    class(func_t), intent(inout)     :: this
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    real(WP), intent(in)             :: x_tol
    real(WP), optional, intent(in)   :: f_x_a
    real(WP), optional, intent(in)   :: f_x_b
    integer, optional, intent(inout) :: n_iter
    logical, optional, intent(in)    :: relative_tol
    real(WP)                         :: x

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
       n_iter_ = 50
    end if

    if(PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Use Brent's method [based on the ALGOL 60 routine 'zero'
    ! published in Brent (1973, "Algorithms for Minimization without
    ! Derivatives", Prentice Hall, Englewood Cliffs] to find a root of
    ! the real function this%eval_r(x)

    ! Set up the initial state

    a = x_a
    b = x_b

    if(PRESENT(f_x_a)) then
       f_a = f_x_a
    else
       f_a = this%eval_r(a)
    endif

    if(PRESENT(f_x_b)) then
       f_b = f_x_b
    else
       f_b = this%eval_r(b)
    endif

    c = b
    f_c = f_b

    ! Check that a root does indeed lie within the bracket

    $ASSERT((f_a >= 0. .AND. f_b <= 0.) .OR. (f_a <= 0. .AND. f_b >= 0.),Root is not bracketed)

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

       f_b = this%eval_r(b)

    end do iterate_loop

    if(PRESENT(n_iter)) then
       n_iter = i
    else
       $ASSERT(i <= n_iter_,Too many iterations)
    endif

    ! Store the result

    x = b

    ! Finish

    return

  end function root_r

!****

  function root_c (this, z_a, z_b, z_tol, f_z_a, f_z_b, n_iter, relative_tol) result (z)

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
    ! this%eval_c(z)

    ! Set up the initial state

    a = z_a
    b = z_b

    if(PRESENT(f_z_a)) then
       f_a = f_z_a
    else
       f_a = this%eval_c(a)
    endif

    if(PRESENT(f_z_b)) then
       f_b = f_z_b
    else
       f_b = this%eval_c(b)
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
       f_b = this%eval_c(b)

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

  end function root_c

!****

  function minimum_r (this, x_a, x_b, x_c, x_tol, relative_tol) result (x)

    class(func_t), intent(inout)  :: this
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    real(WP), intent(in)          :: x_c
    real(WP), intent(in)          :: x_tol
    logical, intent(in), optional :: relative_tol
    real(WP)                      :: x

    ! Find a local minimum of the real function this%eval_r(x)
    
    x = extremum_r(this, x_a, x_b, x_c, x_tol, .TRUE., relative_tol)

    ! Finish

    return

  end function minimum_r

!****

  function maximum_r (this, x_a, x_b, x_c, x_tol, relative_tol) result (x)

    class(func_t), intent(inout)  :: this
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    real(WP), intent(in)          :: x_c
    real(WP), intent(in)          :: x_tol
    logical, intent(in), optional :: relative_tol
    real(WP)                      :: x

    ! Find a local maximum of the real function this%eval_r(x)
    
    x = extremum_r(this, x_a, x_b, x_c, x_tol, .FALSE., relative_tol)

    ! Finish

    return

  end function maximum_r

!****

  function extremum_r (rf, x_a, x_b, x_c, x_tol, minimum, relative_tol) result (x)

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
    ! function this%eval_r(x)

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

    f_1 = rf%eval_r(x_1)*MERGE(1._WP, -1._WP, minimum)
    f_2 = rf%eval_r(x_2)*MERGE(1._WP, -1._WP, minimum)

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
          call shft3(x_0, x_1, x_2, R*x_2+C*x_3)
          call shft2(f_1, f_2, rf%eval_r(x_2)*MERGE(1._WP, -1._WP, minimum))
       else
          call shft3(x_3, x_2, x_1, R*x_1+C*x_0)
          call shft2(f_2, f_1, rf%eval_r(x_1)*MERGE(1._WP, -1._WP, minimum))
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

    subroutine shft2 (a, b, c)

      real(WP), intent(out)   :: a
      real(WP), intent(inout) :: b
      real(WP), intent(in)    :: c

      ! Shift b to a and c to b

      a=b
      b=c

      ! Finish

      return

    end subroutine shft2

    subroutine shft3 (a,b,c,d)

      real(WP), intent(out)   :: a 
      real(WP), intent(inout) :: b 
      real(WP), intent(inout) :: c
      real(WP), intent(in)    :: d

      ! Shift b to a, c to b and d to c

      a=b
      b=c
      c=d

    end subroutine shft3

  end function extremum_r

end module core_func
