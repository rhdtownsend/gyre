! Module   : core_spline
! Purpose  : spline interpolation

$include 'core.inc'
$include 'core_parallel.inc'

module core_spline

  ! Uses

  use core_kinds
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_parallel
  use core_order
  use core_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type spline_t
     !private
     real(WP), allocatable :: x(:)     ! Abscissa
     real(WP), allocatable :: y(:)     ! Ordinate
     real(WP), allocatable :: dy_dx(:) ! First derivatives
     integer               :: n        ! Number of points
   contains
     private
     procedure       :: interp_1_
     procedure       :: interp_v_
     procedure       :: interp_n_
     generic, public :: interp => interp_1_, interp_v_, interp_n_
     procedure       :: deriv_1_
     procedure       :: deriv_v_
     procedure       :: deriv_n_
     generic, public :: deriv => deriv_1_, deriv_v_, deriv_n_
     procedure       :: integ_n_
     generic, public :: integ => integ_n_
  end type spline_t

  ! Interfaces

  interface spline_t
     module procedure spline_t_
     module procedure spline_t_eval_derivs_
     module procedure spline_t_y_func_
  end interface spline_t

  $if ($HDF5)
  interface read
     module procedure read_
  end interface read
  interface write
     module procedure write_
  end interface write
  $endif

  $if ($MPI)
 interface bcast
     module procedure bcast_0_
     module procedure bcast_1_
     module procedure bcast_2_
     module procedure bcast_3_
     module procedure bcast_4_
  end interface bcast
  interface bcast_alloc
     module procedure bcast_alloc_0_
     module procedure bcast_alloc_1_
     module procedure bcast_alloc_2_
     module procedure bcast_alloc_3_
     module procedure bcast_alloc_4_
  end interface bcast_alloc
  $endif

  ! Access specifiers

  private

  public :: spline_t
  $if ($HDF5)
  public :: read
  public :: write
  $endif
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  function spline_t_ (x, y, dy_dx) result (sp)

    real(WP), intent(in) :: x(:)
    real(WP), intent(in) :: y(:)
    real(WP), intent(in) :: dy_dx(:)
    type(spline_t)       :: sp
    
    $CHECK_BOUNDS(SIZE(y),SIZE(x))
    $CHECK_BOUNDS(SIZE(dy_dx),SIZE(x))

    ! Construct the spline_t

    sp%x = x
    sp%y = y

    sp%dy_dx = dy_dx

    sp%n = SIZE(x)

    ! Finish

    return

  end function spline_t_
  
!****

  function spline_t_eval_derivs_ (x, y, deriv_type, dy_dx_a, dy_dx_b) result (sp)

    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: y(:)
    character(*), intent(in)       :: deriv_type
    real(WP), optional, intent(in) :: dy_dx_a
    real(WP), optional, intent(in) :: dy_dx_b
    type(spline_t)                 :: sp

    real(WP) :: dy_dx(SIZE(x))

    $CHECK_BOUNDS(SIZE(y),SIZE(x))

    ! Construct the spline_t, with derivatives calculated according to
    ! deriv_type

    select case (deriv_type)
    case('NATURAL')
       dy_dx = natural_dy_dx_(x, y, dy_dx_a, dy_dx_b)
    case('FINDIFF')
       dy_dx = findiff_dy_dx_(x, y, dy_dx_a, dy_dx_b)
    case('MONO')
       dy_dx = mono_dy_dx_(x, y, dy_dx_a, dy_dx_b)
    case default
       $ABORT(Invalid deriv_type)
    end select

    sp = spline_t(x, y, dy_dx)

    ! Finish

    return

  end function spline_t_eval_derivs_

!****

  function spline_t_y_func_ (x_a, x_b, y_func, y_tol, deriv_type, log_samp, relative_tol, dy_dx_a, dy_dx_b) result (sp)

    real(WP), intent(in)           :: x_a
    real(WP), intent(in)           :: x_b
    interface
       function y_func (x) result (y)
         use core_kinds
         implicit none
         real(WP), intent(in) :: x
         real(WP)             :: y
       end function y_func
    end interface
    real(WP), optional, intent(in) :: y_tol
    character(*), intent(in)       :: deriv_type
    logical, optional, intent(in)  :: log_samp
    logical, optional, intent(in)  :: relative_tol
    real(WP), optional, intent(in) :: dy_dx_a
    real(WP), optional, intent(in) :: dy_dx_b
    type(spline_t)                 :: sp

    real(WP), parameter :: EPS = 4._WP*EPSILON(0._WP)

    logical               :: log_samp_
    logical               :: relative_tol_
    integer               :: n
    real(WP), allocatable :: x(:)
    real(WP), allocatable :: y(:)
    real(WP), allocatable :: x_dbl(:)
    real(WP), allocatable :: y_dbl(:)
    integer               :: j
    real(WP), allocatable :: err(:)
    real(WP), allocatable :: err_thresh(:)

    if(PRESENT(log_samp)) then
       log_samp_ = log_samp
    else
       log_samp_ = .FALSE.
    end if
    
    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Construct the spline_t by sampling the function y(x) until the
    ! residuals drop below y_tol

    ! Initialize n to the smallest possible value

    n = 2

    x = [x_a,x_b]
    y = [y_func(x_a),y_func(x_b)]

    ! Now increase n until the desired tolerance is reached

    n_loop : do

       ! Construct the spline

       sp = spline_t(x, y, deriv_type, dy_dx_a=dy_dx_a, dy_dx_b=dy_dx_b)

       ! Calculate x and y at double the resolution

       allocate(x_dbl(2*n-1))
       allocate(y_dbl(2*n-1))

       x_dbl(1::2) = x

       if(log_samp_) then
          x_dbl(2::2) = 10**(0.5_WP*(LOG10(x(:n-1)) + LOG10(x(2:))))
       else ! linear sampling
          x_dbl(2::2) = 0.5_WP*(x(:n-1) + x(2:))
       end if

       y_dbl(1::2) = y
       
       do j = 1, n-1
          y_dbl(2*j) = y_func(x_dbl(2*j))
       end do

       ! Examine how well the spline fits the double-res data

       err = ABS(sp%interp(x_dbl(2::2)) - y_dbl(2::2))

       if (relative_tol_) then
          err_thresh = (EPS + y_tol)*ABS(y_dbl(2::2))
       else
          err_thresh = EPS*ABS(y_dbl(2::2)) + y_tol
       endif

       if (ALL(err <= err_thresh)) exit n_loop

       ! Loop around

       call MOVE_ALLOC(x_dbl, x)
       call MOVE_ALLOC(y_dbl, y)

       n = 2*n-1

    end do n_loop
       
    ! Finish

    return

  end function spline_t_y_func_

!****

  $if ($HDF5)

  subroutine read_ (hg, sp)

    type(hgroup_t), intent(inout) :: hg
    type(spline_t), intent(out)   :: sp

    real(WP), allocatable :: x(:)
    real(WP), allocatable :: y(:)
    real(WP), allocatable :: dy_dx(:)

    ! Read the spline_t

    call read_dset_alloc(hg, 'x', x)
    call read_dset_alloc(hg, 'y', y)
    call read_dset_alloc(hg, 'dy_dx', dy_dx)

    sp = spline_t(x, y, dy_dx)

    ! Finish

    return

  end subroutine read_

!****

  subroutine write_ (hg, sp)

    type(hgroup_t), intent(inout) :: hg
    type(spline_t), intent(in)    :: sp

    ! Write the spline_t

    call write_attr(hg, 'n', sp%n)

    call write_dset(hg, 'x', sp%x)
    call write_dset(hg, 'y', sp%y)
    call write_dset(hg, 'dy_dx', sp%dy_dx)

    ! Finish

    return

  end subroutine write_

  $endif

!****

  $if ($MPI)

  subroutine bcast_0_ (sp, root_rank)

    class(spline_t), intent(inout) :: sp
    integer, intent(in)            :: root_rank

    ! Broadcast the spline

    call bcast(sp%n, root_rank)

    call bcast_alloc(sp%x, root_rank)
    call bcast_alloc(sp%y, root_rank)
    call bcast_alloc(sp%dy_dx, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST(type(spline_t),1)
  $BCAST(type(spline_t),2)
  $BCAST(type(spline_t),3)
  $BCAST(type(spline_t),4)

  $BCAST_ALLOC(type(spline_t),0)
  $BCAST_ALLOC(type(spline_t),1)
  $BCAST_ALLOC(type(spline_t),2)
  $BCAST_ALLOC(type(spline_t),3)
  $BCAST_ALLOC(type(spline_t),4)

  $endif

!****

  function interp_1_ (this, x) result (y)

    class(spline_t), intent(in) :: this
    real(WP), intent(in)        :: x
    real(WP)                    :: y

    real(WP) :: y_n(1)

    ! Interpolate y at a single point

    y_n = this%interp([x])

    y = y_n(1)

    ! Finish

    return

  end function interp_1_

!****

  function interp_v_ (this, x) result (y)

    class(spline_t), intent(in) :: this
    real(WP), intent(in)        :: x(:)
    real(WP)                    :: y(SIZE(x))

    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: w

    ! Interpolate y at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%x, x(j), i)
       $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

       ! Set up the interpolation weights

       h = this%x(i+1) - this%x(i)
       w = (x(j) - this%x(i))/h

       ! Do the interpolation

       y(j) = this%y(i  )*phi_(1._WP-w) + &
              this%y(i+1)*phi_(w      ) - &
        h*this%dy_dx(i  )*psi_(1._WP-w) + &
        h*this%dy_dx(i+1)*psi_(w      )

    end do x_loop

    ! Finish

    return

  contains

    function phi_ (t)

      real(WP), intent(in) :: t
      real(WP)             :: phi_

      phi_ = 3._WP*t**2 - 2._WP*t**3

      return

    end function phi_

    function psi_ (t)

      real(WP), intent(in) :: t
      real(WP)             :: psi_

      psi_ = t**3 - t**2

      return

    end function psi_

  end function interp_v_

!****

  function interp_n_ (this) result (y)

    class(spline_t), intent(in) :: this
    real(WP)                    :: y(this%n)

    ! Interpolate y at abscissa points

    y = this%y

    ! Finish

  end function interp_n_

!****

  function deriv_1_ (this, x) result (dy_dx)

    class(spline_t), intent(in) :: this
    real(WP), intent(in)        :: x
    real(WP)                    :: dy_dx

    real(WP) :: dy_dx_n(1)

    ! Differentiate y at a single point

    dy_dx_n = this%deriv([x])

    dy_dx = dy_dx_n(1)

    ! Finish

    return

  end function deriv_1_

!****

  function deriv_v_ (this, x) result (dy_dx)

    class(spline_t), intent(in) :: this
    real(WP), intent(in)        :: x(:)
    real(WP)                    :: dy_dx(SIZE(x))

    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: w

    ! Differentiate y at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%x, x(j), i)
       $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

       ! Set up the interpolation weights

       h = this%x(i+1) - this%x(i)
       w = (x(j) - this%x(i))/h

       ! Do the interpolation

       dy_dx(j) = -this%y(i  )*dphi_dt_(1._WP-w)/h + &
                   this%y(i+1)*dphi_dt_(w      )/h + &
               this%dy_dx(i  )*dpsi_dt_(1._WP-w) + &
               this%dy_dx(i+1)*dpsi_dt_(w      )

    end do x_loop

    ! Finish

    return

  contains

    function dphi_dt_ (t)

      real(WP), intent(in) :: t
      real(WP)             :: dphi_dt_

      dphi_dt_ = 6._WP*t - 6._WP*t**2

      return

    end function dphi_dt_

    function dpsi_dt_ (t)

      real(WP), intent(in) :: t
      real(WP)             :: dpsi_dt_

      dpsi_dt_ = 3._WP*t**2 - 2._WP*t

      return

    end function dpsi_dt_

  end function deriv_v_

!****

  function deriv_n_ (this) result (dy_dx)

    class(spline_t), intent(in) :: this
    real(WP)                    :: dy_dx(this%n)

    ! Differentiate y at abscissa points

    dy_dx = this%dy_dx

    ! Finish

    return

  end function deriv_n_

!****

  function integ_n_ (this) result (Y)

    class(spline_t), intent(in) :: this
    real(WP)                    :: Y(this%n)

    integer  :: i
    real(WP) :: h

    ! Calculate the integral of y

    Y(1) = 0._WP

    x_loop : do i = 1,this%n-1
       
       h = this%x(i+1) - this%x(i)

       Y(i+1) = Y(i) + (this%y(i) + this%y(i+1))*h/2._WP - &
                       (this%dy_dx(i+1) - this%dy_dx(i))*h**2/12._WP

    end do x_loop

    ! Finish

    return

  end function integ_n_

!****

  function natural_dy_dx_ (x, y, dy_dx_a, dy_dx_b) result (dy_dx)

    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: y(:)
    real(WP), intent(in), optional :: dy_dx_a
    real(WP), intent(in), optional :: dy_dx_b
    real(WP)                       :: dy_dx(SIZE(x))

    integer  :: n
    real(WP) :: h(SIZE(x)-1)
    real(WP) :: L(SIZE(x)-1)
    real(WP) :: D(SIZE(x))
    real(WP) :: U(SIZE(x)-1)
    real(WP) :: B(SIZE(x),1)
    integer  :: info

    $CHECK_BOUNDS(SIZE(y),SIZE(x))
    
    ! Calcualte the first derivatives for a natural spline (ensuring
    ! the second derivatives are continuous)

    n = SIZE(x)

    h = x(2:) - x(:n-1)

    ! Set up the tridiagonal matrix and RHS

    ! Inner boundary

    D(1) = 1._WP
    U(1) = 0._WP

    if(PRESENT(dy_dx_a)) then
       B(1,1) = dy_dx_a
    else
       B(1,1) = (y(2) - y(1))/h(1)
    endif

    ! Internal points

    L(1:n-2) = 2._WP/h(1:n-2)
    D(2:n-1) = 4._WP/h(1:n-2) + 4._WP/h(2:n-1)
    U(2:n-1) = 2._WP/h(2:n-1)

    B(2:n-1,1) = -6._WP*y(1:n-2)/h(1:n-2)**2 + 6._WP*y(2:n-1)/h(1:n-2)**2 + &
                 6._WP*y(3:n  )/h(2:n-1)**2 - 6._WP*y(2:n-1)/h(2:n-1)**2

    ! Outer boundary

    L(n-1) = 0._WP
    D(n) = 1._WP

    if(PRESENT(dy_dx_b)) then
       B(n,1) = dy_dx_b
    else
       B(n,1) = (y(n) - y(n-1))/h(n-1)
    endif

    ! Solve the tridiagonal system

    call XGTSV(n, 1, L, D, U, B, SIZE(B, 1), info)
    $ASSERT(info == 0,Non-zero return from XTGSV)

    dy_dx = B(:,1)

    ! Finish

    return

  end function natural_dy_dx_

!****

  function findiff_dy_dx_ (x, y, dy_dx_a, dy_dx_b) result (dy_dx)

    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: y(:)
    real(WP), intent(in), optional :: dy_dx_a
    real(WP), intent(in), optional :: dy_dx_b
    real(WP)                       :: dy_dx(SIZE(x))

    integer  :: n
    real(WP) :: h(SIZE(x)-1)
    real(WP) :: s(SIZE(x)-1)

    $CHECK_BOUNDS(SIZE(y),SIZE(x))

    ! Calculate the first derivatives via centered finite differences

    n = SIZE(x)

    h = x(2:) - x(:n-1)

    s = (y(2:) - y(:n-1))/h

    if(PRESENT(dy_dx_a)) then
       dy_dx(1) = dy_dx_a
    else
       dy_dx(1) = s(1)
    endif

    dy_dx(2:n-1) = 0.5_WP*(s(1:n-2) + s(2:n-1))

    if(PRESENT(dy_dx_b)) then
       dy_dx(n) = dy_dx_b
    else
       dy_dx(n) = s(n-1)
    endif

    ! Finish

    return

  end function findiff_dy_dx_

!****

  function mono_dy_dx_ (x, y, dy_dx_a, dy_dx_b) result (dy_dx)

    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: y(:)
    real(WP), intent(in), optional :: dy_dx_a
    real(WP), intent(in), optional :: dy_dx_b
    real(WP)                       :: dy_dx(SIZE(x))

    integer  :: n
    real(WP) :: h(SIZE(x)-1)
    real(WP) :: s(SIZE(x)-1)
    real(WP) :: p(SIZE(x))
    integer  :: i

    $CHECK_BOUNDS(SIZE(y),SIZE(x))

    ! Calculate the first derivatives using the Steffen (1990, A&A,
    ! 239, 443) monontonicity preserving algorithm

    n = SIZE(x)

    h = x(2:) - x(:n-1)

    s = (y(2:) - y(:n-1))/h

    ! Calculate parabolic gradients

    if(PRESENT(dy_dx_a)) then
       p(1) = dy_dx_a
    else
       p(1) = s(1)
    endif

    p(2:n-1) = (s(1:n-2)*h(2:n-1) + s(2:n-1)*h(1:n-2))/(h(1:n-2) + h(2:n-1))

    if(PRESENT(dy_dx_b)) then
       p(n) = dy_dx_b
    else
       p(n) = s(n-1)
    endif

    ! Calculate monotonic gradients

    dy_dx(1) = p(1)

    do i = 2,n-1
       dy_dx(i) = (SIGN(1._WP, s(i-1)) + SIGN(1._WP, s(i)))* &
                  MIN(ABS(s(i-1)), ABS(s(i)), 0.5*ABS(p(i)))
    end do

    dy_dx(n) = p(n)

    ! Finish

    return

  end function mono_dy_dx_

end module core_spline
