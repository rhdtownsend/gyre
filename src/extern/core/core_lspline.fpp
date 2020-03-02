! Module   : core_lspline
! Purpose  : linear spline interpolation

$include 'core.inc'
$include 'core_parallel.inc'

module core_lspline

  ! Uses

  use core_kinds
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_parallel
  use core_order
  use core_linalg
  use core_interp

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (interp_t) :: lspline_t
     private
     real(WP), allocatable :: knot_x(:)  ! Knot abscissae
     real(WP), allocatable :: knot_y(:)  ! Knot ordinates
     real(WP), allocatable :: knot_dy(:) ! Knot derivatives
   contains
     private
     procedure :: x_n_
     procedure :: y_1_
     procedure :: y_v_
     procedure :: y_n_
     procedure :: dy_1_
     procedure :: dy_v_
     procedure :: dy_n_
     procedure :: iy_1_1_
     procedure :: iy_1_v_
     procedure :: iy_v_1_
     procedure :: iy_v_v_
     procedure :: attribs_
  end type lspline_t

  ! Interfaces

  interface lspline_t
     module procedure lspline_t_
     module procedure lspline_t_eval_derivs_
     module procedure lspline_t_y_func_
  end interface lspline_t

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

  public :: lspline_t
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

  function lspline_t_ (knot_x, knot_y, knot_dy) result (ls)

    real(WP), intent(in) :: knot_x(:)
    real(WP), intent(in) :: knot_y(:)
    real(WP), intent(in) :: knot_dy(:)
    type(lspline_t)      :: ls

    $CHECK_BOUNDS(SIZE(knot_y),SIZE(knot_x))
    $CHECK_BOUNDS(SIZE(knot_dy),SIZE(knot_x))

    ! Construct the lspline_t

    ls%knot_x = knot_x
    ls%knot_y = knot_y
    ls%knot_dy = knot_dy

    ls%n = SIZE(knot_x)

    ! Finish

    return

  end function lspline_t_
  
!****

  function lspline_t_eval_derivs_ (knot_x, knot_y, deriv_type, knot_dy_a, knot_dy_b) result (ls)

    real(WP), intent(in)           :: knot_x(:)
    real(WP), intent(in)           :: knot_y(:)
    character(*), intent(in)       :: deriv_type
    real(WP), optional, intent(in) :: knot_dy_a
    real(WP), optional, intent(in) :: knot_dy_b
    type(lspline_t)                :: ls

    real(WP) :: knot_dy(SIZE(knot_x))

    $CHECK_BOUNDS(SIZE(knot_y),SIZE(knot_x))

    ! Construct the lspline_t, with knot derivatives calculated
    ! according to deriv_type

    select case (deriv_type)
    case('FINDIFF')
       knot_dy = findiff_dy_(knot_x, knot_y, knot_dy_a, knot_dy_b)
    case default
       $ABORT(Invalid deriv_type)
    end select

    ls = lspline_t(knot_x, knot_y, knot_dy)

    ! Finish

    return

  end function lspline_t_eval_derivs_

!****

  function lspline_t_y_func_ (knot_x_a, knot_x_b, y_func, y_tol, deriv_type, relative_tol, knot_dy_a, knot_dy_b) result (ls)

    real(WP), intent(in)           :: knot_x_a
    real(WP), intent(in)           :: knot_x_b
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
    logical, optional, intent(in)  :: relative_tol
    real(WP), optional, intent(in) :: knot_dy_a
    real(WP), optional, intent(in) :: knot_dy_b
    type(lspline_t)                :: ls

    real(WP), parameter :: EPS = 4._WP*EPSILON(0._WP)

    logical               :: relative_tol_
    integer               :: n
    real(WP), allocatable :: knot_x(:)
    real(WP), allocatable :: knot_y(:)
    real(WP), allocatable :: knot_x_dbl(:)
    real(WP), allocatable :: knot_y_dbl(:)
    integer               :: j
    real(WP), allocatable :: err(:)
    real(WP), allocatable :: err_thresh(:)

    if (PRESENT(relative_tol)) then
       relative_tol_ = relative_tol
    else
       relative_tol_ = .FALSE.
    endif

    ! Construct the lspline_t by sampling the function y(x) until the
    ! errors at knots drop below y_tol

    ! Initialize n to the smallest possible value

    n = 2

    knot_x = [knot_x_a,knot_x_b]
    knot_y = [y_func(knot_x_a),y_func(knot_x_b)]

    ! Now increase n until the desired tolerance is reached

    n_loop : do

       ! Construct the lspline

       ls = lspline_t(knot_x, knot_y, deriv_type, knot_dy_a=knot_dy_a, knot_dy_b=knot_dy_b)

       ! Calculate x and y at double the resolution

       allocate(knot_x_dbl(2*n-1))
       allocate(knot_y_dbl(2*n-1))

       knot_x_dbl(1::2) = knot_x
       knot_x_dbl(2::2) = 0.5_WP*(knot_x(:n-1) + knot_x(2:))

       knot_y_dbl(1::2) = knot_y

       do j = 1, n-1
          knot_y_dbl(2*j) = y_func(knot_x_dbl(2*j))
       end do

       ! Examine how well the lspline fits the double-res data

       err = ABS(ls%y(knot_x_dbl(2::2)) - knot_y_dbl(2::2))

       if (relative_tol_) then
          err_thresh = (EPS + y_tol)*ABS(knot_y_dbl(2::2))
       else
          err_thresh = EPS*ABS(knot_y_dbl(2::2)) + y_tol
       endif

       if (ALL(err <= err_thresh)) exit n_loop

       ! Loop around

       call MOVE_ALLOC(knot_x_dbl, knot_x)
       call MOVE_ALLOC(knot_y_dbl, knot_y)

       n = 2*n-1

    end do n_loop
       
    ! Finish

    return

  end function lspline_t_y_func_

!****

  $if ($HDF5)

  subroutine read_ (hg, ls)

    type(hgroup_t), intent(inout) :: hg
    type(lspline_t), intent(out)  :: ls

    real(WP), allocatable :: knot_x(:)
    real(WP), allocatable :: knot_y(:)
    real(WP), allocatable :: knot_dy(:)

    ! Read the lspline_t

    call read_dset_alloc(hg, 'knot_x', knot_x)
    call read_dset_alloc(hg, 'knot_y', knot_y)
    call read_dset_alloc(hg, 'knot_dy', knot_dy)

    ls = lspline_t(knot_x, knot_y, knot_dy)

    ! Finish

    return

  end subroutine read_

!****

  subroutine write_ (hg, ls)

    type(hgroup_t), intent(inout) :: hg
    type(lspline_t), intent(in)   :: ls

    ! Write the lspline_t

    call write_attr(hg, 'n', ls%n)

    call write_dset(hg, 'knot_x', ls%knot_x)
    call write_dset(hg, 'knot_y', ls%knot_y)
    call write_dset(hg, 'knot_dy', ls%knot_dy)

    ! Finish

    return

  end subroutine write_

  $endif

!****

  $if ($MPI)

  subroutine bcast_0_ (ls, root_rank)

    class(lspline_t), intent(inout) :: ls
    integer, intent(in)             :: root_rank

    ! Broadcast the lspline

    call bcast(ls%n, root_rank)

    call bcast_alloc(ls%knot_x, root_rank)
    call bcast_alloc(ls%knot_y, root_rank)
    call bcast_alloc(ls%knot_dy, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST(type(lspline_t),1)
  $BCAST(type(lspline_t),2)
  $BCAST(type(lspline_t),3)
  $BCAST(type(lspline_t),4)

  $BCAST_ALLOC(type(lspline_t),0)
  $BCAST_ALLOC(type(lspline_t),1)
  $BCAST_ALLOC(type(lspline_t),2)
  $BCAST_ALLOC(type(lspline_t),3)
  $BCAST_ALLOC(type(lspline_t),4)

  $endif

!****

  function x_n_ (this) result (x)

    class(lspline_t), intent(in) :: this
    real(WP)                     :: x(this%n)

    ! Return x at knot points

    x = this%knot_x

    ! Finish

    return

  end function x_n_

!****

  function y_1_ (this, x) result (y)

    class(lspline_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP)                     :: y

    real(WP) :: y_v(1)

    ! Interpolate y at a point

    y_v = this%y([x])

    y = y_v(1)

    ! Finish

    return

  end function y_1_

!****

  function y_v_ (this, x) result (y)

    class(lspline_t), intent(in) :: this
    real(WP), intent(in)         :: x(:)
    real(WP)                     :: y(SIZE(x))

    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: t

    ! Interpolate y at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%knot_x, x(j), i)
       $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

       ! Set up the interpolation variables

       h = this%knot_x(i+1) - this%knot_x(i)
       t = (x(j) - this%knot_x(i))/h

       ! Calculate the interpolant

       y(j) = this%knot_y(i  )*phi_(1._WP-t) + &
              this%knot_y(i+1)*phi_(      t)

    end do x_loop

    ! Finish

    return

  contains

    function phi_ (t)

      real(WP), intent(in) :: t
      real(WP)             :: phi_

      phi_ = t

      return

    end function phi_

  end function y_v_

!****

  function y_n_ (this) result (y)

    class(lspline_t), intent(in) :: this
    real(WP)                     :: y(this%n)

    ! Interpolate y at knot points

    y = this%knot_y

    ! Finish

  end function y_n_

!****

  function dy_1_ (this, x) result (dy)

    class(lspline_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP)                     :: dy

    real(WP) :: dy_v(1)

    ! Differentiate y at a point

    dy_v = this%dy([x])

    dy = dy_v(1)

    ! Finish

    return

  end function dy_1_

!****

  function dy_v_ (this, x) result (dy)

    class(lspline_t), intent(in) :: this
    real(WP), intent(in)         :: x(:)
    real(WP)                     :: dy(SIZE(x))

    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: t

    ! Differentiate y at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%knot_x, x(j), i)
       $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

       if (x(j) == this%knot_x(i)) then

          dy(j) = this%knot_dy(i)

       elseif (x(j) == this%knot_x(i+1)) then

          dy(j) = this%knot_dy(i+1)

       else

          ! Set up the interpolation variables

          h = this%knot_x(i+1) - this%knot_x(i)
          t = (x(j) - this%knot_x(i))/h

          ! Calculate the derivative

          dy(j) = -this%knot_y(i  )*dphi_(1._WP-t)/h + &
                   this%knot_y(i+1)*dphi_(      t)/h

       endif

    end do x_loop

    ! Finish

    return

  contains

    function dphi_ (t)

      real(WP), intent(in) :: t
      real(WP)             :: dphi_

      dphi_ = 1._WP

      return

    end function dphi_

  end function dy_v_

!****

  function dy_n_ (this) result (dy)

    class(lspline_t), intent(in) :: this
    real(WP)                     :: dy(this%n)

    ! Differentiate y at knot points

    dy = this%knot_dy

    ! Finish

    return

  end function dy_n_

!****

  function iy_1_1_ (this, x_a, x_b) result (iy)

    class(lspline_t), intent(in) :: this
    real(WP), intent(in)         :: x_a
    real(WP), intent(in)         :: x_b
    real(WP)                     :: iy

    real(WP) :: iy_v(1)

    ! Intergrate y between one point and another

    iy_v = this%iy([x_a], [x_b])

    iy = iy_v(1)

    ! Finish

    return

  end function iy_1_1_

!****

  function iy_1_v_ (this, x_a, x_b) result (iy)

    class(lspline_t), intent(in) :: this
    real(WP), intent(in)         :: x_a
    real(WP), intent(in)         :: x_b(:)
    real(WP)                     :: iy(SIZE(x_b))

    ! Intergrate y between one point and a vector of points

    iy = this%iy(SPREAD(x_a, DIM=1, NCOPIES=SIZE(x_b)), x_b)

    ! Finish

    return

  end function iy_1_v_

!****

  function iy_v_1_ (this, x_a, x_b) result (iy)

    class(lspline_t), intent(in) :: this
    real(WP), intent(in)         :: x_a(:)
    real(WP), intent(in)         :: x_b
    real(WP)                     :: iy(SIZE(x_a))

    ! Intergrate y between a vector of points and a point

    iy = this%iy(x_a, SPREAD(x_b, DIM=1, NCOPIES=SIZE(x_a)))

    ! Finish

    return

  end function iy_v_1_

!****

  function iy_v_v_ (this, x_a, x_b) result (iy)

    class(lspline_t), intent(in) :: this
    real(WP), intent(in)         :: x_a(:)
    real(WP), intent(in)         :: x_b(:)
    real(WP)                     :: iy(SIZE(x_a))

    integer  :: i_a
    integer  :: i_b
    integer  :: j
    integer  :: i
    real(WP) :: h
    real(WP) :: t_1
    real(WP) :: t_2

    $CHECK_BOUNDS(SIZE(x_b),SIZE(x_a))

    ! Integrate y between two vectors of points

    i_a = 1
    i_b = 1

    x_loop : do j = 1,SIZE(x_a)

       ! Update the bracketing indices

       call locate(this%knot_x, x_a(j), i_a)
       $ASSERT(i_a > 0 .AND. i_a < this%n,Out-of-bounds interpolation)

       call locate(this%knot_x, x_b(j), i_b)
       $ASSERT(i_b > 0 .AND. i_b < this%n,Out-of-bounds interpolation)

       ! Calculate the integral

       iy(j) = 0._WP

       do i = i_a, i_b

          ! Set up the interpolation variables

          h = this%knot_x(i+1) - this%knot_x(i)

          if (i == i_a) then
             t_1 = (x_a(j) - this%knot_x(i))/h
          else
             t_1 = 0._WP
          endif

          if (i == i_b) then
             t_2 = (x_b(j) - this%knot_x(i))/h
          else
             t_2 = 1._WP
          endif

          ! Update the sum

          iy(j) = iy(j) - this%knot_y(i  )*(iphi_(1._WP-t_2) - iphi_(1._WP-t_1))*h + &
                          this%knot_y(i+1)*(iphi_(      t_2) - iphi_(      t_1))*h

       end do

    end do x_loop

    ! Finish

    return

  contains

    function iphi_ (t)

      real(WP), intent(in) :: t
      real(WP)             :: iphi_

      iphi_ = t**2/2._WP

      return

    end function iphi_

  end function iy_v_v_

!****

  subroutine attribs_ (this, x_min, x_max, n)

    class(lspline_t), intent(in)    :: this
    real(WP), optional, intent(out) :: x_min
    real(WP), optional, intent(out) :: x_max
    integer, optional, intent(out)  :: n

    ! Return attributes of the lspline

    if (PRESENT(x_min)) x_min = this%knot_x(1)
    if (PRESENT(x_max)) x_max = this%knot_x(this%n)
    if (PRESENT(n)) n = this%n

    ! Finish

    return

  end subroutine attribs_

!****

  function findiff_dy_ (x, y, dy_a, dy_b) result (dy)

    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: y(:)
    real(WP), intent(in), optional :: dy_a
    real(WP), intent(in), optional :: dy_b
    real(WP)                       :: dy(SIZE(x))

    integer  :: n
    real(WP) :: h(SIZE(x)-1)
    real(WP) :: s(SIZE(x)-1)

    $CHECK_BOUNDS(SIZE(y),SIZE(x))

    ! Calculate first derivatives via centered finite differences

    n = SIZE(x)

    h = x(2:) - x(:n-1)

    s = (y(2:) - y(:n-1))/h

    if(PRESENT(dy_a)) then
       dy(1) = dy_a
    else
       dy(1) = s(1)
    endif

    dy(2:n-1) = 0.5_WP*(s(1:n-2) + s(2:n-1))

    if (PRESENT(dy_b)) then
       dy(n) = dy_b
    else
       dy(n) = s(n-1)
    endif

    ! Finish

    return

  end function findiff_dy_

end module core_lspline
