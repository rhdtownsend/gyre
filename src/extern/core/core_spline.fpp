! Module   : core_spline
! Purpose  : spline interpolation

$include 'core.inc'

module core_spline

  ! Uses

  use core_kinds
  $if($HDF5)
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
     private
     real(WP), allocatable :: x(:)     ! Abscissa
     real(WP), allocatable :: y(:)     ! Ordinate
     real(WP), allocatable :: dy_dx(:) ! First derivatives
     integer               :: n        ! Number of points
   contains
     private
     procedure         :: init_base
     procedure         :: init_eval_derivs
     generic, public   :: init => init_base, init_eval_derivs
     $if($GFORTRAN_PR57922)
     procedure, public :: final
     $endif
     $if($HDF5)
     procedure, public :: read
     procedure, public :: write
     $endif
     procedure         :: interp_1
     procedure         :: interp_v
     procedure         :: interp_n
     generic, public   :: interp => interp_1, interp_v, interp_n
     procedure         :: deriv_1
     procedure         :: deriv_v
     procedure         :: deriv_n
     generic, public   :: deriv => deriv_1, deriv_v, deriv_n
     procedure         :: integ_n
     generic, public   :: integ => integ_n
  end type spline_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_sp
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: spline_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init_base (this, x, y, dy_dx)

    class(spline_t), intent(out) :: this
    real(WP), intent(in)         :: x(:)
    real(WP), intent(in)         :: y(:)
    real(WP), intent(in)         :: dy_dx(:)

    $CHECK_BOUNDS(SIZE(y),SIZE(x))
    $CHECK_BOUNDS(SIZE(dy_dx),SIZE(x))

    ! Initialize the spline

    this%x = x
    this%y = y

    this%dy_dx = dy_dx

    this%n = SIZE(x)

    ! Finish

    return

  end subroutine init_base

!****

  subroutine init_eval_derivs (this, x, y, deriv_type, dy_dx_a, dy_dx_b)

    class(spline_t), intent(out)   :: this
    real(WP), intent(in)           :: x(:)
    real(WP), intent(in)           :: y(:)
    character(LEN=*), intent(in)   :: deriv_type
    real(WP), intent(in), optional :: dy_dx_a
    real(WP), intent(in), optional :: dy_dx_b

    real(WP) :: dy_dx(SIZE(x))

    $CHECK_BOUNDS(SIZE(y),SIZE(x))

    ! Initialize the spline

    select case (deriv_type)
    case('NATURAL')
       dy_dx = natural_dy_dx(x, y, dy_dx_a, dy_dx_b)
    case('FINDIFF')
       dy_dx = findiff_dy_dx(x, y, dy_dx_a, dy_dx_b)
    case('MONO')
       dy_dx = mono_dy_dx(x, y, dy_dx_a, dy_dx_b)
    case default
       $ABORT(Invalid deriv_type)
    end select

    call this%init(x, y, dy_dx)

    ! Finish

    return

  end subroutine init_eval_derivs

!****

  $if($GFORTRAN_PR57922)

  subroutine final (this)

    class(spline_t), intent(inout) :: this

    ! Finalize the spline

    deallocate(this%x)
    deallocate(this%y)
    deallocate(this%dy_dx)

    ! Finish

    return

  end subroutine final

  $endif

!****

  $if($HDF5)

  subroutine read (this, hg)

    class(spline_t), intent(out)  :: this
    type(hgroup_t), intent(inout) :: hg

    real(WP), allocatable :: x(:)
    real(WP), allocatable :: y(:)
    real(WP), allocatable :: dy_dx(:)

    ! Read the spline

    call read_dset_alloc(hg, 'x', x)
    call read_dset_alloc(hg, 'y', y)
    call read_dset_alloc(hg, 'dy_dx', dy_dx)

    call this%init(x, y, dy_dx)

    ! Finish

    return

  end subroutine read

!****

  subroutine write (this, hg)

    class(spline_t), intent(in)   :: this
    type(hgroup_t), intent(inout) :: hg

    ! Write the spline

    call write_attr(hg, 'n', this%n)

    call write_dset(hg, 'x', this%x)
    call write_dset(hg, 'y', this%y)
    call write_dset(hg, 'dy_dx', this%dy_dx)

    ! Finish

    return

  end subroutine write

  $endif

!****

  $if($MPI)

  subroutine bcast_sp (sp, root_rank)

    class(spline_t), intent(inout) :: sp
    integer, intent(in)            :: root_rank

    ! Broadcast the spline

    call bcast(sp%n, root_rank)

    call bcast_alloc(sp%x, root_rank)
    call bcast_alloc(sp%y, root_rank)
    call bcast_alloc(sp%dy_dx, root_rank)

    ! Finish

    return

  end subroutine bcast_sp

  $endif

!****

  function interp_1 (this, x) result (y)

    class(spline_t), intent(in) :: this
    real(WP), intent(in)        :: x
    real(WP)                    :: y

    real(WP) :: y_n(1)

    ! Interpolate y at a single point

    y_n = this%interp([x])

    y = y_n(1)

    ! Finish

    return

  end function interp_1

!****

  function interp_v (this, x) result (y)

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

       y(j) = this%y(i  )*phi(1._WP-w) + &
              this%y(i+1)*phi(w      ) - &
        h*this%dy_dx(i  )*psi(1._WP-w) + &
        h*this%dy_dx(i+1)*psi(w      )

    end do x_loop

    ! Finish

    return

  contains

    function phi (t)

      real(WP), intent(in) :: t
      real(WP)             :: phi

      phi = 3._WP*t**2 - 2._WP*t**3

      return

    end function phi

    function psi (t)

      real(WP), intent(in) :: t
      real(WP)             :: psi

      psi = t**3 - t**2

      return

    end function psi

  end function interp_v

!****

  function interp_n (this) result (y)

    class(spline_t), intent(in) :: this
    real(WP)                    :: y(this%n)

    ! Interpolate y at abscissa points

    y = this%y

    ! Finish

  end function interp_n

!****

  function deriv_1 (this, x) result (dy_dx)

    class(spline_t), intent(in) :: this
    real(WP), intent(in)        :: x
    real(WP)                    :: dy_dx

    real(WP) :: dy_dx_n(1)

    ! Differentiate y at a single point

    dy_dx_n = this%deriv([x])

    dy_dx = dy_dx_n(1)

    ! Finish

    return

  end function deriv_1

!****

  function deriv_v (this, x) result (dy_dx)

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

       dy_dx(j) = -this%y(i  )*dphi_dt(1._WP-w)/h + &
                   this%y(i+1)*dphi_dt(w      )/h + &
               this%dy_dx(i  )*dpsi_dt(1._WP-w) + &
               this%dy_dx(i+1)*dpsi_dt(w      )

    end do x_loop

    ! Finish

    return

  contains

    function dphi_dt (t)

      real(WP), intent(in) :: t
      real(WP)             :: dphi_dt

      dphi_dt = 6._WP*t - 6._WP*t**2

      return

    end function dphi_dt

    function dpsi_dt (t)

      real(WP), intent(in) :: t
      real(WP)             :: dpsi_dt

      dpsi_dt = 3._WP*t**2 - 2._WP*t

      return

    end function dpsi_dt

  end function deriv_v

!****

  function deriv_n (this) result (dy_dx)

    class(spline_t), intent(in) :: this
    real(WP)                    :: dy_dx(this%n)

    ! Differentiate y at abscissa points

    dy_dx = this%dy_dx

    ! Finish

    return

  end function deriv_n

!****

  function integ_n (this) result (Y)

    class(spline_t), intent(in) :: this
    real(WP)                    :: Y(this%n)

    integer  :: i
    real(WP) :: h

    ! Calculate the integral of y

    Y(1) = 0._WP

    x_loop : do i = 1,this%n-1
       
       h = this%x(i+1) - this%x(i)

       Y(i+1) = Y(i) + (this%y(i) + this%y(i+1))*h/2._WP + &
                       (this%dy_dx(i+1) - this%dy_dx(i))*h**2/12._WP

    end do x_loop

    ! Finish

    return

  end function integ_n

!****

  function natural_dy_dx (x, y, dy_dx_a, dy_dx_b) result (dy_dx)

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

  end function natural_dy_dx

!****

  function findiff_dy_dx (x, y, dy_dx_a, dy_dx_b) result (dy_dx)

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

  end function findiff_dy_dx

!****

  function mono_dy_dx (x, y, dy_dx_a, dy_dx_b) result (dy_dx)

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

  end function mono_dy_dx

end module core_spline
