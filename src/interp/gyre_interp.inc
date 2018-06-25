! Module   : gyre_interp
! Purpose  : piecewise cubic interpolation (template)
!
! Copyright 2015-2018 Rich Townsend
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
$include 'core_parallel.inc'

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_interp

  ! Uses

  use core_kinds
  $if ($HDF5)
  use core_hgroup
  $endif
  use core_linalg
  use core_parallel
  use core_order

  use gyre_deriv

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type ${T}_interp_t
     private
     real(WP), allocatable  :: x_(:)     ! Abscissa
     $TYPE(WP), allocatable :: f_(:)     ! Ordinate
     $TYPE(WP), allocatable :: df_dx_(:) ! First derivatives
     $TYPE(WP), allocatable :: int_f_(:) ! First integrals
     integer                :: n         ! Number of points
   contains
     private
     procedure         :: x_n_
     generic, public   :: x => x_n_
     procedure, public :: x_min => x_min_
     procedure, public :: x_max => x_max_
     procedure         :: f_1_
     procedure         :: f_v_
     procedure         :: f_n_
     generic, public   :: f => f_1_, f_v_, f_n_
     procedure         :: df_dx_1_
     procedure         :: df_dx_v_
     procedure         :: df_dx_n_
     generic, public   :: df_dx => df_dx_1_, df_dx_v_, df_dx_n_
     procedure         :: int_f_1_
     procedure         :: int_f_v_
     procedure         :: int_f_n_
     generic, public   :: int_f => int_f_1_, int_f_v_, int_f_n_
  end type ${T}_interp_t

  ! Interfaces

  interface ${T}_interp_t
     module procedure ${T}_interp_t_
     module procedure ${T}_interp_t_eval_derivs_
  end interface ${T}_interp_t

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

  public :: ${T}_interp_t
  public :: ${T}_interp_t_eval_derivs_
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

  function ${T}_interp_t_ (x, f, df_dx) result (in)

    real(WP), intent(in)  :: x(:)
    $TYPE(WP), intent(in) :: f(:)
    $TYPE(WP), intent(in) :: df_dx(:)
    type(${T}_interp_t)   :: in

    integer  :: n
    integer  :: i
    real(WP) :: h
    
    $CHECK_BOUNDS(SIZE(f),SIZE(x))
    $CHECK_BOUNDS(SIZE(df_dx),SIZE(x))

    $ASSERT(ALL(x(2:) > x(:SIZE(x)-1)),Non-monotonic abscissa)

    ! Construct the interp_t

    ! Store values

    in%x_ = x
    in%f_ = f

    in%df_dx_ = df_dx

    in%n = SIZE(x)

    ! Evaluate the integral

    allocate(in%int_f_(in%n))

    in%int_f_(1) = 0._WP

    x_loop : do i = 1, in%n-1

       h = x(i+1) - x(i)

       in%int_f_(i+1) = in%int_f_(i) - f(i  )*(dphi_dt_(0._WP, -1) - dphi_dt_(1._WP, -1))*h +    &
                                       f(i+1)*(dphi_dt_(1._WP, -1) - dphi_dt_(0._WP, -1))*h +    &
                                   df_dx(i  )*(dpsi_dt_(0._WP, -1) - dpsi_dt_(1._WP, -1))*h**2 + &
                                   df_dx(i+1)*(dpsi_dt_(1._WP, -1) - dpsi_dt_(0._WP, -1))*h**2

    end do x_loop

    ! Finish

    return

  end function ${T}_interp_t_
  
  !****

  function ${T}_interp_t_eval_derivs_ (x, f, deriv_type, df_dx_a, df_dx_b) result (in)

    real(WP), intent(in)            :: x(:)
    $TYPE(WP), intent(in)           :: f(:)
    character(*), intent(in)        :: deriv_type
    $TYPE(WP), optional, intent(in) :: df_dx_a
    $TYPE(WP), optional, intent(in) :: df_dx_b
    type(${T}_interp_t)             :: in

    $CHECK_BOUNDS(SIZE(f),SIZE(x))

    ! Construct the interp_t, with derivatives calculated according to
    ! deriv_type

    in = ${T}_interp_t(x, f, deriv(x, f, deriv_type, df_dx_a, df_dx_b))

    ! Finish

    return

  end function ${T}_interp_t_eval_derivs_

  !****

  $if ($HDF5)

  subroutine read_ (hg, in)

    type(hgroup_t), intent(inout)    :: hg
    type(${T}_interp_t), intent(out) :: in

    real(WP), allocatable  :: x(:)
    $TYPE(WP), allocatable :: f(:)
    $TYPE(WP), allocatable :: df_dx(:)

    ! Read the interp_t

    call read_dset_alloc(hg, 'x', x)
    call read_dset_alloc(hg, 'f', f)
    call read_dset_alloc(hg, 'df_dx', df_dx)

    in = ${T}_interp_t(x, f, df_dx)

    ! Finish

    return

  end subroutine read_

  !****

  subroutine write_ (hg, in)

    type(hgroup_t), intent(inout)   :: hg
    type(${T}_interp_t), intent(in) :: in

    ! Write the interp_t

    call write_attr(hg, 'n', in%n)

    call write_dset(hg, 'x', in%x_)
    call write_dset(hg, 'f', in%f_)
    call write_dset(hg, 'df_dx', in%df_dx_)

    ! Finish

    return

  end subroutine write_

  $endif

  !****

  $if ($MPI)

  subroutine bcast_0_ (in, root_rank)

    class(${T}_interp_t), intent(inout) :: in
    integer, intent(in)                 :: root_rank

    ! Broadcast the interp

    call bcast(in%n, root_rank)

    call bcast_alloc(in%x_, root_rank)
    call bcast_alloc(in%f_, root_rank)
    call bcast_alloc(in%df_dx_, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST(type(${T}_interp_t),1)
  $BCAST(type(${T}_interp_t),2)
  $BCAST(type(${T}_interp_t),3)
  $BCAST(type(${T}_interp_t),4)

  $BCAST_ALLOC(type(${T}_interp_t),0)
  $BCAST_ALLOC(type(${T}_interp_t),1)
  $BCAST_ALLOC(type(${T}_interp_t),2)
  $BCAST_ALLOC(type(${T}_interp_t),3)
  $BCAST_ALLOC(type(${T}_interp_t),4)

  $endif

  !****

  function x_n_ (this) result (x)

    class(${T}_interp_t), intent(in) :: this
    real(WP)                         :: x(this%n)

    ! Return the abscissa points

    x = this%x_

    return

  end function x_n_

  !****

  function x_min_ (this) result (x_min)

    class(${T}_interp_t), intent(in) :: this
    real(WP)                         :: x_min

    ! Return the minimum abscissa point

    x_min = this%x_(1)

    ! FInish

    return

  end function x_min_

  !****

  function x_max_ (this) result (x_max)

    class(${T}_interp_t), intent(in) :: this
    real(WP)                         :: x_max

    ! Return the maximum abscissa point

    x_max = this%x_(this%n)

    ! FInish

    return

  end function x_max_

  !****

  function f_1_ (this, x) result (f)

    class(${T}_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x
    $TYPE(WP)                        :: f

    ! Interpolate f at a single point

    f = this%df_dx(x, 0)

    ! Finish

    return

  end function f_1_

  !****

  function f_v_ (this, x) result (f)

    class(${T}_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x(:)
    $TYPE(WP)                        :: f(SIZE(x))

    ! Interpolate f at a vector of points

    f = this%df_dx(x, 0)

    ! Finish

    return

  end function f_v_

  !****

  function f_n_ (this) result (f)

    class(${T}_interp_t), intent(in) :: this
    $TYPE(WP)                        :: f(this%n)

    ! Interpolate f at abscissa points

    f = this%df_dx(0)

    ! Finish

  end function f_n_

  !****

  function df_dx_1_ (this, x, n) result (df_dx)

    class(${T}_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x
    integer, optional, intent(in)    :: n
    $TYPE(WP)                        :: df_dx

    integer  :: n_
    integer  :: i
    real(WP) :: h
    real(WP) :: w

    if (PRESENT(n)) then
       n_ = n
    else
       n_ = 1
    endif

    ! Interpolate the first derivative of f, or the n'th derivative if
    ! n is present, at a single point

    ! Set up the bracketing index

    i = 1

    call locate(this%x_, x, i)
    $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

    ! Set up the interpolation weights

    h = this%x_(i+1) - this%x_(i)
    w = (x - this%x_(i))/h

    ! Do the interpolation

    df_dx = (-1._WP)**n_ *     this%f_(i  )*dphi_dt_(1._WP-w, n_)/h**n_ +     &
                               this%f_(i+1)*dphi_dt_(      w, n_)/h**n_ -     &
            (-1._WP)**n_ * this%df_dx_(i  )*dpsi_dt_(1._WP-w, n_)/h**(n_-1) + &
                           this%df_dx_(i+1)*dpsi_dt_(      w, n_)/h**(n_-1)

    ! Finish

    return

  end function df_dx_1_

  !****

  function df_dx_v_ (this, x, n) result (df_dx)

    class(${T}_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x(:)
    integer, optional, intent(in)    :: n
    $TYPE(WP)                        :: df_dx(SIZE(x))

    integer  :: n_
    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: w

    if (PRESENT(n)) then
       n_ = n
    else
       n_ = 1
    endif

    ! Interpolate the first derivative of f, or the n'th derivative if
    ! n is present, at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%x_, x(j), i)
       $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

       ! Set up the interpolation weights

       h = this%x_(i+1) - this%x_(i)
       w = (x(j) - this%x_(i))/h

       ! Do the interpolation

       df_dx(j) = (-1._WP)**n_ *     this%f_(i  )*dphi_dt_(1._WP-w, n_)/h**n_ +     &
                                     this%f_(i+1)*dphi_dt_(w      , n_)/h**n_ -     &
                  (-1._WP)**n_ * this%df_dx_(i  )*dpsi_dt_(1._WP-w, n_)/h**(n_-1) + &
                                 this%df_dx_(i+1)*dpsi_dt_(w      , n_)/h**(n_-1)

    end do x_loop

    ! Finish

    return

  end function df_dx_v_

  !****

  function df_dx_n_ (this, n) result (df_dx)

    class(${T}_interp_t), intent(in) :: this
    integer, optional, intent(in)    :: n
    $TYPE(WP)                        :: df_dx(this%n)

    integer :: n_

    if (PRESENT(n)) then
       n_ = n
    else
       n_ = 1
    endif

    ! Interpolate the first derivative of f, or the n'th derivative if
    ! n is present, at abscissa points

    select case (n_)
    case (0)
       df_dx = this%f_
    case (1)
       df_dx = this%df_dx_
    case default
       df_dx = this%df_dx(this%x_, n_)
    end select
       
    ! Finish

    return

  end function df_dx_n_

  !****

  function int_f_1_ (this, x) result (int_f)

    class(${T}_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x
    $TYPE(WP)                        :: int_f

    integer  :: i
    real(WP) :: h
    real(WP) :: w
    
    ! Interpolate the integral of f at a single point

    ! Set up the bracketing index

    i = 1

    call locate(this%x_, x, i)
    $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

    ! Set up the interpolation weights

    h = this%x_(i+1) - this%x_(i)
    w = (x - this%x_(i))/h

    ! Do the interpolation

    int_f = this%int_f_(i) - this%f_(i  )*(dphi_dt_(1._WP-w, -1) - dphi_dt_(1._WP, -1))*h +    &
                             this%f_(i+1)*(dphi_dt_(      w, -1) - dphi_dt_(0._WP, -1))*h +    &
                         this%df_dx_(i  )*(dpsi_dt_(1._WP-w, -1) - dpsi_dt_(1._WP, -1))*h**2 + &
                         this%df_dx_(i+1)*(dpsi_dt_(      w, -1) - dpsi_dt_(0._WP, -1))*h**2

    ! Finish

    return

  end function int_f_1_
    
  !****

  function int_f_v_ (this, x) result (int_f)

    class(${T}_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x(:)
    $TYPE(WP)                        :: int_f(SIZE(x))

    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: w
    
    ! Interpolate the integral of f at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%x_, x(j), i)
       $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

       ! Set up the interpolation weights

       h = this%x_(i+1) - this%x_(i)
       w = (x(j) - this%x_(i))/h

       ! Do the interpolation

       int_f(j) = this%int_f_(i) - this%f_(i  )*(dphi_dt_(1._WP-w, -1) - dphi_dt_(1._WP, -1))*h +    &
                                   this%f_(i+1)*(dphi_dt_(      w, -1) - dphi_dt_(0._WP, -1))*h +    &
                               this%df_dx_(i  )*(dpsi_dt_(1._WP-w, -1) - dpsi_dt_(1._WP, -1))*h**2 + &
                               this%df_dx_(i+1)*(dpsi_dt_(      w, -1) - dpsi_dt_(0._WP, -1))*h**2

    end do x_loop

    ! Finish

    return

  end function int_f_v_

  !****

  function int_f_n_ (this) result (int_f)

    class(${T}_interp_t), intent(in) :: this
    $TYPE(WP)                        :: int_f(this%n)

    ! Interpolate the integral of f at abscissa points

    int_f = this%df_dx_
       
    ! Finish

    return

  end function int_f_n_

  !****

  function dphi_dt_ (t, n) result (dphi_dt)

    real(WP), intent(in) :: t
    integer, intent(in)  :: n
    real(WP)             :: dphi_dt

    ! Evaluate the n'th derivative of the phi basis function

    select case (n)
    case (-1)
       dphi_dt = t**3 - t**4/2._WP
    case (0)
       dphi_dt = 3._WP*t**2 - 2._WP*t**3
    case (1)
       dphi_dt = 6._WP*t - 6._WP*t**2
    case (2)
       dphi_dt = 6._WP - 12._WP*t
    case (3)
       dphi_dt = -12._WP
    case (4:)
       dphi_dt = 0._WP
    case default
       $ABORT(Invalid n)
    end select

    ! Finish

    return

  end function dphi_dt_

  !****

  function dpsi_dt_ (t, n) result (dpsi_dt)

    real(WP), intent(in) :: t
    integer, intent(in)  :: n
    real(WP)             :: dpsi_dt

    ! Evaluate the n'th derivative of the psi basis function

    select case (n)
    case (-1)
       dpsi_dt = t**4/4._WP - t**3/3._WP
    case (0)
       dpsi_dt = t**3 - t**2
    case (1)
       dpsi_dt = 3._WP*t**2 - 2._WP*t
    case (2)
       dpsi_dt = 6._WP*t - 2._WP
    case (3)
       dpsi_dt = 6._WP
    case (4:)
       dpsi_dt = 0._WP
    case default
       $ABORT(Invalid n)
    end select

    ! Finish

    return

  end function dpsi_dt_

end module gyre_${T}_interp

$endsub

