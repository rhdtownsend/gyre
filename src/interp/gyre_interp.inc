! Module   : gyre_interp
! Purpose  : piecewise cubic interpolation (template)
!
! Copyright 2015-2017 Rich Townsend
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

  use gyre_calc

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type ${T}_interp_t
     private
     real(WP), allocatable  :: x_(:)     ! Abscissa
     $TYPE(WP), allocatable :: f_(:)     ! Ordinate
     $TYPE(WP), allocatable :: df_dx_(:) ! First derivatives
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
     procedure         :: int_f_n_
     generic, public   :: int_f => int_f_n_
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
    
    $CHECK_BOUNDS(SIZE(f),SIZE(x))
    $CHECK_BOUNDS(SIZE(df_dx),SIZE(x))

    $ASSERT(ALL(x(2:) > x(:SIZE(x)-1)),Non-monotonic abscissa)

    ! Construct the interp_t

    in%x_ = x
    in%f_ = f

    in%df_dx_ = df_dx

    in%n = SIZE(x)

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

    integer  :: i
    real(WP) :: h
    real(WP) :: w

    ! Interpolate f at a single point

    ! Set up the bracketing index

    i = 1

    call locate(this%x_, x, i)
    $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

    ! Set up the interpolation weights

    h = this%x_(i+1) - this%x_(i)
    w = (x - this%x_(i))/h

    ! Do the interpolation

    f =    this%f_(i  )*phi_(1._WP-w) + &
           this%f_(i+1)*phi_(w      ) - &
     h*this%df_dx_(i  )*psi_(1._WP-w) + &
     h*this%df_dx_(i+1)*psi_(w      )

    ! Finish

    return

  end function f_1_

  !****

  function f_v_ (this, x) result (f)

    class(${T}_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x(:)
    $TYPE(WP)                        :: f(SIZE(x))

    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: w

    ! Interpolate f at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%x_, x(j), i)
       $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

       ! Set up the interpolation weights

       h = this%x_(i+1) - this%x_(i)
       w = (x(j) - this%x_(i))/h

       ! Do the interpolation

       f(j) = this%f_(i  )*phi_(1._WP-w) + &
              this%f_(i+1)*phi_(w      ) - &
        h*this%df_dx_(i  )*psi_(1._WP-w) + &
        h*this%df_dx_(i+1)*psi_(w      )

    end do x_loop

    ! Finish

    return

  end function f_v_

  !****

  function f_n_ (this) result (f)

    class(${T}_interp_t), intent(in) :: this
    $TYPE(WP)                        :: f(this%n)

    ! Return f at abscissa points

    f = this%f_

    ! Finish

  end function f_n_

  !****

  function df_dx_1_ (this, x) result (df_dx)

    class(${T}_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x
    $TYPE(WP)                        :: df_dx

    integer  :: i
    real(WP) :: h
    real(WP) :: w

    ! Differentiate f at a single point

    ! Set up the bracketing index

    i = 1

    call locate(this%x_, x, i)
    $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

    ! Set up the interpolation weights

    h = this%x_(i+1) - this%x_(i)
    w = (x - this%x_(i))/h

    ! Do the interpolation

    df_dx =    -this%f_(i  )*dphi_dt_(1._WP-w)/h + &
                this%f_(i+1)*dphi_dt_(w      )/h + &
            this%df_dx_(i  )*dpsi_dt_(1._WP-w) + &
            this%df_dx_(i+1)*dpsi_dt_(w      )

    ! Finish

    return

  end function df_dx_1_

  !****

  function df_dx_v_ (this, x) result (df_dx)

    class(${T}_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x(:)
    $TYPE(WP)                        :: df_dx(SIZE(x))

    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: w

    ! Differentiate f at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%x_, x(j), i)
       $ASSERT(i > 0 .AND. i < this%n,Out-of-bounds interpolation)

       ! Set up the interpolation weights

       h = this%x_(i+1) - this%x_(i)
       w = (x(j) - this%x_(i))/h

       ! Do the interpolation

       df_dx(j) = -this%f_(i  )*dphi_dt_(1._WP-w)/h + &
                   this%f_(i+1)*dphi_dt_(w      )/h + &
               this%df_dx_(i  )*dpsi_dt_(1._WP-w) + &
               this%df_dx_(i+1)*dpsi_dt_(w      )

    end do x_loop

    ! Finish

    return

  end function df_dx_v_

  !****

  function df_dx_n_ (this) result (df_dx)

    class(${T}_interp_t), intent(in) :: this
    $TYPE(WP)                   :: df_dx(this%n)

    ! Return df_dx at abscissa points

    df_dx = this%df_dx_

    ! Finish

    return

  end function df_dx_n_

  !****

  function int_f_n_ (this) result (int_f)

    class(${T}_interp_t), intent(in) :: this
    $TYPE(WP)                        :: int_f(this%n)

    integer  :: i
    real(WP) :: h

    ! Integrate f across the full domain

    int_f(1) = 0._WP

    x_loop : do i = 1,this%n-1
       
       h = this%x_(i+1) - this%x_(i)

       int_f(i+1) = int_f(i) + (this%f_(i) + this%f_(i+1))*h/2._WP - &
                               (this%df_dx_(i+1) - this%df_dx_(i))*h**2/12._WP

    end do x_loop

    ! Finish

    return

  end function int_f_n_

  !****

  function phi_ (t)

    real(WP), intent(in) :: t
    real(WP)             :: phi_

    ! Evaluate the phi basis function

    phi_ = 3._WP*t**2 - 2._WP*t**3

    return

  end function phi_

  !****

  function psi_ (t)

    real(WP), intent(in) :: t
    real(WP)             :: psi_

    ! Evaluate the psi basis function 

    psi_ = t**3 - t**2

    return

  end function psi_

  !****

  function dphi_dt_ (t)

    real(WP), intent(in) :: t
    real(WP)             :: dphi_dt_

    ! Evaluate the first derivative of the phi basis function

    dphi_dt_ = 6._WP*t - 6._WP*t**2

    return

  end function dphi_dt_

  !****

  function dpsi_dt_ (t)

    real(WP), intent(in) :: t
    real(WP)             :: dpsi_dt_

    ! Evaluate the first derivative of the psi basis function

    dpsi_dt_ = 3._WP*t**2 - 2._WP*t

    return

  end function dpsi_dt_

end module gyre_${T}_interp

$endsub

