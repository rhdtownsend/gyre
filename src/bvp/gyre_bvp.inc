! Incfile  : gyre_bvp
! Purpose  : parametric boundary value problems (template)
!
! Copyright 2013-2017 Rich Townsend
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_bvp

  ! Uses

  use core_kinds

  use gyre_bound
  use gyre_diff
  use gyre_ext
  use gyre_num_par
  use gyre_state
  use gyre_status
  use gyre_sysmtx
  use gyre_sysmtx_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ${T}_bvp_t
     private
     class(${T}_diff_t), allocatable   :: df(:)
     class(${T}_bound_t), allocatable  :: bd
     class(${T}_sysmtx_t), allocatable :: sm
     integer, public                   :: n_k
     integer, public                   :: n_e
     integer, public                   :: n_i
     integer, public                   :: n_o
     logical                           :: factored
   contains 
     private
     procedure, public :: build
     procedure, public :: factor
     procedure, public :: det
     procedure, public :: soln_vec_hom
     procedure, public :: soln_vec_inhom
     procedure, public :: resd_vec
  end type ${T}_bvp_t

  ! Interfaces

  interface ${T}_bvp_t
     module procedure ${T}_bvp_t
  end interface ${T}_bvp_t

  ! Access specifiers

  private

  public :: ${T}_bvp_t
  public :: ${T}_bvp_t_

contains

  function ${T}_bvp_t_ (bd, df, nm_p) result (bp)

    class(${T}_diff_t), intent(in)   :: df(:)
    class(${T}_bound_t), intent(in)  :: bd
    type(num_par_t), intent(in)      :: nm_p
    type(${T}_bvp_t)                 :: bp

    integer :: n_k
    integer :: n_e
    integer :: n_i
    integer :: n_o

    ! Perform basic validations

    n_k = SIZE(df) + 1
    n_e = bd%n_e

    n_i = bd%n_i
    n_o = bd%n_o

    $CHECK_BOUNDS(n_i+n_o,n_e)

    $ASSERT_DEBUG(ALL(df%n_e == n_e),Dimension mismatch)

    ! Construct the bvp_t

    allocate(bp%df(n_k-1), SOURCE=df)

    allocate(bp%bd, SOURCE=bd)

    allocate(bp%sm, SOURCE=${T}_sysmtx_t(n_k-1, n_e, n_i, n_o, nm_p))

    bp%n_k = n_k
    bp%n_e = n_e
    bp%n_i = n_i
    bp%n_o = n_o

    bp%factored = .FALSE.

    ! Finish
    
    return

  end function ${T}_bvp_t_

  !****

  subroutine build (this, st)

    class(${T}_bvp_t), target, intent(inout) :: this
    class(${T}_state_t), intent(in)          :: st

    $TYPE(WP)        :: B_i(this%n_i,this%n_e)
    $TYPE(WP)        :: B_o(this%n_o,this%n_e)
    $TYPE(WP)        :: E_l(this%n_e,this%n_e)
    $TYPE(WP)        :: E_r(this%n_e,this%n_e)
    $TYPE(WP)        :: scl_i(this%n_i)
    $TYPE(WP)        :: scl_o(this%n_o)
    type(${T}_ext_t) :: scl
    integer          :: k
    
    ! Build the bvp for the specified state

    ! Set up boundary conditions

    call this%bd%build_i(st, B_i, scl_i)
    call this%sm%set_B_i(B_i, scl_i)

    call this%bd%build_o(st, B_o, scl_o)
    call this%sm%set_B_o(B_o, scl_o)

    ! Set up difference equations

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, scl) SCHEDULE (DYNAMIC)
    sub_loop : do k = 1, this%n_k-1
       call this%df(k)%build(st, E_l, E_r, scl)
       call this%sm%set_E(k, E_l, E_r, scl)
    end do sub_loop

    ! Reset the factored flag

    this%factored = .FALSE.

    ! Finish

    return

  end subroutine build

  !****

  subroutine factor (this)

    class(${T}_bvp_t), intent(inout) :: this

    ! Factorize the sysmtx

    call this%sm%factor()

    this%factored = .TRUE.

    ! Finish

    return

  end subroutine factor

  !****

  function det (this)

    class(${T}_bvp_t), intent(inout) :: this
    type(${T}_ext_t)                 :: det

    $ASSERT(this%factored,Matrix has not been factorized)

    ! Evaluate the determinant of the sysmtx

    det = this%sm%det()

    ! Finish

    return

  end function det

  !****

  function soln_vec_hom (this) result (y)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP)                        :: y(this%n_e,this%n_k)

    $TYPE(WP) :: v(this%n_e*this%n_k)

    $ASSERT_DEBUG(this%factored,Matrix has not been factorized)

    ! Evaluate the solution vector y of the homogeneous system

    v = this%sm%soln_vec_hom()

    y = RESHAPE(v, SHAPE(y))

    ! Finish

    return

  end function soln_vec_hom

  !****

  function soln_vec_inhom (this, w_i, w_o) result (y)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: w_i(:)
    $TYPE(WP), intent(in)            :: w_o(:)
    $TYPE(WP)                        :: y(this%n_e,this%n_k)

    $TYPE(WP) :: v(this%n_e*this%n_k)

    $CHECK_BOUNDS(SIZE(w_i),this%n_i)
    $CHECK_BOUNDS(SIZE(w_o),this%n_o)

    $ASSERT_DEBUG(this%factored,Matrix has not been factorized)

    ! Evaluate the solution vector y of the inhomogeneous system

    v = this%sm%soln_vec_inhom(w_i, w_o)

    y = RESHAPE(v, SHAPE(y))
    
    ! Finish

    return

  end function soln_vec_inhom

  !****

  function resd_vec (this, y) result (w)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: y(:,:)
    $TYPE(WP)                        :: w(this%n_e*this%n_k)

    $TYPE(WP) :: v(this%n_e*this%n_k)

    $CHECK_BOUNDS(SIZE(y,1),this%n_e)
    $CHECK_BOUNDS(SIZE(y,2),this%n_k)

    $ASSERT_DEBUG(.NOT. this%factored,Matrix has already been factorized)

    ! Evaluate the residuals vector

    v = RESHAPE(y, SHAPE(v))

    w = this%sm%resd_vec(v)

    ! Finish

    return

  end function resd_vec

end module gyre_${T}_bvp

$endsub
