! Incfile : bvp
! Purpose : parametric boundary value problems (template)
!
! Copyright 2013-2022 Rich Townsend
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

module ${T}_bvp_m

  ! Uses

  use kinds_m

  use bound_m
  use diff_m
  use ext_m
  use num_par_m
  use state_m
  use status_m
  use sysmtx_m
  use sysmtx_factory_m

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ${T}_bvp_t
     private
     class(${T}_diff_t), allocatable   :: df(:)
     class(${T}_bound_t), allocatable  :: bd
     class(${T}_sysmtx_t), allocatable :: sm
     integer, public                   :: n
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
     procedure, public :: resd_vec_hom
     procedure, public :: resd_vec_inhom
  end type ${T}_bvp_t

  ! Interfaces

  interface ${T}_bvp_t
     module procedure ${T}_bvp_t_
  end interface ${T}_bvp_t

  ! Access specifiers

  private

  public :: ${T}_bvp_t

contains

  function ${T}_bvp_t_ (bd, df, nm_p) result (bp)

    class(${T}_diff_t), intent(in)   :: df(:)
    class(${T}_bound_t), intent(in)  :: bd
    type(num_par_t), intent(in)      :: nm_p
    type(${T}_bvp_t)                 :: bp

    integer :: n
    integer :: n_e
    integer :: n_i
    integer :: n_o

    ! Perform basic validations

    n = SIZE(df) + 1
    n_e = bd%n_e

    n_i = bd%n_i
    n_o = bd%n_o

    $CHECK_BOUNDS(n_i+n_o,n_e)

    $ASSERT_DEBUG(ALL(df%n_e == n_e),Dimension mismatch)

    ! Construct the bvp_t

    allocate(bp%df(n-1), SOURCE=df)

    allocate(bp%bd, SOURCE=bd)

    allocate(bp%sm, SOURCE=${T}_sysmtx_t(n-1, n_e, n_i, n_o, nm_p))

    bp%n = n
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
    integer          :: j
    
    ! Build the bvp for the specified state

    ! Set up boundary conditions

    call this%bd%build_i(st, B_i, scl_i)
    call this%sm%set_B_i(B_i, scl_i)

    call this%bd%build_o(st, B_o, scl_o)
    call this%sm%set_B_o(B_o, scl_o)

    ! Set up difference equations

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, scl) SCHEDULE (DYNAMIC)
    sub_loop : do j = 1, this%n-1
       call this%df(j)%build(st, E_l, E_r, scl)
       call this%sm%set_E(j, E_l, E_r, scl)
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
    $TYPE(WP)                        :: y(this%n_e,this%n)

    $TYPE(WP) :: u(this%n_e*this%n)

    $ASSERT_DEBUG(this%factored,Matrix has not been factorized)

    ! Evaluate the solution vector y of the homogeneous system

    u = this%sm%soln_vec_hom()

    y = RESHAPE(u, SHAPE(y))

    ! Finish

    return

  end function soln_vec_hom

  !****

  function soln_vec_inhom (this, z_i, z_o) result (y)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: z_i(:)
    $TYPE(WP), intent(in)            :: z_o(:)
    $TYPE(WP)                        :: y(this%n_e,this%n)

    $TYPE(WP) :: u(this%n_e*this%n)

    $CHECK_BOUNDS(SIZE(z_i),this%n_i)
    $CHECK_BOUNDS(SIZE(z_o),this%n_o)

    $ASSERT_DEBUG(this%factored,Matrix has not been factorized)

    ! Evaluate the solution vector y of the inhomogeneous system

    u = this%sm%soln_vec_inhom(z_i, z_o)

    y = RESHAPE(u, SHAPE(y))
    
    ! Finish

    return

  end function soln_vec_inhom

  !****

  function resd_vec_hom (this, y) result (dz)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: y(:,:)
    $TYPE(WP)                        :: dz(this%n_e,this%n)

    $TYPE(WP) :: u(this%n_e*this%n)
    $TYPE(WP) :: dv(this%n_e*this%n)

    $CHECK_BOUNDS(SIZE(y,1),this%n_e)
    $CHECK_BOUNDS(SIZE(y,2),this%n)

    $ASSERT_DEBUG(.NOT. this%factored,Matrix has already been factorized)

    ! Evaluate the residuals vector dy of the homogeneous system

    u = RESHAPE(y, SHAPE(u))

    dv = this%sm%resd_vec_hom(u)

    dz = RESHAPE(dv, SHAPE(dz))

    ! Finish

    return

  end function resd_vec_hom

  !****

  function resd_vec_inhom (this, y, z_i, z_o) result (dz)

    class(${T}_bvp_t), intent(inout) :: this
    $TYPE(WP), intent(in)            :: y(:,:)
    $TYPE(WP), intent(in)            :: z_i(:)
    $TYPE(WP), intent(in)            :: z_o(:)
    $TYPE(WP)                        :: dz(this%n_e,this%n)

    $TYPE(WP) :: u(this%n_e*this%n)
    $TYPE(WP) :: dv(this%n_e*this%n)

    $CHECK_BOUNDS(SIZE(y,1),this%n_e)
    $CHECK_BOUNDS(SIZE(y,2),this%n)

    $CHECK_BOUNDS(SIZE(z_i),this%n_i)
    $CHECK_BOUNDS(SIZE(z_o),this%n_o)

    $ASSERT_DEBUG(.NOT. this%factored,Matrix has already been factorized)

    ! Evaluate the residuals vector dy of the inhomogeneous system

    u = RESHAPE(y, SHAPE(u))

    dv = this%sm%resd_vec_inhom(u, z_i, z_o)

    dz = RESHAPE(dv, SHAPE(dz))

    ! Finish

    return

  end function resd_vec_inhom

end module ${T}_bvp_m

$endsub
