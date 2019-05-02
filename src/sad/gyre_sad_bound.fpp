! Incfile  : gyre_sad_bound
! Purpose  : static adiabatic boundary conditions
!
! Copyright 2019 Rich Townsend
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

module gyre_sad_bound

  ! Uses

  use core_kinds

  use gyre_atmos
  use gyre_bound
  use gyre_context
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_sad_trans
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: REGULAR_TYPE = 1
  integer, parameter :: VACUUM_TYPE = 1

  integer, parameter :: J_U = 1

  integer, parameter :: J_LAST = J_U
  
  ! Derived-type definitions

  type, extends (r_bound_t) :: sad_bound_t
     private
     type(context_t), pointer :: cx => null()
     type(sad_trans_t)        :: tr
     real(WP), allocatable    :: coeff(:,:)
     integer                  :: type_i
     integer                  :: type_o
     integer                  :: l
   contains 
     private
     procedure         :: stencil_
     procedure, public :: build_i
     procedure         :: build_regular_i_
     procedure, public :: build_o
     procedure         :: build_vacuum_o_
  end type sad_bound_t

  ! Interfaces

  interface sad_bound_t
     module procedure sad_bound_t_
  end interface sad_bound_t

  ! Access specifiers

  private

  public :: sad_bound_t

  ! Procedures

contains

  function sad_bound_t_ (cx, md_p, os_p) result (bd)

    type(context_t), pointer, intent(in) :: cx
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(sad_bound_t)                    :: bd

    class(model_t), pointer :: ml
    type(point_t)           :: pt_i
    type(point_t)           :: pt_o

    ! Construct the sad_bound_t

    bd%cx => cx
    
    bd%tr = sad_trans_t(cx, md_p, os_p)

    ml => cx%model()
    pt_i = cx%point_i()
    pt_o = cx%point_o()

    select case (os_p%inner_bound)
    case ('REGULAR')
       $ASSERT(pt_i%x == 0._WP,Boundary condition invalid for x /= 0)
       bd%type_i = REGULAR_TYPE
    case default
       $ABORT(Invalid inner_bound)
    end select

    select case (os_p%outer_bound)
    case ('VACUUM')
       bd%type_o = VACUUM_TYPE
    case default
       $ABORT(Invalid outer_bound)
    end select

    bd%l = md_p%l

    call bd%stencil_(pt_i, pt_o)

    bd%n_i = 1
    bd%n_o = 1

    bd%n_e = 2

    ! Finish

    return
    
  end function sad_bound_t_

  !****

  subroutine stencil_ (this, pt_i, pt_o)

    class(sad_bound_t), intent(inout) :: this
    type(point_t), intent(in)         :: pt_i
    type(point_t), intent(in)         :: pt_o

    class(model_t), pointer :: ml

    ! Calculate coefficients at the stencil points

    ml => this%cx%model()

    call check_model(ml, [I_U])

    allocate(this%coeff(2,J_LAST))

    ! Inner boundary

    ! Outer boundary

    select case (this%type_o)
    case (VACUUM_TYPE)
       this%coeff(2,J_U) = ml%coeff(I_U, pt_o)
    case default
       $ABORT(Invalid type_o)
    end select

    ! Set up stencil for the tr component

    call this%tr%stencil([pt_i,pt_o])

    ! Finish

    return

  end subroutine stencil_

  !****

  subroutine build_i (this, st, B, scl)

    class(sad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_i)
    
    ! Evaluate the inner boundary conditions

    select case (this%type_i)
    case (REGULAR_TYPE)
       call this%build_regular_i_(st, B, scl)
    case default
       $ABORT(Invalid type_i)
    end select

    ! Apply the variables transformation

    call this%tr%trans_cond(B, 1, st)

    ! Finish

    return

  end subroutine build_i

  !****

  subroutine build_regular_i_ (this, st, B, scl)

    class(sad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    real(WP) :: l

    $CHECK_BOUNDS(SIZE(B, 1),this%n_i)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(scl),this%n_i)

    ! Evaluate the inner boundary conditions (regular-enforcing)

    l = this%l

    ! Set up the boundary conditions

    B(1,1) = l
    B(1,2) = -1._WP

    scl = 1._WP

    ! Finish

    return

  end subroutine build_regular_i_

  !****

  subroutine build_o (this, st, B, scl)

    class(sad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)
    
    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions

    select case (this%type_o)
    case (VACUUM_TYPE)
       call this%build_vacuum_o_(st, B, scl)
    case default
       $ABORT(Invalid type_o)
    end select

    ! Apply the variables transformation

    call this%tr%trans_cond(B, 2, st)

    ! Finish

    return

  end subroutine build_o
  
  !****

  subroutine build_vacuum_o_ (this, st, B, scl)

    class(sad_bound_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: B(:,:)
    real(WP), intent(out)          :: scl(:)

    integer :: l

    $CHECK_BOUNDS(SIZE(B, 1),this%n_o)
    $CHECK_BOUNDS(SIZE(B, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(scl),this%n_o)

    ! Evaluate the outer boundary conditions (vacuum)

    associate ( &
         U => this%coeff(2,J_U))

      l = this%l

      ! Set up the boundary conditions

      B(1,1) = l + 1._WP - U
      B(1,2) = 1._WP
      
      scl = 1._WP

    end associate
    
    ! Finish

    return

  end subroutine build_vacuum_o_

end module gyre_sad_bound
