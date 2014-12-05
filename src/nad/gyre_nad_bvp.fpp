! Module   : gyre_nad_bvp
! Purpose  : boundary-value solver (nonadiabatic)
!
! Copyright 2013-2014 Rich Townsend
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

module gyre_nad_bvp

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_model
  use gyre_mode
  use gyre_ext
  use gyre_ivp
  use gyre_sysmtx
  use gyre_modepar
  use gyre_oscpar
  use gyre_numpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_bvp_t) :: nad_bvp_t
  end type nad_bvp_t

  ! Interfaces

  interface nad_bvp_t
     module procedure nad_bvp_t_
  end interface nad_bvp_t

  ! Access specifiers

  private

  public :: nad_bvp_t

  ! Procedures

contains

  function nad_bvp_t_ (x, ml, mp, op, np) result (bp)

    use gyre_nad_jacob
    use gyre_nad_bound

    use gyre_nad_magnus_ivp
    use gyre_nad_findiff_ivp

    use gyre_block_sysmtx

    real(WP), intent(in)                :: x(:)
    class(model_t), pointer, intent(in) :: ml
    type(modepar_t), intent(in)         :: mp
    type(oscpar_t), intent(in)          :: op
    type(numpar_t), intent(in)          :: np
    type(nad_bvp_t), target             :: bp

    type(nad_jacob_t)              :: jc
    integer                        :: n
    real(WP)                       :: x_i
    real(WP)                       :: x_o
    type(nad_bound_t)              :: bd
    class(c_ivp_t), allocatable    :: iv
    class(c_sysmtx_t), allocatable :: sm

    ! Construct the nad_bvp_t

    ! Initialize the jacobian

    select case (op%variables_type)
    case ('DZIEM')
       jc = nad_jacob_t(ml, mp, 'DZIEM')
    case ('JCD')
       jc = nad_jacob_t(ml, mp, 'JCD')
    case default
       $ABORT(Invalid variables_type)
    end select

    ! Initialize the boundary conditions

    n = SIZE(x)

    x_i = x(1)
    x_o = x(n)

    select case (op%outer_bound_type)
    case ('ZERO')
       bd = nad_bound_t(ml, jc, mp, x_i, x_o, 'REGULAR', 'ZERO')
    case ('DZIEM')
       bd = nad_bound_t(ml, jc, mp, x_i, x_o, 'REGULAR', 'DZIEM')
    case ('UNNO')
       bd = nad_bound_t(ml, jc, mp, x_i, x_o, 'REGULAR', 'UNNO')
    case ('JCD')
       bd = nad_bound_t(ml, jc, mp, x_i, x_o, 'REGULAR', 'JCD')
    case default
       $ABORT(Invalid bound_type)
    end select

    ! Initialize the IVP solver

    select case (np%ivp_solver_type)
    case ('MAGNUS_GL2')
       allocate(iv, SOURCE=nad_magnus_ivp_t(ml, jc, 'GL2'))
    case ('FINDIFF')
       allocate(iv, SOURCE=nad_findiff_ivp_t(ml, jc))
    case default
       $ABORT(Invalid ivp_solver_type)
    end select

    ! Initialize the system matrix

    if (np%use_banded) then
       $ABORT(Not yet implemented)
    else
       allocate(sm, SOURCE=c_block_sysmtx_t(n-1, jc%n_e, bd%n_i, bd%n_o))
    endif

    ! Initialize the bvp_t

    bp%c_bvp_t = c_bvp_t(x, ml, jc, bd, iv, sm)

    ! Finish

    return

  end function nad_bvp_t_

end module gyre_nad_bvp
