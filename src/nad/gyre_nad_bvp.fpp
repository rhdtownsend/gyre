! Module   : gyre_nad_bvp
! Purpose  : boundary-value solver (nonadiabatic)
!
! Copyright 2013-2015 Rich Townsend
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
  use gyre_ext
  use gyre_ivp
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_nad_eqns
  use gyre_num_par
  use gyre_osc_par
  use gyre_sysmtx
  use gyre_sysmtx_factory
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_bvp_t) :: nad_bvp_t
     class(model_t), pointer :: ml => null()
     type(nad_eqns_t)        :: eq
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

  function nad_bvp_t_ (x, ml, mp, op, np, omega_min, omega_max) result (bp)

    use gyre_nad_magnus_ivp
    use gyre_nad_findiff_ivp
    use gyre_colloc_ivp

    use gyre_nad_eqns
    use gyre_nad_bound

    real(WP), intent(in)                :: x(:)
    class(model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)        :: mp
    type(osc_par_t), intent(in)         :: op
    type(num_par_t), intent(in)         :: np
    real(WP), intent(in)                :: omega_min
    real(WP), intent(in)                :: omega_max
    type(nad_bvp_t), target             :: bp

    class(c_rot_t), allocatable    :: rt
    type(nad_eqns_t)               :: eq
    integer                        :: n
    real(WP)                       :: x_i
    real(WP)                       :: x_o
    type(nad_bound_t)              :: bd
    class(c_ivp_t), allocatable    :: iv
    class(c_sysmtx_t), allocatable :: sm

    ! Construct the nad_bvp_t

    ! Initialize the rotational effects

    allocate(rt, SOURCE=c_rot_t(ml, mp, op))
 
    ! Initialize the equations

    eq = nad_eqns_t(ml, rt, op)

    ! Initialize the boundary conditions

    n = SIZE(x)

    x_i = x(1)
    x_o = x(n)

    bd = nad_bound_t(ml, rt, eq, op, x_i, x_o)

    ! Initialize the IVP solver

    select case (np%ivp_solver)
    case ('MAGNUS_GL2')
       allocate(iv, SOURCE=nad_magnus_ivp_t(ml, eq, 'GL2'))
    case ('MAGNUS_GL4')
       allocate(iv, SOURCE=nad_magnus_ivp_t(ml, eq, 'GL4'))
    case ('MAGNUS_GL6')
       allocate(iv, SOURCE=nad_magnus_ivp_t(ml, eq, 'GL6'))
    case ('COLLOC_GL2')
       allocate(iv, SOURCE=c_colloc_ivp_t(eq, 'GL2'))
    case ('COLLOC_GL4')
       allocate(iv, SOURCE=c_colloc_ivp_t(eq, 'GL4'))
    case ('COLLOC_GL6')
       allocate(iv, SOURCE=c_colloc_ivp_t(eq, 'GL6'))
    case ('FINDIFF')
       allocate(iv, SOURCE=nad_findiff_ivp_t(ml, eq))
   case default
       $ABORT(Invalid ivp_solver)
    end select

    ! Initialize the system matrix

    allocate(sm, SOURCE=c_sysmtx_t(n-1, eq%n_e, bd%n_i, bd%n_o, np))

    ! Initialize the bvp_t

    bp%c_bvp_t = c_bvp_t(x, bd, iv, sm, omega_min, omega_max)

    bp%ml => ml
    bp%eq = eq

    ! Finish

    return

  end function nad_bvp_t_

end module gyre_nad_bvp
