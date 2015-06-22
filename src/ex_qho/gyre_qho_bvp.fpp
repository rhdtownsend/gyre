! Module   : gyre_qho_bvp
! Purpose  : boundary-value solver (adiabatic)
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

module gyre_qho_bvp

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_ivp
  use gyre_ivp_factory
  use gyre_num_par
  use gyre_sysmtx
  use gyre_sysmtx_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_bvp_t) :: qho_bvp_t
  end type qho_bvp_t

  ! Interfaces

  interface qho_bvp_t
     module procedure qho_bvp_t_
  end interface qho_bvp_t

  ! Access specifiers

  private

  public :: qho_bvp_t

  ! Procedures

contains

  function qho_bvp_t_ (x, c, np, omega_min, omega_max) result (bp)

    use gyre_qho_eqns
    use gyre_qho_bound

    real(WP), intent(in)        :: x(:)
    real(WP), intent(in)        :: c
    type(num_par_t), intent(in) :: np
    real(WP), intent(in)        :: omega_min
    real(WP), intent(in)        :: omega_max
    type(qho_bvp_t), target     :: bp

    type(qho_eqns_t)               :: eq
    integer                        :: n
    real(WP)                       :: x_i
    real(WP)                       :: x_o
    type(qho_bound_t)              :: bd
    class(r_ivp_t), allocatable    :: iv
    class(r_sysmtx_t), allocatable :: sm

    ! Construct the qho_bvp_t

    ! Initialize the equations

    eq = qho_eqns_t(c)

    ! Initialize the boundary conditions

    n = SIZE(x)

    x_i = x(1)
    x_o = x(n)

    bd = qho_bound_t(x_i, x_o)

    ! Initialize the IVP solver

    allocate(iv, SOURCE=r_ivp_t(eq, np))

    ! Initialize the system matrix

    allocate(sm, SOURCE=r_sysmtx_t(n-1, eq%n_e, bd%n_i, bd%n_o, np))

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t(x, eq, bd, iv, sm, omega_min, omega_max)

    ! Finish

    return

  end function qho_bvp_t_

end module gyre_qho_bvp

