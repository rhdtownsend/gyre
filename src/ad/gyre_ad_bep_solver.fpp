! Module   : gyre_ad_bep_solver
! Purpose  : adiabatic bounary eigenvalue problem solver
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

module gyre_ad_bep_solver

  ! Uses

  use core_kinds

  use gyre_ad_bound
  use gyre_bep_solver
  use gyre_model
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_bep_solver_t) :: ad_bep_solver_t
     class(model_t), pointer :: ml => null()
     type(mode_par_t)        :: md_p
     type(osc_par_t)         :: os_p
     type(ad_vars_t)         :: vr
     integer, allocatable    :: s(:)
     real(WP), allocatable   :: x(:)
  end type ad_bep_solver_t

  ! Interfaces

  interface ad_bep_solver_t
     module procedure ad_bep_solver_t_
  end interface ad_bep_solver_t

  interface sol_t
     module procedure sol_t_
  end interface sol_t

  ! Access specifiers

  private

  public :: ad_bep_solver_t
  public :: sol_t

  ! Procedures

contains

  function ad_bep_solver_t_ (ml, omega, md_p, nm_p, os_p) result (as)

    class(model_t), pointer, intent(in) :: ml
    real(WP), intent(in)                :: omega(:)
    type(mode_par_t), intent(in)        :: md_p
    type(num_par_t), intent(in)         :: nm_p
    type(osc_par_t), intent(in)         :: os_p
    type(ad_solver_t)                   :: as

    integer, allocatable         :: s(:)
    real(WP), allocatable        :: x(:)
    integer                      :: n_k
    type(ad_bound_t)             :: bd_i
    type(ad_bound_t)             :: bd_o
    type(ad_diff_t), allocatable :: df(:)

    ! Construct the ad_solver_t

    ! Build the grid

    call build_grid(ml, omega, gr_p, md_p, os_p, s, x)

    n_k = SIZE(s)

    ! Initialize the boundary conditions

    bd_i = ad_bound_t(ml, .TRUE., md_p, os_p)
    bd_i = ad_bound_t(ml, .FALSE., md_p, os_p)

    ! Initialize the difference equations

    allocate(df(n_k-1))

    do k = 1, n_k-1
       df(k) = ad_diff_t(ml, s(k), x(k), s(k+1), x(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bep_solver_t

    bs%r_bep_solver_t = r_bep_solver_t(bd_i, bd_o, df, MINVAL(omega), MAXVAL(omega), nm_p) 

    ! Other initializations

    bs%s = s
    bs%x = x

    bs%md_p = md_p
    bs%os_p = os_p

    ! Finish

    return

  end function ad_bep_solver_t_

  !****

  function sol_t_ (bs, omega) result (sl)

    class(ad_bep_solver_t), intent(inout) :: this
    real(WP), intent(in)                  :: omega
    class(sol_t)                          :: sl

    real(WP)      :: y(4,this%n_k)
    type(r_ext_t) :: discrim
    integer       :: k
    complex(WP)   :: y_c(6,this%n_k)

    ! Construct a sol_t from the ad_bvp_t

    ! Solve the BEP

    call this%solve(omega, y, discrim)

    ! Calculate the canonical 6-variable solution

    !$OMP PARALLEL DO 
    do k = 1, this%n_k
       y_c(1:4,k) = MATMUL(this%vr%B(this%s(k), this%x(k), omega), y(:,k))
       y_c(5:6,k) = 0._WP
    end do

    ! Construct the sol_t

    sl = sol_t(this%ml, this%s, this%x, y_c, this%md_p, this%os_p)

    ! Finish

    return

  end function sol_t_

end module gyre_ad_bvp

