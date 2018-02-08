! Module   : gyre_nad_bvp
! Purpose  : nonadiabatic boundary value problem solver
!
! Copyright 2013-2018 Rich Townsend
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
  use gyre_context
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_model
  use gyre_mode
  use gyre_mode_par
  use gyre_nad_bound
  use gyre_nad_diff
  use gyre_nad_trans
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_bvp_t) :: nad_bvp_t
     private
     type(context_t), pointer :: cx => null()
     type(grid_t)             :: gr
     type(nad_trans_t)        :: tr
     type(mode_par_t)         :: md_p
     type(osc_par_t)          :: os_p
  end type nad_bvp_t

  ! Interfaces

  interface nad_bvp_t
     module procedure nad_bvp_t_
  end interface nad_bvp_t

  interface mode_t
     module procedure mode_t_
  end interface mode_t

  ! Access specifiers

  private

  public :: nad_bvp_t
  public :: mode_t

  ! Procedures

contains

  function nad_bvp_t_ (cx, gr, md_p, nm_p, os_p) result (bp)

    class(context_t), pointer, intent(in) :: cx
    type(grid_t), intent(in)              :: gr
    type(mode_par_t), intent(in)          :: md_p
    type(num_par_t), intent(in)           :: nm_p
    type(osc_par_t), intent(in)           :: os_p
    type(nad_bvp_t)                       :: bp

    type(nad_bound_t)             :: bd
    integer                       :: k
    type(nad_diff_t), allocatable :: df(:)

    ! Construct the nad_bvp_t

    ! Initialize the boundary conditions

    bd = nad_bound_t(cx, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    !$OMP PARALLEL DO
    do k = 1, gr%n_k-1
       df(k) = nad_diff_t(cx, gr%pt(k), gr%pt(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%c_bvp_t = c_bvp_t_(bd, df, nm_p) 

    ! Other initializations

    bp%cx => cx
    bp%gr = gr

    bp%tr = nad_trans_t(cx, md_p, os_p)
    call bp%tr%stencil(gr%pt)

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function nad_bvp_t_

  !****

  function mode_t_ (bp, st, j) result (md)

    class(nad_bvp_t), intent(inout) :: bp
    type(c_state_t), intent(in)     :: st
    integer, intent(in)             :: j
    type(mode_t)                    :: md

    complex(WP)   :: y(6,bp%n_k)
    type(c_ext_t) :: discrim
    integer       :: k

    ! Calculate the solution vector

    call bp%build(st)

    y = bp%soln_vec_hom()
    discrim = bp%det()

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do k = 1, bp%n_k
       call bp%tr%trans_vars(y(:,k), k, st)
    end do

    ! Construct the mode_t

    md = mode_t(st, y, discrim, bp%cx, bp%gr, bp%md_p, bp%os_p, j)

    ! Finish

    return

  end function mode_t_

end module gyre_nad_bvp
