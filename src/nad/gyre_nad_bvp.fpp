! Module   : gyre_nad_bvp
! Purpose  : nonadiabatic boundary value problem solver
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

module gyre_nad_bvp

  ! Uses

  use core_kinds

  use gyre_nad_bound
  use gyre_nad_diff
  use gyre_nad_share
  use gyre_nad_trans
  use gyre_bvp
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_model
  use gyre_mode
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_point

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_bvp_t) :: nad_bvp_t
     private
     type(nad_share_t), pointer :: sh => null()
     type(grid_t)               :: gr
     type(nad_trans_t)          :: tr
     type(mode_par_t)           :: md_p
     type(osc_par_t)            :: os_p
   contains
     private
     final             :: finalize_
     procedure, public :: set_omega_r
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

  function nad_bvp_t_ (ml, gr, md_p, nm_p, os_p) result (bp)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    type(mode_par_t), intent(in)        :: md_p
    type(num_par_t), intent(in)         :: nm_p
    type(osc_par_t), intent(in)         :: os_p
    type(nad_bvp_t)                     :: bp

    type(point_t)                 :: pt_i
    type(point_t)                 :: pt_o
    type(nad_share_t), pointer    :: sh
    type(nad_bound_t)             :: bd
    integer                       :: k
    type(nad_diff_t), allocatable :: df(:)

    ! Construct the nad_bvp_t

    pt_i = gr%pt(1)
    pt_o = gr%pt(gr%n_k)

    ! Initialize the shared data

    allocate(sh)

    sh = nad_share_t(ml, pt_i, pt_o, md_p, os_p)

    ! Initialize the boundary conditions

    bd = nad_bound_t(sh, pt_i, pt_o, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    !$OMP PARALLEL DO
    do k = 1, gr%n_k-1
       df(k) = nad_diff_t(sh, pt_i, gr%pt(k), gr%pt(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%c_bvp_t = c_bvp_t_(bd, df, nm_p) 

    ! Other initializations

    bp%sh => sh
    bp%gr = gr

    bp%tr = nad_trans_t(sh, pt_i, md_p, os_p)
    call bp%tr%stencil(gr%pt)

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function nad_bvp_t_

  !****

  subroutine finalize_ (this)

    type(nad_bvp_t), intent(inout) :: this

    ! Finalize the nad_bvp_t

    deallocate(this%sh)

    ! Finish

    return

  end subroutine finalize_

  !****

  subroutine set_omega_r (this, omega_r)

    class(nad_bvp_t), intent(inout) :: this
    real(WP), intent(in)            :: omega_r

    ! Set the real frequency to be used in rotation evaluations

    call this%sh%set_omega_r(omega_r)

    ! Finish

    return

  end subroutine set_omega_r

  !****

  function mode_t_ (bp, omega, j) result (md)

    class(nad_bvp_t), intent(inout) :: bp
    complex(WP), intent(in)         :: omega
    integer, intent(in)             :: j
    type(mode_t)                    :: md

    complex(WP)   :: y(6,bp%n_k)
    type(c_ext_t) :: discrim
    integer       :: k
    complex(WP)   :: y_c(6,bp%n_k)

    ! Calculate the solution vector

    call bp%build(omega)

    y = bp%soln_vec_hom()
    discrim = bp%det()

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do k = 1, bp%n_k

       call bp%tr%trans_vars(y(:,k), k, omega)

       y_c(:,k) = y(:,k)

    end do

    ! Construct the mode_t

    md = mode_t(omega, y_c, discrim, bp%sh%ml, bp%gr, bp%md_p, bp%os_p, j)

    ! Finish

    return

  end function mode_t_
  
end module gyre_nad_bvp
