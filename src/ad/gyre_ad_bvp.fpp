! Module   : gyre_ad_bvp
! Purpose  : adiabatic bounary value problem solver
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

module gyre_ad_bvp

  ! Uses

  use core_kinds

  use gyre_ad_bound
  use gyre_ad_diff
  use gyre_ad_trans
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

  type, extends (r_bvp_t) :: ad_bvp_t
     class(model_t), pointer :: ml => null()
     type(grid_t)            :: gr
     type(ad_trans_t)        :: tr
     type(mode_par_t)        :: md_p
     type(osc_par_t)         :: os_p
  end type ad_bvp_t

  ! Interfaces

  interface ad_bvp_t
     module procedure ad_bvp_t_
  end interface ad_bvp_t

  interface mode_t
     module procedure mode_t_
  end interface mode_t

  ! Access specifiers

  private

  public :: ad_bvp_t
  public :: mode_t

  ! Procedures

contains

  function ad_bvp_t_ (ml, gr, md_p, nm_p, os_p) result (bp)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    type(mode_par_t), intent(in)        :: md_p
    type(num_par_t), intent(in)         :: nm_p
    type(osc_par_t), intent(in)         :: os_p
    type(ad_bvp_t)                      :: bp

    type(point_t)                :: pt_i
    type(point_t)                :: pt_o
    type(ad_bound_t)             :: bd
    integer                      :: k
    type(ad_diff_t), allocatable :: df(:)

    ! Construct the ad_bvp_t

    pt_i = gr%pt(1)
    pt_o = gr%pt(gr%n_k)

    ! Initialize the boundary conditions

    bd = ad_bound_t(ml, pt_i, pt_o, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    !$OMP PARALLEL DO
    do k = 1, gr%n_k-1
       df(k) = ad_diff_t(ml, pt_i, gr%pt(k), gr%pt(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t_(bd, df, nm_p) 

    ! Other initializations
 
    bp%ml => ml
    bp%gr = gr

    bp%tr = ad_trans_t(ml, pt_i, md_p, os_p)
    call bp%tr%stencil(gr%pt)

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function ad_bvp_t_

  !****

  function mode_t_ (bp, omega, j) result (md)

    class(ad_bvp_t), intent(inout) :: bp
    real(WP), intent(in)           :: omega
    integer, intent(in)            :: j
    type(mode_t)                   :: md

    real(WP)      :: y(4,bp%n_k)
    type(r_ext_t) :: discrim
    integer       :: k
    complex(WP)   :: y_c(6,bp%n_k)

    ! Calculate the solution vector

    call bp%build(omega)

    y = bp%soln_vec()
    discrim = bp%det()

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do k = 1, bp%n_k

       call bp%tr%trans_vars(y(:,k), k, omega, from=.FALSE.)

       y_c(1:4,k) = y(:,k)
       y_c(5:6,k) = 0._WP

    end do

    ! Construct the mode_t

    md = mode_t(CMPLX(omega, KIND=WP), y_c, c_ext_t(discrim), bp%ml, bp%gr, bp%md_p, bp%os_p, j)

    ! Finish

    return

  end function mode_t_

end module gyre_ad_bvp
