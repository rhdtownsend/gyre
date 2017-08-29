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
  use gyre_ad_eqns
  use gyre_ad_share
  use gyre_ad_trans
  use gyre_bvp
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_interp
  use gyre_model
  use gyre_mode
  use gyre_mode_par
  use gyre_nad_eqns
  use gyre_nad_share
  use gyre_nad_trans
  use gyre_num_par
  use gyre_osc_par
  use gyre_point

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_bvp_t) :: ad_bvp_t
     private
     type(ad_share_t), pointer  :: ad_sh => null()
     type(nad_share_t), pointer :: nad_sh => null()
     type(nad_eqns_t)           :: nad_eq
     type(nad_trans_t)          :: nad_tr
     type(grid_t)               :: gr
     type(mode_par_t)           :: md_p
     type(osc_par_t)            :: os_p
   contains
     private
     final :: finalize_
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
    type(ad_share_t), pointer    :: ad_sh
    type(nad_share_t), pointer   :: nad_sh
    type(ad_bound_t)             :: bd
    integer                      :: k
    type(ad_diff_t), allocatable :: df(:)

    ! Construct the ad_bvp_t

    pt_i = gr%pt(1)
    pt_o = gr%pt(gr%n_k)

    ! Initialize the shared data

    allocate(ad_sh)
    allocate(nad_sh)

    ad_sh = ad_share_t(ml, pt_i, md_p, os_p)
    nad_sh = nad_share_t(ml, pt_i, pt_o, md_p, os_p)

    ! Initialize the boundary conditions

    bd = ad_bound_t(ad_sh, pt_i, pt_o, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    !$OMP PARALLEL DO
    do k = 1, gr%n_k-1
       df(k) = ad_diff_t(ad_sh, pt_i, gr%pt(k), gr%pt(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t_(bd, df, nm_p) 

    ! Other initializations

    bp%ad_sh => ad_sh
    bp%nad_sh => nad_sh

    bp%nad_eq = nad_eqns_t(nad_sh, pt_i, md_p, os_p)
    call bp%nad_eq%stencil(gr%pt)

    bp%nad_tr = nad_trans_t(nad_sh, pt_i, md_p, os_p)
    call bp%nad_tr%stencil(gr%pt)

    bp%gr = gr

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function ad_bvp_t_

  !****

  subroutine finalize_ (this)

    type(ad_bvp_t), intent(inout) :: this

    ! Finalize the ad_bvp_t

    deallocate(this%ad_sh)
    deallocate(this%nad_sh)

    ! Finish

    return

  end subroutine finalize_

  !****

  function mode_t_ (bp, omega, j) result (md)

    class(ad_bvp_t), intent(inout) :: bp
    real(WP), intent(in)           :: omega
    integer, intent(in)            :: j
    type(mode_t)                   :: md

    real(WP)         :: y(4,bp%n_k)
    type(r_ext_t)    :: discrim
    integer          :: k
    complex(WP)      :: y_c(6,bp%n_k)
    integer          :: s
    complex(WP)      :: xA(6,6)
    complex(WP)      :: xA_5(6,bp%n_k)
    complex(WP)      :: xA_6(6,bp%n_k)
    complex(WP)      :: dy_6(bp%n_k)
    type(c_interp_t) :: in

    ! Calculate the solution vector

    call bp%build(omega)

    y = bp%soln_vec_hom()
    discrim = bp%det()

    ! Set up complex eigenfunctions

    y_c(1:4,:) = y

    if (bp%os_p%quasiad_eigfuncs) then

       ! First, evaluate components of the non-adiabatic RHS matrix
       ! corresponding to the energy conservation and transport
       ! equations

       call bp%nad_sh%set_omega_r(omega)

       !$OMP PARALLEL DO PRIVATE (xA)
       do k = 1, bp%n_k

          xA = bp%nad_eq%xA(k, CMPLX(omega, KIND=WP))

          xA_5(:,k) = xA(5,:)
          xA_6(:,k) = xA(6,:)

       end do

       ! Evaluate the luminosity perturbation eigenfunction

       where (bp%gr%pt%x /= 0._WP)
          y_c(6,:) = -(xA_5(1,:)*y_c(1,:) + xA_5(2,:)*y_c(2,:) + xA_5(3,:)*y_c(3,:) + xA_5(4,:)*y_c(4,:))/xA_5(6,:)
       elsewhere
          y_c(6,:) = 0._WP
       end where

       ! Evaluate the gradient of the luminosity perturbation
       ! eigenfunction, segment-by-segment (this is an overkill
       ! approach, shouldn't need to use c_interp_t)

       seg_loop : do s = bp%gr%s_i(), bp%gr%s_o()
          associate (k_i => bp%gr%k_i(s), k_o => bp%gr%k_o(s))
            in = c_interp_t(bp%gr%pt(k_i:k_o)%x, y_c(6,k_i:k_o), 'SPLINE')
            dy_6(k_i:k_o) = bp%gr%pt(k_i:k_o)%x*in%df_dx(bp%gr%pt(k_i:k_o)%x)
          end associate
       end do seg_loop

       ! Evaluate the entropy perturbation eigenfunction

       where (bp%gr%pt%x /= 0._WP)
          y_c(5,:) = (dy_6 - (xA_6(1,:)*y_c(1,:) + xA_6(2,:)*y_c(2,:) + xA_6(3,:)*y_c(3,:) + &
                              xA_6(4,:)*y_c(4,:) + xA_6(6,:)*y_c(6,:)))/xA_6(5,:)
       elsewhere
          y_c(5,:) = 0._WP
       end where

    else

       y_c(5:6,:) = 0._WP

    endif

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do k = 1, bp%n_k
       call bp%nad_tr%trans_vars(y_c(:,k), k, CMPLX(omega, KIND=WP), from=.FALSE.)
    end do

    ! Construct the mode_t

    md = mode_t(CMPLX(omega, KIND=WP), y_c, c_ext_t(discrim), bp%ad_sh%ml, bp%gr, bp%md_p, bp%os_p, j)

    ! Finish

    return

  end function mode_t_

end module gyre_ad_bvp
