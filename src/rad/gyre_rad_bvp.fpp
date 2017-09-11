! Module   : gyre_rad_bvp
! Purpose  : adiabatic radial bounary value problem solver
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

module gyre_rad_bvp

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_context
  use gyre_ext
  use gyre_grid
  use gyre_model
  use gyre_mode
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_qad_eval
  use gyre_rad_bound
  use gyre_rad_diff
  use gyre_rad_trans
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_bvp_t) :: rad_bvp_t
     private
     type(context_t), pointer :: cx => null()
     type(grid_t)             :: gr
     type(rad_trans_t)        :: tr
     type(qad_eval_t)         :: qe
     type(mode_par_t)         :: md_p
     type(osc_par_t)          :: os_p
  end type rad_bvp_t

  ! Interfaces

  interface rad_bvp_t
     module procedure rad_bvp_t_
  end interface rad_bvp_t

  interface mode_t
     module procedure mode_t_
  end interface mode_t

  ! Access specifiers

  private

  public :: rad_bvp_t
  public :: mode_t

  ! Procedures

contains

  function rad_bvp_t_ (cx, gr, md_p, nm_p, os_p) result (bp)

    type(context_t), pointer, intent(in) :: cx
    type(grid_t), intent(in)             :: gr
    type(mode_par_t), intent(in)         :: md_p
    type(num_par_t), intent(in)          :: nm_p
    type(osc_par_t), intent(in)          :: os_p
    type(rad_bvp_t)                      :: bp

    type(point_t)                 :: pt_i
    type(point_t)                 :: pt_o
    type(rad_bound_t)             :: bd
    integer                       :: k
    type(rad_diff_t), allocatable :: df(:)
    type(osc_par_t)               :: qad_os_p

    ! Construct the rad_bvp_t

    if (os_p%cowling_approx) then
       $WARN(cowling_approx is ignored in 2nd-order radial equations)
    endif
    
    pt_i = gr%pt(1)
    pt_o = gr%pt(gr%n_k)

    ! Initialize the boundary conditions

    bd = rad_bound_t(cx, pt_i, pt_o, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    !$OMP PARALLEL DO
    do k = 1, gr%n_k-1
       df(k) = rad_diff_t(cx, pt_i, gr%pt(k), gr%pt(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t_(bd, df, nm_p) 

    ! Other initializations

    bp%cx => cx
    bp%gr = gr

    bp%tr = rad_trans_t(cx, pt_i, md_p, os_p)
    call bp%tr%stencil(gr%pt)

    if (os_p%quasiad_eigfuncs) then
       qad_os_p = os_p
       qad_os_p%variables_set = 'GYRE'
       bp%qe = qad_eval_t(cx, gr, md_p, qad_os_p)
    endif

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function rad_bvp_t_

  !****

  function mode_t_ (bp, omega, j) result (md)

    class(rad_bvp_t), intent(inout) :: bp
    real(WP), intent(in)            :: omega
    integer, intent(in)             :: j
    type(mode_t)                    :: md

    real(WP)      :: y(2,bp%n_k)
    type(r_ext_t) :: discrim
    integer       :: k
    real(WP)      :: y_g(4,bp%n_k)
    real(WP)      :: U
    complex(WP)   :: y_c(6,bp%n_k)

    ! Calculate the solution vector

    call bp%build(omega)

    y = bp%soln_vec_hom()
    discrim = bp%det()

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do k = 1, bp%n_k
       call bp%tr%trans_vars(y(:,k), k, omega, from=.FALSE.)
    end do

    ! Set up gravitational eigenfunctions

    y_g(1:2,:) = y

    !$OMP PARALLEL DO PRIVATE (U)
    do k = 1, bp%n_k

       associate (pt => bp%gr%pt(k))
         U = bp%cx%ml%coeff(I_U, pt)
       end associate

       y_g(4,k) = -U*y(1,k)
       
    end do

    y_g(3,:) = 0._WP

    ! Set up complex eigenfunctions

    if (bp%os_p%quasiad_eigfuncs) then

       y_c = bp%qe%y_qad(omega, y_g)

    else

       y_c(1:4,:) = y_g
       y_c(5:6,:) = 0._WP

       call bp%cx%set_omega_ad(omega)

    endif

    ! Construct the mode_t

    md = mode_t(CMPLX(omega, KIND=WP), y_c, c_ext_t(discrim), bp%cx, bp%gr, bp%md_p, bp%os_p, j)

    ! Finish

    return

  end function mode_t_

end module gyre_rad_bvp
