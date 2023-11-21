! Module  : cad_bvp_m
! Purpose : adiabatic bounary value problem solver (complex variables)
!
! Copyright 2013-2022 Rich Townsend
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

module cad_bvp_m

  ! Uses

  use kinds_m

  use cad_bound_m
  use cad_diff_m
  use cad_eqns_m
  use cad_trans_m
  use bvp_m
  use context_m
  use ext_m
  use grid_m
  use grid_factory_m
  use interp_m
  use model_m
  use mode_par_m
  use num_par_m
  use osc_par_m
  use point_m
  use state_m
  use wave_m

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_bvp_t) :: cad_bvp_t
     private
     type(context_t), pointer :: cx => null()
     type(cad_trans_t)        :: tr
     type(grid_t)             :: gr
     type(mode_par_t)         :: md_p
     type(osc_par_t)          :: os_p
  end type ad_bvp_t

  ! Interfaces

  interface cad_bvp_t
     module procedure cad_bvp_t_
  end interface cad_bvp_t

  interface wave_t
     module procedure wave_t_hom_
     module procedure wave_t_inhom_
  end interface wave_t

  ! Access specifiers

  private

  public :: cad_bvp_t
  public :: wave_t

  ! Procedures

contains

  function cad_bvp_t_ (cx, gr, md_p, nm_p, os_p) result (bp)

    type(context_t), pointer, intent(in) :: cx
    type(grid_t), intent(in)             :: gr
    type(mode_par_t), intent(in)         :: md_p
    type(num_par_t), intent(in)          :: nm_p
    type(osc_par_t), intent(in)          :: os_p
    type(cad_bvp_t)                      :: bp

    type(cad_bound_t)             :: bd
    integer                       :: j
    type(cad_diff_t), allocatable :: df(:)
    type(osc_par_t)               :: qad_os_p

    ! Construct the ad_bvp_t

    ! Initialize the boundary conditions

    bd = cad_bound_t(cx, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n-1))

    !$OMP PARALLEL DO
    do j = 1, gr%n-1
       df(j) = cad_diff_t(cx, gr%pt(j), gr%pt(j+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%r_bvp_t = c_bvp_t_(bd, df, nm_p)

    ! Other initializations

    bp%cx => cx
    bp%gr = gr

    bp%tr = cad_trans_t(cx, md_p, os_p)
    call bp%tr%stencil(gr%pt)

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function ad_bvp_t_

  !****

  function wave_t_hom_ (bp, st, id) result (wv)

    class(ad_bvp_t), intent(inout) :: bp
    type(r_state_t), intent(in)    :: st
    integer, intent(in)            :: id
    type(wave_t)                   :: wv

    real(WP) :: y(4,bp%n_p)
    integer  :: j

    ! Calculate the homogeneous solution vector

    call bp%build(st)
    call bp%factor()

    y = bp%soln_vec_hom()

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do j = 1, bp%n
       call bp%tr%trans_vars(y(:,j), j, st, from=.FALSE.)
    end do

    ! Construct the wave_t

    wv = wave_t_y_(bp, st, y, id)

    ! Finish

    return

  end function wave_t_hom_

  !****

  function wave_t_inhom_ (bp, st, z_i, z_o, id) result (wv)

    class(cad_bvp_t), intent(inout) :: bp
    type(c_state_t), intent(in)     :: st
    complex(WP), intent(in)         :: z_i(:)
    complex(WP), intent(in)         :: z_o(:)
    integer, intent(in)             :: id
    type(wave_t)                    :: wv

    complex(WP) :: y(4,bp%n_p)
    integer     :: j

    $CHECK_BOUNDS(SIZE(z_i),bp%n_i)
    $CHECK_BOUNDS(SIZE(z_o),bp%n_o)

    ! Calculate the inhomogeneous solution vector

    call bp%build(st)
    call bp%factor()

    y = bp%soln_vec_inhom(z_i, z_o)

    ! Convert to canonical form

    !$OMP PARALLEL DO
    do j = 1, bp%n
       call bp%tr%trans_vars(y(:,j), j, st, from=.FALSE.)
    end do

    ! Construct the wave_t

    wv = wave_t_y_(bp, st, y, id)

    ! Finish

    return

  end function wave_t_inhom_

  !****

  function wave_t_y_ (bp, st, y, id) result (wv)

    class(cad_bvp_t), intent(inout) :: bp
    type(c_state_t), intent(in)     :: st
    complex(WP), intent(in)         :: y(:,:)
    integer, intent(in)             :: id
    type(wave_t)                    :: wv

    type(c_state_t) :: st_c
    complex(WP)     :: y_c(6,bp%n_p)
    type(c_ext_t)   :: discrim

    $CHECK_BOUNDS(SIZE(y, 1),bp%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),bp%n_p)

    ! Set up complex eigenfunctions

    st_c = c_state_t(CMPLX(st%omega, KIND=WP), st%omega)

    y_c(1:4,:) = y
    y_c(5:6,:) = 0._WP

    ! Construct the wave_t

    discrim = c_ext_t(bp%det())

    wv = wave_t(st_c, y_c, discrim, bp%cx, bp%gr, bp%md_p, bp%os_p, id)

    ! Finish

    return

  end function wave_t_y_

end module cad_bvp_m
