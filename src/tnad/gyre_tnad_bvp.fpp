! Module   : gyre_tnad_bvp
! Purpose  : nonadiabatic (+turbulent convection) boundary value problem solver
!
! Copyright 2021-2022 Rich Townsend & The GYRE Team
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

module gyre_tnad_bvp

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_context
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_model
  use gyre_mode_par
  use gyre_nad_bound
  use gyre_nad_eqns
  use gyre_nad_trans
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_state
  use gyre_tnad_diff
  use gyre_tnad_eqns
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_bvp_t) :: tnad_bvp_t
     private
     type(context_t), pointer :: cx => null()
     type(grid_t)             :: gr
     type(tnad_eqns_t)        :: eq
     type(nad_trans_t)        :: tr
     type(mode_par_t)         :: md_p
     type(num_par_t)          :: nm_p
     type(osc_par_t)          :: os_p
  end type tnad_bvp_t

  ! Interfaces

  interface tnad_bvp_t
     module procedure tnad_bvp_t_
  end interface tnad_bvp_t

  interface wave_t
     module procedure wave_t_hom_
     module procedure wave_t_inhom_
  end interface wave_t

  ! Access specifiers

  private

  public :: tnad_bvp_t
  public :: wave_t

  ! Procedures

contains

  function tnad_bvp_t_ (cx, gr, md_p, nm_p, os_p) result (bp)

    class(context_t), pointer, intent(in) :: cx
    type(grid_t), intent(in)              :: gr
    type(mode_par_t), intent(in)          :: md_p
    type(num_par_t), intent(in)           :: nm_p
    type(osc_par_t), intent(in)           :: os_p
    type(tnad_bvp_t)                      :: bp
    
    type(nad_bound_t)              :: bd
    integer                        :: p
    type(tnad_diff_t), allocatable :: df(:)

    ! Construct the tnad_bvp_t

    ! Initialize the boundary conditions

    bd = nad_bound_t(cx, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_p-1))

    !$OMP PARALLEL DO
    do p = 1, gr%n_p-1
       df(p) = tnad_diff_t(cx, gr%pt(p), gr%pt(p+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%c_bvp_t = c_bvp_t(bd, df, nm_p) 

    ! Other initializations

    bp%cx => cx
    bp%gr = gr

    bp%eq = tnad_eqns_t(cx, md_p, os_p)
    call bp%eq%stencil(gr%pt)

    bp%tr = nad_trans_t(cx, md_p, os_p)
    call bp%tr%stencil(gr%pt)

    bp%md_p = md_p
    bp%nm_p = nm_p
    bp%os_p = os_p

    ! Finish

    return

  end function tnad_bvp_t_

  !****

  function wave_t_hom_ (bp, st, j) result (wv)

    class(tnad_bvp_t), intent(inout) :: bp
    type(c_state_t), intent(in)      :: st
    integer, intent(in)              :: j
    type(wave_t)                     :: wv

    complex(WP) :: y(6,bp%n_p)
    integer     :: p

    ! Calculate the homogeneous solution vector

    call bp%build(st)
    call bp%factor()

    y = bp%soln_vec_hom()

    !$OMP PARALLEL DO
    do p = 1, bp%n_p
       call bp%tr%trans_vars(y(:,p), p, st, from=.FALSE.)
    end do

    ! Construct the wave_t

    wv = wave_t_y_(bp, st, y, j)

    ! Finish

    return

  end function wave_t_hom_
  
  !****

  function wave_t_inhom_ (bp, st, z_i, z_o, j) result (wv)

    class(tnad_bvp_t), intent(inout) :: bp
    type(c_state_t), intent(in)      :: st
    complex(WP), intent(in)          :: z_i(:)
    complex(WP), intent(in)          :: z_o(:)
    integer, intent(in)              :: j
    type(wave_t)                     :: wv

    complex(WP) :: y(6,bp%n_p)
    integer     :: p

    $CHECK_BOUNDS(SIZE(z_i),bp%n_i)
    $CHECK_BOUNDS(SIZE(z_o),bp%n_o)

    ! Calculate the inhomogeneous solution vector

    call bp%build(st)
    call bp%factor()

    y = bp%soln_vec_inhom(z_i, z_o)

    !$OMP PARALLEL DO
    do p = 1, bp%n_p
       call bp%tr%trans_vars(y(:,p), p, st, from=.FALSE.)
    end do

    ! Construct the wave_t

    wv = wave_t_y_(bp, st, y, j)

    ! Finish

    return

  end function wave_t_inhom_

  !****

  function wave_t_y_ (bp, st, y, j) result (wv)

    class(tnad_bvp_t), intent(inout) :: bp
    type(c_state_t), intent(in)      :: st
    complex(WP), intent(in)          :: y(:,:)
    integer, intent(in)              :: j
    type(wave_t)                     :: wv

    integer       :: p
    complex(WP)   :: y_c(bp%n_e,bp%n_p)
    type(c_ext_t) :: discrim

    $CHECK_BOUNDS(SIZE(y, 1),bp%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),bp%n_p)

    ! Subtract off the turbulent force variable

    !$OMP PARALLEL DO
    do p = 1, bp%gr%n_p

       y_c(1,p) = y(1,p)
       y_c(2,p) = y(2,p) - bp%eq%y_trb(p, st, y(:,p))
       y_c(3:6,p) = y(3:6,p)

    end do

    ! Construct the wave_t

    discrim = bp%det()

    wv = wave_t(st, y_c, discrim, bp%cx, bp%gr, bp%md_p, bp%nm_p, bp%os_p, j)

    ! Finish

    return

  end function wave_t_y_

end module gyre_tnad_bvp
