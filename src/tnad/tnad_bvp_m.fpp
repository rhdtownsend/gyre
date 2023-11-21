! Module  : tnad_bvp_m
! Purpose : nonadiabatic (+turbulent convection) boundary value problem solver
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

module tnad_bvp_m

  ! Uses

  use kinds_m

  use bvp_m
  use context_m
  use ext_m
  use grid_m
  use grid_factory_m
  use model_m
  use mode_par_m
  use nad_bound_m
  use nad_eqns_m
  use nad_trans_m
  use num_par_m
  use osc_par_m
  use point_m
  use state_m
  use tnad_diff_m
  use tnad_eqns_m
  use wave_m

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
    integer                        :: j
    type(tnad_diff_t), allocatable :: df(:)

    ! Construct the tnad_bvp_t

    ! Initialize the boundary conditions

    bd = nad_bound_t(cx, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n-1))

    !$OMP PARALLEL DO
    do j = 1, gr%n-1
       df(j) = tnad_diff_t(cx, gr%pt(j), gr%pt(j+1), md_p, nm_p, os_p)
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

  function wave_t_hom_ (bp, st, id) result (wv)

    class(tnad_bvp_t), intent(inout) :: bp
    type(c_state_t), intent(in)      :: st
    integer, intent(in)              :: id
    type(wave_t)                     :: wv

    complex(WP) :: y(6,bp%n)
    integer     :: j

    ! Calculate the homogeneous solution vector

    call bp%build(st)
    call bp%factor()

    y = bp%soln_vec_hom()

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

    class(tnad_bvp_t), intent(inout) :: bp
    type(c_state_t), intent(in)      :: st
    complex(WP), intent(in)          :: z_i(:)
    complex(WP), intent(in)          :: z_o(:)
    integer, intent(in)              :: id
    type(wave_t)                     :: wv

    complex(WP) :: y(6,bp%n)
    integer     :: j

    $CHECK_BOUNDS(SIZE(z_i),bp%n_i)
    $CHECK_BOUNDS(SIZE(z_o),bp%n_o)

    ! Calculate the inhomogeneous solution vector

    call bp%build(st)
    call bp%factor()

    y = bp%soln_vec_inhom(z_i, z_o)

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

    class(tnad_bvp_t), intent(inout) :: bp
    type(c_state_t), intent(in)      :: st
    complex(WP), intent(in)          :: y(:,:)
    integer, intent(in)              :: id
    type(wave_t)                     :: wv

    integer       :: j
    complex(WP)   :: y_c(bp%n_e,bp%n)
    type(c_ext_t) :: discrim

    $CHECK_BOUNDS(SIZE(y, 1),bp%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),bp%n)

    ! Subtract off the turbulent force variable

    !$OMP PARALLEL DO
    do j = 1, bp%gr%n

       y_c(1,j) = y(1,j)
       y_c(2,j) = y(2,j) - bp%eq%y_trb(j, st, y(:,j))
       y_c(3:6,j) = y(3:6,j)

    end do

    ! Construct the wave_t

    discrim = bp%det()

    wv = wave_t(st, y_c, discrim, bp%cx, bp%gr, bp%md_p, bp%nm_p, bp%os_p, id)

    ! Finish

    return

  end function wave_t_y_

end module tnad_bvp_m
