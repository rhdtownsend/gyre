! Module   : gyre_nad_bep
! Purpose  : nonadiabatic boundary eigenvalue problem solver
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_nad_bep

  ! Uses

  use core_kinds

  use gyre_nad_bound
  use gyre_nad_diff
  use gyre_nad_vars
  use gyre_bep
  use gyre_ext
  use gyre_grid
  use gyre_grid_par
  use gyre_model
  use gyre_mode
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_bep_t) :: nad_bep_t
     class(model_t), pointer :: ml => null()
     type(mode_par_t)        :: md_p
     type(osc_par_t)         :: os_p
     type(nad_vars_t)        :: vr
     integer, allocatable    :: s(:)
     real(WP), allocatable   :: x(:)
  end type nad_bep_t

  ! Interfaces

  interface nad_bep_t
     module procedure nad_bep_t_
  end interface nad_bep_t

  interface mode_t
     module procedure mode_t_
  end interface mode_t

  ! Access specifiers

  private

  public :: nad_bep_t
  public :: mode_t

  ! Procedures

contains

  function nad_bep_t_ (ml, omega, gr_p, md_p, nm_p, os_p) result (bp)

    class(model_t), pointer, intent(in) :: ml
    real(WP), intent(in)                :: omega(:)
    type(grid_par_t), intent(in)        :: gr_p
    type(mode_par_t), intent(in)        :: md_p
    type(num_par_t), intent(in)         :: nm_p
    type(osc_par_t), intent(in)         :: os_p
    type(nad_bep_t)                     :: bp

    integer, allocatable          :: s(:)
    real(WP), allocatable         :: x(:)
    integer                       :: n_k
    type(nad_bound_t)             :: bd
    integer                       :: k
    type(nad_diff_t), allocatable :: df(:)
    real(WP)                      :: omega_min
    real(WP)                      :: omega_max

    ! Construct the nad_bep_t

    ! Build the grid

    call build_grid(ml, omega, gr_p, md_p, os_p, s, x)

    n_k = SIZE(s)

    ! Initialize the boundary conditions

    bd = nad_bound_t(ml, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(n_k-1))

    do k = 1, n_k-1
       df(k) = nad_diff_t(ml, s(k), x(k), s(k+1), x(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bep_t

    if (nm_p%restrict_roots) then
       omega_min = MINVAL(omega)
       omega_max = MAXVAL(omega)
    else
       omega_min = -HUGE(0._WP)
       omega_max = HUGE(0._WP)
    endif
    
    bp%c_bep_t = c_bep_t_(bd, df, omega_min, omega_max, nm_p) 

    ! Other initializations

    bp%ml => ml

    bp%vr = nad_vars_t(ml, md_p, os_p)

    bp%s = s
    bp%x = x

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function nad_bep_t_

  !****

  function mode_t_ (bp, omega) result (md)

    class(nad_bep_t), intent(inout) :: bp
    complex(WP), intent(in)         :: omega
    type(mode_t)                    :: md

    complex(WP)   :: y(6,bp%n_k)
    type(c_ext_t) :: discrim
    integer       :: k
    complex(WP)   :: y_c(6,bp%n_k)

    ! Construct a mode_t from the ad_bep_t

    ! Solve the BEP

    call bp%solve(omega, y, discrim)

    ! Calculate the canonical 6-variable solution

    !$OMP PARALLEL DO 
    do k = 1, bp%n_k
       y_c(:,k) = MATMUL(bp%vr%H(bp%s(k), bp%x(k), omega), y(:,k))
    end do

    ! Construct the mode_t

    md = mode_t(bp%ml, bp%s, bp%x, y_c, CMPLX(omega, KIND=WP), bp%n_k, bp%md_p, bp%os_p)

    md%discrim = discrim

    ! Finish

    return

  end function mode_t_

end module gyre_nad_bep
