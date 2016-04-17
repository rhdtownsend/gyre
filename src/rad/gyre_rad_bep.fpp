! Module   : gyre_rad_bvp
! Purpose  : adiabatic radial bounary eigenvalue problem solver
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

module gyre_rad_bep

  ! Uses

  use core_kinds

  use gyre_bep
  use gyre_ext
  use gyre_grid
  use gyre_grid_build
  use gyre_grid_par
  use gyre_model
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_rad_bound
  use gyre_rad_diff
  use gyre_rad_eqns
  use gyre_rad_vars
  use gyre_soln
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_bep_t) :: rad_bep_t
     class(model_t), pointer :: ml => null()
     type(grid_t)            :: gr
     type(rad_eqns_t)        :: eq
     type(mode_par_t)        :: md_p
     type(osc_par_t)         :: os_p
     type(rad_vars_t)        :: vr
  end type rad_bep_t

  ! Interfaces

  interface rad_bep_t
     module procedure rad_bep_t_
  end interface rad_bep_t

  interface soln_t
     module procedure soln_t_
  end interface soln_t

  ! Access specifiers

  private

  public :: rad_bep_t
  public :: soln_t

  ! Procedures

contains

  function rad_bep_t_ (ml, omega, gr_p, md_p, nm_p, os_p) result (bp)

    class(model_t), pointer, intent(in) :: ml
    real(WP), intent(in)                :: omega(:)
    type(grid_par_t), intent(in)        :: gr_p
    type(mode_par_t), intent(in)        :: md_p
    type(num_par_t), intent(in)         :: nm_p
    type(osc_par_t), intent(in)         :: os_p
    type(rad_bep_t)                     :: bp

    type(grid_t)                  :: gr
    integer                       :: n_k
    type(rad_bound_t)             :: bd
    integer                       :: k
    type(rad_diff_t), allocatable :: df(:)
    real(WP)                      :: omega_min
    real(WP)                      :: omega_max

    ! Construct the rad_bep_t

    ! Build the grid

    call build_grid(ml, omega, gr_p, md_p, os_p, gr, verbose=.TRUE.)

    n_k = gr%n_k()

    ! Initialize the boundary conditions

    bd = rad_bound_t(ml, gr, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(n_k-1))

    do k = 1, n_k-1
       df(k) = rad_diff_t(ml, gr%pt(k), gr%pt(k+1), md_p, nm_p, os_p)
    end do

    ! Initialize the bep_t

    if (nm_p%restrict_roots) then
       omega_min = MINVAL(omega)
       omega_max = MAXVAL(omega)
    else
       omega_min = -HUGE(0._WP)
       omega_max = HUGE(0._WP)
    endif
    
    bp%r_bep_t = r_bep_t_(bd, df, omega_min, omega_max, nm_p) 

    ! Other initializations

    bp%ml => ml
    bp%gr = gr

    bp%eq = rad_eqns_t(ml, md_p, os_p)
    bp%vr = rad_vars_t(ml, md_p, os_p)

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function rad_bep_t_

  !****

  function soln_t_ (bp, omega) result (sl)

    class(rad_bep_t), intent(inout) :: bp
    real(WP), intent(in)            :: omega
    type(soln_t)                    :: sl

    real(WP)      :: y(2,bp%n_k)
    type(r_ext_t) :: discrim
    integer       :: k
    type(point_t) :: pt
    real(WP)      :: xA(2,2)
    real(WP)      :: dy_dx(2,bp%n_k)
    real(WP)      :: H(2,2)
    real(WP)      :: dH(2,2)
    complex(WP)   :: y_c(2,bp%n_k)
    complex(WP)   :: dy_c_dx(2,bp%n_k)
    complex(WP)   :: y_g(bp%n_k)
    complex(WP)   :: dy_g_dx(bp%n_k)
    real(WP)      :: U
    real(WP)      :: dU
    integer       :: i

    ! Calculate the solution vector y

    call bp%solve(omega, y, discrim)

    ! Calculate its derivatives (nb: the vacuum check prevents errors, but
    ! leads to incorrect values for dy_dx)

    !$OMP PARALLEL DO PRIVATE (pt, xA)
    do k = 1, bp%n_k

       pt = bp%gr%pt(k)

       if (bp%ml%vacuum(pt)) then
          xA = 0._WP
       else
          xA = bp%eq%xA(pt, omega)
       endif

       if (pt%x /= 0._WP) then
          dy_dx(:,k) = MATMUL(xA, y(:,k))/pt%x
       else
          dy_dx(:,k) = 0._WP
       endif

    end do

    ! Convert to canonical form (nb: the vacuum check prevents errors, but
    ! leads to incorrect values for dy_c_dx)

    !$OMP PARALLEL DO PRIVATE (pt, H, dH)
    do k = 1, bp%n_k

       pt = bp%gr%pt(k)

       H = bp%vr%H(pt, omega)

       y_c(:,k) = MATMUL(H, y(:,k))

       if (bp%ml%vacuum(pt)) then
          dH = 0._WP
       else
          dH = bp%vr%dH(pt, omega)
       endif

       if (pt%x /= 0._WP) then
          dy_c_dx(:,k) = MATMUL(dH/pt%x, y(:,k)) + MATMUL(H, dy_dx(:,k))
       else
          dy_c_dx(:,k) = 0._WP
       endif

    end do

    ! Calculate the gravity perturbation y_g and its derivative (nb:
    ! the vacuum check prevents errors, but leads to incorrect values
    ! for dy_g_dx)

    !$OMP PARALLEL DO PRIVATE (pt, U, dU)
    do k = 1, bp%n_k

       pt = bp%gr%pt(k)

       U = bp%ml%U(pt)

       y_g(k) = -U*y_c(1,k)

       if (bp%ml%vacuum(pt)) then
          dU = 0._WP
       else
          dU = bp%ml%dU(pt)
       endif

       if (pt%x /= 0._WP) then
          dy_g_dx(k) = -U*dy_c_dx(1,k) - U*dU*y_c(1,k)/pt%x
       else
          dy_g_dx(k) = 0._WP
       endif

    end do

    ! Construct the soln_t

    sl = soln_t(bp%gr, CMPLX(omega, KIND=WP), c_ext_t(discrim))

    do i = 1, 2
       call sl%set_y(i, y_c(i,:), dy_c_dx(i,:))
    end do

    call sl%set_y(4, y_g, dy_g_dx)

    call sl%set_y(5, SPREAD(CMPLX(0._WP, KIND=WP), 1, bp%n_k), &
                     SPREAD(CMPLX(0._WP, KIND=WP), 1, bp%n_k))

    ! Finish

    return

  end function soln_t_

end module gyre_rad_bep
