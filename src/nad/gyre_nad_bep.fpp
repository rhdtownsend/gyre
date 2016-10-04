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
  use gyre_nad_eqns
  use gyre_nad_vars
  use gyre_bep
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_model
  use gyre_mode
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_soln

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_bep_t) :: nad_bep_t
     class(model_t), pointer :: ml => null()
     type(grid_t)            :: gr
     type(nad_eqns_t)        :: eq
     type(mode_par_t)        :: md_p
     type(osc_par_t)         :: os_p
     type(nad_vars_t)        :: vr
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

  function nad_bep_t_ (ml, gr, omega, md_p, nm_p, os_p) result (bp)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    real(WP), intent(in)                :: omega(:)
    type(mode_par_t), intent(in)        :: md_p
    type(num_par_t), intent(in)         :: nm_p
    type(osc_par_t), intent(in)         :: os_p
    type(nad_bep_t)                     :: bp

    type(nad_bound_t)             :: bd
    integer                       :: k
    type(nad_diff_t), allocatable :: df(:)
    real(WP)                      :: omega_min
    real(WP)                      :: omega_max

    ! Construct the nad_bep_t

    ! Initialize the boundary conditions

    bd = nad_bound_t(ml, gr, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    do k = 1, gr%n_k-1
       df(k) = nad_diff_t(ml, gr, k, md_p, nm_p, os_p)
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
    bp%gr = gr

    bp%eq = nad_eqns_t(ml, gr, md_p, os_p)
    bp%vr = nad_vars_t(ml, gr, md_p, os_p)

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function nad_bep_t_

  !****

  function mode_t_ (bp, j, omega) result (md)

    class(nad_bep_t), intent(inout) :: bp
    integer, intent(in)             :: j
    complex(WP), intent(in)         :: omega
    type(mode_t)                    :: md

    type(soln_t) :: sl

    ! Construct the mode_t

    sl = soln_t_(bp, omega)

    md = mode_t(bp%ml, sl, j, bp%md_p, bp%os_p)

    ! Finish

    return

  end function mode_t_

  !****

  function soln_t_ (bp, omega) result (sl)

    class(nad_bep_t), intent(inout) :: bp
    complex(WP), intent(in)         :: omega
    type(soln_t)                    :: sl

    complex(WP)   :: y(6,bp%n_k)
    type(c_ext_t) :: discrim
    integer       :: k
    complex(WP)   :: xA(6,6)
    complex(WP)   :: dy_dx(6,bp%n_k)
    complex(WP)   :: H(6,6)
    complex(WP)   :: dH(6,6)
    complex(WP)   :: y_c(6,bp%n_k)
    complex(WP)   :: dy_c_dx(6,bp%n_k)
    integer       :: i

    ! Calculate the solution vector y

    call bp%solve(omega, y, discrim)

    ! Calculate its derivatives (nb: the vacuum check prevents errors, but
    ! leads to incorrect values for dy_dx)

    !$OMP PARALLEL DO PRIVATE (xA)
    do k = 1, bp%n_k

       associate (pt => bp%gr%pt(k))

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

       end associate
         
    end do

    ! Convert to canonical form (nb: the vacuum check prevents errors, but
    ! leads to incorrect values for dy_c_dx)

    !$OMP PARALLEL DO PRIVATE (H, dH)
    do k = 1, bp%n_k

       associate (pt => bp%gr%pt(k))

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

       end associate

    end do

    ! Construct the soln_t

    sl = soln_t(bp%gr, omega, discrim)

    do i = 1, 6
       call sl%set_y(i, y_c(i,:), dy_c_dx(i,:))
    end do

    ! Finish

    return

  end function soln_t_

end module gyre_nad_bep
