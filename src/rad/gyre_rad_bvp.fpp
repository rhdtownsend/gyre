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
  use gyre_ext
  use gyre_grid
  use gyre_model
  use gyre_mode
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

  type, extends (r_bvp_t) :: rad_bvp_t
     class(model_t), pointer :: ml => null()
     type(grid_t)            :: gr
     type(rad_eqns_t)        :: eq
     type(rad_vars_t)        :: vr
     type(mode_par_t)        :: md_p
     type(osc_par_t)         :: os_p
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

  function rad_bvp_t_ (ml, gr, md_p, nm_p, os_p) result (bp)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    type(mode_par_t), intent(in)        :: md_p
    type(num_par_t), intent(in)         :: nm_p
    type(osc_par_t), intent(in)         :: os_p
    type(rad_bvp_t)                     :: bp

    type(rad_bound_t)             :: bd
    integer                       :: k
    type(rad_diff_t), allocatable :: df(:)

    ! Construct the rad_bvp_t

    if (os_p%cowling_approx) then
       $WARN(cowling_approx is ignored in 2nd-order radial equations)
    endif
    
    ! Initialize the boundary conditions

    bd = rad_bound_t(ml, gr, md_p, os_p)

    ! Initialize the difference equations

    allocate(df(gr%n_k-1))

    do k = 1, gr%n_k-1
       df(k) = rad_diff_t(ml, gr, k, md_p, nm_p, os_p)
    end do

    ! Initialize the bvp_t

    bp%r_bvp_t = r_bvp_t_(bd, df, nm_p) 

    ! Other initializations

    bp%ml => ml
    bp%gr = gr

    bp%eq = rad_eqns_t(ml, gr, md_p, os_p)
    bp%vr = rad_vars_t(ml, gr, md_p, os_p)

    bp%md_p = md_p
    bp%os_p = os_p

    ! Finish

    return

  end function rad_bvp_t_

  !****

  function mode_t_ (bp, j, omega) result (md)

    class(rad_bvp_t), intent(inout) :: bp
    integer, intent(in)             :: j
    real(WP), intent(in)            :: omega
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

    class(rad_bvp_t), intent(inout) :: bp
    real(WP), intent(in)            :: omega
    type(soln_t)                    :: sl

    real(WP)      :: y(2,bp%n_k)
    type(r_ext_t) :: discrim
    integer       :: k
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

    ! Calculate the solution vector and discriminant

    call bp%build(omega)

    y = bp%soln_vec()
    discrim = bp%det()

    ! Calculate its derivatives (nb: the vacuum check prevents errors, but
    ! leads to incorrect values for dy_dx)

    !$OMP PARALLEL DO PRIVATE (xA)
    do k = 1, bp%n_k

       associate (pt => bp%gr%pt(k))

         if (bp%ml%is_vacuum(pt)) then
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

         if (bp%ml%is_vacuum(pt)) then
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

    ! Calculate the gravity perturbation y_g and its derivative (nb:
    ! the vacuum check prevents errors, but leads to incorrect values
    ! for dy_g_dx)

    !$OMP PARALLEL DO PRIVATE (U, dU)
    do k = 1, bp%n_k

       associate (pt => bp%gr%pt(k))

         U = bp%ml%coeff(I_U, pt)
         
         y_g(k) = -U*y_c(1,k)

         if (bp%ml%is_vacuum(pt)) then
            dU = 0._WP
         else
            dU = bp%ml%dcoeff(I_U, pt)
         endif

         if (pt%x /= 0._WP) then
            dy_g_dx(k) = -U*dy_c_dx(1,k) - U*dU*y_c(1,k)/pt%x
         else
            dy_g_dx(k) = 0._WP
         endif

       end associate

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

end module gyre_rad_bvp
