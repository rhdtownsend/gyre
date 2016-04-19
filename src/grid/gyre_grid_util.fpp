! Module   : gyre_grid_util
! Purpose  : grid utilities
!
! Copyright 2016 Rich Townsend
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

module gyre_grid_util

  ! Uses

  use core_kinds
  use core_func

  use gyre_grid
  use gyre_grid_par
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions (used internally for root-finding)

  type, extends (func_t) :: gamma_func_t
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     real(WP)                    :: omega
     integer                     :: s
   contains
     procedure :: eval_c_
  end type gamma_func_t

  ! Access specifiers

  private

  public :: find_turn

  ! Procedures

contains

  subroutine find_turn (ml, gr, omega, md_p, os_p, k_turn, x_turn)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    real(WP), intent(in)                :: omega
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    integer, intent(out)                :: k_turn
    real(WP), intent(out)               :: x_turn

    class(r_rot_t), allocatable :: rt
    real(WP)                    :: gamma_a
    real(WP)                    :: gamma_b
    integer                     :: k
    type(gamma_func_t)          :: gf
    type(point_t)               :: pt_a
    type(point_t)               :: pt_b

    ! Find the cell index and abscissa of the inner turning point at
    ! frequency omega

    k_turn = gr%n_k
    x_turn = HUGE(0._WP)

    allocate(rt, SOURCE=r_rot_t(ml, gr, md_p, os_p))

    gamma_b = gamma_(ml, rt, gr%pt(1), omega)

    turn_loop : do k = 1, gr%n_k-1

       ! Check for a sign change in gamma

       gamma_a = gamma_b
       gamma_b = gamma_(ml, rt, gr%pt(k+1), omega)

       if (gamma_a > 0._WP .AND. gamma_b <= 0._WP) then

          k_turn = k

          pt_a = gr%pt(k)
          pt_b = gr%pt(k+1)

          if (pt_a%s == pt_b%s) then

             if (ABS(gamma_a) < EPSILON(0._WP)*ABS(gamma_b)) then

                x_turn = pt_a%x

             elseif (ABS(gamma_b) < EPSILON(0._WP)*ABS(gamma_a)) then

                x_turn = pt_b%x

             else

                gf%ml => ml
                allocate(gf%rt, SOURCE=rt)
                gf%s = pt_a%s
                gf%omega = omega
                
                x_turn = gf%root(pt_a%x, pt_b%x, 0._WP)

             endif

          else

             x_turn = pt_a%x

          end if

          exit turn_loop

       endif

    end do turn_loop

    ! Finish

    return

  end subroutine find_turn

  !****

  function gamma_ (ml, rt, pt, omega) result (gamma)

    class(model_t), pointer, intent(in) :: ml
    class(r_rot_t), intent(in)          :: rt
    type(point_t), intent(in)           :: pt
    real(WP), intent(in)                :: omega
    real(WP)                            :: gamma

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: g_4
    real(WP) :: g_2
    real(WP) :: g_0

    ! Calculate the propagation discriminant gamma (> 0 : propagation,
    ! < 0 : evanescence)

    U = ml%U(pt)

    if (U > 0._WP) then

       V_g = ml%V_2(pt)*pt%x**2/ml%Gamma_1(pt)
       As = ml%As(pt)
       c_1 = ml%c_1(pt)

       lambda = rt%lambda(pt, omega)

       g_4 = -4._WP*V_g*c_1
       g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
       g_0 = -4._WP*lambda*As/c_1

       gamma = (g_4*omega**4 + g_2*omega**2 + g_0)/omega**2

    else

       gamma = HUGE(0._WP)

    endif

    ! Finish

    return

  end function gamma_

  !****

  function eval_c_ (this, z) result (gamma)

    class(gamma_func_t), intent(inout) :: this
    complex(WP), intent(in)            :: z
    complex(WP)                        :: gamma

    type(point_t) :: pt

    ! Evaluate the gamma_func_t

    pt = point_t(this%s, REAL(z))

    gamma = gamma_(this%ml, this%rt, pt, this%omega)

    ! Finish

    return

  end function eval_c_

end module gyre_grid_util
