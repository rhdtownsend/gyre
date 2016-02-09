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

  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
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

  subroutine find_turn (ml, s, x, omega, md_p, os_p, k_turn, x_turn)

    class(model_t), pointer, intent(in) :: ml
    integer, intent(in)                 :: s(:)
    real(WP), intent(in)                :: x(:)
    real(WP), intent(in)                :: omega
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    integer, intent(out)                :: k_turn
    real(WP), intent(out)               :: x_turn

    integer                     :: n_k
    class(r_rot_t), allocatable :: rt
    real(WP)                    :: gamma_a
    real(WP)                    :: gamma_b
    integer                     :: k
    type(gamma_func_t)          :: gf

    ! Find the cell index and location of the inner turning point at
    ! frequency omega

    n_k = SIZE(s)

    k_turn = n_k
    x_turn = HUGE(0._WP)

    allocate(rt, SOURCE=r_rot_t(ml, md_p, os_p))

    gamma_b = gamma_(ml, rt, s(1), x(1), omega)

    turn_loop : do k = 1, n_k-1

       gamma_a = gamma_b
       gamma_b = gamma_(ml, rt, s(k+1), x(k+1), omega)

       if (gamma_a > 0._WP .AND. gamma_b <= 0._WP) then

          k_turn = k

          if (s(k) == s(k+1)) then

             if (ABS(gamma_a) < EPSILON(0._WP)*ABS(gamma_b)) then

                x_turn = x(k_turn)

             elseif (ABS(gamma_b) < EPSILON(0._WP)*ABS(gamma_a)) then

                x_turn = x(k_turn+1)

             else

                gf%ml => ml
                allocate(gf%rt, SOURCE=rt)
                gf%s = s(k)
                gf%omega = omega
                
                x_turn = gf%root(x(k), x(k+1), 0._WP)

             endif

          else

             x_turn = x(k_turn)

          end if

          exit turn_loop

       endif

    end do turn_loop

    ! Finish

    return

  end subroutine find_turn

  !****

  function gamma_ (ml, rt, s, x, omega) result (gamma)

    class(model_t), pointer, intent(in) :: ml
    class(r_rot_t), intent(in)          :: rt
    integer, intent(in)                 :: s
    real(WP), intent(in)                :: x
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

    U = ml%U(s, x)

    if (U > 0._WP) then

       V_g = ml%V_2(s, x)*x**2/ml%Gamma_1(s, x)
       As = ml%As(s, x)
       c_1 = ml%c_1(s, x)

       lambda = rt%lambda(s, x, omega)

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

    ! Evaluate the gamma_func_t

    gamma = gamma_(this%ml, this%rt, this%s, REAL(z), this%omega)

    ! Finish

    return

  end function eval_c_

end module gyre_grid_util
