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
     procedure :: eval_c_ => eval_gamma_func_
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

    type(gamma_func_t) :: gf
    real(WP)           :: gf_a
    real(WP)           :: gf_b
    integer            :: k

    ! Find the cell index and location of the inner turning point at
    ! frequency omega

    x_turn = HUGE(0._WP)

    gf%ml => ml
    allocate(gf%rt, SOURCE=r_rot_t(ml, md_p, os_p))
    gf%omega = omega

    gf%s = s(1)
    gf_b = gf%eval(x(1))

    turn_loop : do k = 1, SIZE(x)-2

       gf_a = gf_b

       gf%s = s(k+1)
       gf_b = gf%eval(x(k+1))
       
       if (gf_a > 0._WP .AND. gf_b <= 0._WP) then

          k_turn = k

          if (s(k) == s(k+1)) then

             if (ABS(gf_a) < EPSILON(0._WP)*ABS(gf_b)) then
                x_turn = x(k_turn)
             elseif (ABS(gf_b) < EPSILON(0._WP)*ABS(gf_a)) then
                x_turn = x(k_turn+1)
             else
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

  function eval_gamma_func_ (this, z) result (gamma)

    class(gamma_func_t), intent(inout) :: this
    complex(WP), intent(in)            :: z
    complex(WP)                        :: gamma

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: g_4
    real(WP) :: g_2
    real(WP) :: g_0

    ! Calculate the propagation discriminant

    associate (s => this%s, &
               x => REAL(z))

      V_g = this%ml%V_2(s, x)*x**2/this%ml%Gamma_1(s, x)
      As = this%ml%As(s, x)
      U = this%ml%U(s, x)
      c_1 = this%ml%c_1(s, x)

      lambda = this%rt%lambda(s, x, this%omega)

      g_4 = -4._WP*V_g*c_1
      g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
      g_0 = -4._WP*lambda*As/c_1

      gamma = (g_4*this%omega**4 + g_2*this%omega**2 + g_0)/this%omega**2

    end associate

    ! Finish

    return

  end function eval_gamma_func_

end module gyre_grid_util
