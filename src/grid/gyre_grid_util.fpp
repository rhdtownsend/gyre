! Module   : gyre_grid_util
! Purpose  : grid utilities
!
! Copyright 2016-2020 Rich Townsend & The GYRE Team
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

  use gyre_context
  use gyre_grid
  use gyre_model
  use gyre_point
  use gyre_rot
  use gyre_rot_factory
  use gyre_state
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions (used internally for root-finding)

  type, extends (func_t) :: gamma_func_t
     type(context_t), pointer  :: cx
     class(r_state_t), pointer :: st
     integer                   :: s
     real (WP)                 :: a_g
     real (WP)                 :: a_p
   contains
     procedure :: eval_c_
  end type gamma_func_t

  ! Access specifiers

  private

  public :: find_turn

  ! Procedures

contains

  subroutine find_turn (cx, gr, st, k_turn, x_turn, a_g, a_p)

    type(context_t), target, intent(in)  :: cx
    type(grid_t), intent(in)             :: gr
    class(r_state_t), target, intent(in) :: st
    real(WP), intent(in), optional       :: a_g, a_p

    integer, intent(out)                 :: k_turn
    real(WP), intent(out)                :: x_turn

    real(WP)           :: gamma_a
    real(WP)           :: gamma_b
    real(WP)           :: alpha_gamma, alpha_pi
    integer            :: k
    type(gamma_func_t) :: gf
    type(point_t)      :: pt_a
    type(point_t)      :: pt_b

    if (present(a_g)) then
            alpha_gamma = a_g
    else
            alpha_gamma = 1._WP
    endif
    if (present(a_p)) then
            alpha_pi = a_p
    else
            alpha_pi = 1._WP
    endif

    $ASSERT_DEBUG(cx%point_i() == gr%pt_i(),Context and grid are not conformable)
    $ASSERT_DEBUG(cx%point_o() == gr%pt_o(),Context and grid are not conformable)

    ! Find the cell index and abscissa (in grid gr) of the inner
    ! turning point, where the local solution for state st first
    ! becomes propagative

    k_turn = gr%n_k
    x_turn = HUGE(0._WP)

    gamma_b = gamma_(cx, cx%point_i(), st, alpha_gamma, alpha_pi)

    if (gamma_b <= 0._WP) then

       ! Inner point is already propagative

       associate (pt_i => cx%point_i())
         k_turn = 1
         x_turn = pt_i%x
       end associate

    else

       turn_loop : do k = 1, gr%n_k-1

          ! Check for a sign change in gamma

          gamma_a = gamma_b
          gamma_b = gamma_(cx, gr%pt(k+1), st, alpha_gamma, alpha_pi)

          if (gamma_a > 0._WP .AND. gamma_b <= 0._WP) then

             k_turn = k
             
             pt_a = gr%pt(k)
             pt_b = gr%pt(k+1)
             
             if (pt_a%s == pt_b%s) then

                if (abs(gamma_a) < EPSILON(0._WP)*abs(gamma_b)) then

                   x_turn = pt_a%x

                elseif (abs(gamma_b) < EPSILON(0._WP)*abs(gamma_a)) then

                   x_turn = pt_b%x

                else
                   
                   gf%cx => cx
                   gf%s = pt_a%s
                   gf%st => st
                   gf%a_g = alpha_gamma
                   gf%a_p = alpha_pi

                   x_turn = gf%root(pt_a%x, pt_b%x, 0._WP)

                endif

             else

                x_turn = pt_a%x

             end if

             exit turn_loop

          endif

       end do turn_loop

    endif

    ! Finish

    return

  end subroutine find_turn

  !****

  function gamma_ (cx, pt, st, a_g, a_p) result (gamma)

    type(context_t), intent(in)  :: cx
    type(point_t), intent(in)    :: pt
    class(r_state_t), intent(in) :: st
    real(WP), intent(in),optional:: a_g, a_p
    real(WP)                     :: gamma

    real(WP) :: V
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: Gamma_1
    real(WP) :: Omega_rot
    real(WP) :: omega_c
    real(WP) :: lambda
    real(WP) :: g_4
    real(WP) :: g_2
    real(WP) :: g_0

    real(WP) :: alpha_gamma, alpha_pi

    if (present(a_g)) then
            alpha_gamma = a_g
    else
            alpha_gamma = 1._WP
    end if
    if (present(a_p)) then
            alpha_pi = a_p
    else
            alpha_pi = 1._WP
    end if

    ! Calculate the propagation discriminant gamma (< 0 : propagation,
    ! > 0 : evanescence)

    associate (ml => cx%model())

      if (ml%is_vacuum(pt)) then

         gamma = HUGE(0._WP)

      else

         V = ml%coeff(I_V_2, pt)*pt%x**2
         As = ml%coeff(I_As, pt)
         U = ml%coeff(I_U, pt)
         c_1 = ml%coeff(I_C_1, pt)
         Gamma_1 = ml%coeff(I_GAMMA_1, pt)

         Omega_rot = cx%Omega_rot(pt)

         omega_c = cx%omega_c(Omega_rot, st)

         lambda = cx%lambda(Omega_rot, st)

         g_4 = -4._WP*V/Gamma_1*c_1 * alpha_gamma
         g_2 = (As - V/Gamma_1 - U + 4._WP)**2 + 4._WP*V/Gamma_1*As * alpha_gamma * alpha_pi + 4._WP*lambda
         g_0 = -4._WP*lambda*As/c_1 * alpha_pi

         if (g_0 /= 0._WP) then
            gamma = (g_4*omega_c**4 + g_2*omega_c**2 + g_0)/omega_c**2
         else
            gamma = g_4*omega_c**2 + g_2
         endif

      endif

    end associate

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

    gamma = gamma_(this%cx, pt, this%st, this%a_g, this%a_p)

    ! Finish

    return

  end function eval_c_

end module gyre_grid_util
