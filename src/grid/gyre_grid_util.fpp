! Module   : gyre_grid_util
! Purpose  : grid utilities
!
! Copyright 2016-2021 Rich Townsend & The GYRE Team
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

  use gyre_context
  use gyre_grid
  use gyre_model
  use gyre_num_par
  use gyre_osc_par
  use gyre_point
  use gyre_root
  use gyre_rot
  use gyre_state
  use gyre_status
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: find_turn

  ! Procedures

contains

  subroutine find_turn (cx, gr, st, nm_p, os_p, k_turn, x_turn)

    type(context_t), target, intent(in)  :: cx
    type(grid_t), intent(in)             :: gr
    class(r_state_t), target, intent(in) :: st
    type(num_par_t), intent(in)          :: nm_p
    type(osc_par_t), intent(in)          :: os_p
    integer, intent(out)                 :: k_turn
    real(WP), intent(out)                :: x_turn

    real(WP)      :: gamma_a
    real(WP)      :: gamma_b
    integer       :: k
    type(point_t) :: pt_a
    type(point_t) :: pt_b
    integer       :: status

    $ASSERT_DEBUG(cx%point_i() == gr%pt_i(),Context and grid are not conformable)
    $ASSERT_DEBUG(cx%point_o() == gr%pt_o(),Context and grid are not conformable)

    ! Find the cell index and abscissa (in grid gr) of the inner
    ! turning point, where the local solution for state st first
    ! becomes propagative

    k_turn = gr%n_k
    x_turn = HUGE(0._WP)

    gamma_b = gamma_(cx, cx%point_i(), st, os_p%alpha_gam, os_p%alpha_pi)

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
          gamma_b = gamma_(cx, gr%pt(k+1), st, os_p%alpha_gam, os_p%alpha_pi)

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

                   call solve_root(eval_gamma_, pt_a%x, pt_b%x, 0._WP, nm_p, x_turn, status)

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

  contains

    subroutine eval_gamma_ (x, gamma, status)

      real(WP), intent(in)  :: x
      real(WP), intent(out) :: gamma
      integer, intent(out)  :: status

      ! Evaluate gamma

      gamma = gamma_(cx, point_t(pt_a%s, x), st, os_p%alpha_gam, os_p%alpha_pi)

      status = STATUS_OK

      ! Finish

    end subroutine eval_gamma_

  end subroutine find_turn

  !****

  function gamma_ (cx, pt, st, alpha_gam, alpha_pi) result (gamma)

    type(context_t), intent(in)  :: cx
    type(point_t), intent(in)    :: pt
    class(r_state_t), intent(in) :: st
    real(WP), intent(in)         :: alpha_gam
    real(WP), intent(in)         :: alpha_pi
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

         g_4 = -4._WP*V/Gamma_1*c_1*alpha_gam
         g_2 = (As - V/Gamma_1 - U + 4._WP)**2 + 4._WP*V/Gamma_1*As*alpha_gam*alpha_pi + 4._WP*lambda
         g_0 = -4._WP*lambda*As/c_1*alpha_pi

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

end module gyre_grid_util
