! Module   : gyre_grid_util
! Purpose  : grid utilities
!
! Copyright 2016-2017 Rich Townsend
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

  use gyre_grid
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rot
  use gyre_rot_factory
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: find_turn

  ! Procedures

contains

  subroutine find_turn (ml, gr, omega, md_p, os_p, k_turn)

    class(model_t), pointer, intent(in) :: ml
    type(grid_t), intent(in)            :: gr
    real(WP), intent(in)                :: omega
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    integer, intent(out)                :: k_turn

    class(r_rot_t), allocatable :: rt
    integer                     :: k
    real(WP)                    :: V_g
    real(WP)                    :: As
    real(WP)                    :: U
    real(WP)                    :: c_1
    real(WP)                    :: lambda
    real(WP)                    :: g_0
    real(WP)                    :: g_2
    real(WP)                    :: g_4
    real(WP)                    :: gamma

    ! Find the cell index of the inner turning point, where the local
    ! solution at frequency omega first becomes propagative

    ! Set up rt for every point in the grid

    allocate(rt, SOURCE=r_rot_t(ml, gr%pt(1), md_p, os_p))

    call rt%stencil(gr%pt)

    ! Scan through the grid until the propagation discriminant becomes positive

    k_turn = gr%n_k

    gamma_loop : do k = 1, gr%n_k-1

       associate (pt => gr%pt(k))

         if (ml%is_vacuum(pt)) then

            gamma = HUGE(0._WP)

         else

            V_g = ml%coeff(I_V_2, pt)*pt%x**2/ml%coeff(I_GAMMA_1, pt)
            As = ml%coeff(I_As, pt)
            U = ml%coeff(I_U, pt)
            c_1 = ml%coeff(I_C_1, pt)

            lambda = rt%lambda(k, omega)

            g_4 = -4._WP*V_g*c_1
            g_2 = (As - V_g - U + 4._WP)**2 + 4._WP*V_g*As + 4._WP*lambda
            g_0 = -4._WP*lambda*As/c_1

            gamma = (g_4*omega**4 + g_2*omega**2 + g_0)/omega**2

         end if

       end associate

       if (gamma <= 0._WP) then
          k_turn = k
          exit gamma_loop
       endif

    end do gamma_loop
       
    ! Finish

    return

  end subroutine find_turn

end module gyre_grid_util
