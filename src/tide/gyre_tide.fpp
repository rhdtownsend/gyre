! Program  : gyre_tide
! Purpose  : dynamical tide modeling
!
! Copyright 2018 Rich Townsend
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

module gyre_tide

  ! Uses

  use core_kinds

  use gyre_constants
  use gyre_context
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_mode_par
  use gyre_model
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
  use gyre_state
  use gyre_tide_par
  use gyre_tide_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: eval_tide

  ! Procedures

contains

  subroutine eval_tide (ml, process_wave, os_p, nm_p, gr_p, td_p)

    class(model_t), pointer, intent(in) :: ml
    interface
       subroutine process_wave (wv)
         use gyre_wave
         type(wave_t), intent(in) :: wv
       end subroutine process_wave
    end interface
    type(osc_par_t), intent(in)         :: os_p
    type(num_par_t), intent(in)         :: nm_p
    type(grid_par_t), intent(in)        :: gr_p
    type(tide_par_t), intent(in)        :: td_p

    real(WP)                  :: Omega_orb
    real(WP)                  :: eps_tide
    integer                   :: n_cx
    integer                   :: i
    type(context_t), pointer  :: cx(:) => null()
    real(WP), allocatable     :: omega(:)
    type(wave_t)              :: wv
    integer                   :: l
    integer                   :: m
    integer                   :: k
    type(mode_par_t)          :: md_p
    type(grid_t)              :: gr
    type(nad_bvp_t)           :: bp
    real(WP)                  :: Upsilon_lmk
    complex(WP)               :: w_i(3)
    complex(WP)               :: w_o(3)
    type(c_state_t)           :: st

    ! Calculate the orbital frequency and tidal strength

    Omega_orb = SQRT((1._WP + td_p%q)*td_p%R_a**3)

    eps_tide = (td_p%R_a)**3*td_p%q

    ! Set up contexts

    n_cx = td_p%l_max*(td_p%l_max+2) - 3

    allocate(cx(n_cx))

    i = 0

    cx_l_loop : do l = 2, td_p%l_max
       cx_m_loop : do m = -l, l

          i = i + 1
          $ASSERT_DEBUG(i <= n_cx,Context array overflow)
          
          md_p = mode_par_t(i, l=l, m=m, &
                            n_pg_min=-HUGE(0), n_pg_max=HUGE(0), &
                            rossby=.FALSE., tag='')

          cx(i) = context_t(ml, gr_p, md_p, os_p)

       end do cx_m_loop
    end do cx_l_loop

    ! Set up the inertial-frame forcing frequencies array

    omega = [(k*Omega_orb, k=-td_p%k_max,td_p%k_max)]

    ! Create the grid

    gr = grid_t(cx, omega, gr_p)

    ! Loop over l and m, calculating the tidal contribution

    i = 0

    l_loop : do l = 2, td_p%l_max
       m_loop : do m = -l, l

          i = i + 1

          ! Create the bvp_t 

          bp = nad_bvp_t(cx(i), gr, md_p, nm_p, os_p)

          ! Loop over k

          k_loop : do k = -td_p%k_max, td_p%k_max

             ! Calculate the tidal potential coefficient

             Upsilon_lmk = -eps_tide*td_p%R_a**(l-2)*beta_lm(l, m)*X_lmk(td_p%e, -(l+1),-m,-k)/(4._WP*PI)

             if (Upsilon_lmk /= 0._WP) then

                ! Set up the inhomogeneous boundary conditions

                w_i = 0._WP
         
                w_o = 0._WP
                w_o(2) = (2*l+1)*Upsilon_lmk
         
                ! Solve for the wave function

                st = c_state_t(CMPLX(omega(k), KIND=WP), omega(k))

                wv = wave_t(bp, st, w_i, w_o)

                ! Process it

                call process_wave(wv)

             endif

          end do k_loop

       end do m_loop

    end do l_loop

    ! Clean up

    deallocate(cx)

    ! Finish

    return

  end subroutine eval_tide

end module gyre_tide
