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

  subroutine eval_tide (ml, process_wave, td_p, os_p, nm_p, gr_p)

    class(model_t), pointer, intent(in) :: ml
    interface
       subroutine process_wave (wv)
         use gyre_wave
         type(wave_t), intent(in) :: wv
       end subroutine process_wave
    end interface
    type(tide_par_t), intent(in)        :: td_p
    type(osc_par_t), intent(in)         :: os_p
    type(num_par_t), intent(in)         :: nm_p
    type(grid_par_t), intent(in)        :: gr_p

    real(WP), allocatable     :: omega(:)
    type(context_t), pointer  :: cx => null()
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

    ! Set up the inertial-frame forcing frequencies array

    omega = [(k*td_p%Omega_orb, k=-td_p%k_max,td_p%k_max)]

    ! Allocate other stuff

    allocate(cx)

    ! Loop over l and m

    l_loop : do l = 2, td_p%l_max

       m_loop : do m = -l, l

          ! Set up the mode parameters

          md_p = mode_par_t(0, l=l, m=m, &
                            n_pg_min=-HUGE(0), n_pg_max=HUGE(0), &
                            rossby=.FALSE., tag='')
          
          ! Create the full grid

          gr = grid_t(ml, omega, gr_p, md_p, os_p)
            
          ! Set up the context

          cx = context_t(ml, gr%pt_i(), gr%pt_o(), md_p, os_p)

          ! Create the bvp_t 

          bp = nad_bvp_t(cx, gr, md_p, nm_p, os_p)

          ! Loop over k

          k_loop : do k = -td_p%k_max, td_p%k_max

             ! Calculate the tidal potential coefficient

             Upsilon_lmk = -td_p%eps_T*td_p%R_a**(l-2)*beta_lm(l, m)*X_lmk(td_p%ec, -(l+1),-m,-k)/(4._WP*PI)

             ! Set up the inhomogeneous boundary conditions

             w_i = 0._WP
         
             w_o = 0._WP
             w_o(2) = (2*l+1)*Upsilon_lmk
         
             ! Solve for the wave function

             st = c_state_t(CMPLX(omega(k), KIND=WP), omega(k))

             wv = wave_t(bp, st, w_i, w_o)

             ! Process it

             call process_wave(wv)

          end do k_loop

       end do m_loop

    end do l_loop

    ! Finish

    return

  end subroutine eval_tide

end module gyre_tide
