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
  use gyre_grid_spec
  use gyre_mode_par
  use gyre_model
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
  use gyre_state
  use gyre_tide_par
  use gyre_tide_util
  use gyre_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: NO_TIDE = 0
  integer, parameter :: DYNAMIC_TIDE = 1
  integer, parameter :: STATIC_TIDE = 2
  integer, parameter :: MIXED_TIDE = 3

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

    real(WP)                       :: Omega_orb
    real(WP)                       :: eps_tide
    integer                        :: k_max
    real(WP), allocatable          :: omega(:)
    type(grid_t)                   :: ml_gr
    integer                        :: l_max
    type(mode_par_t), allocatable  :: md_p(:,:)
    type(context_t), pointer       :: cx(:,:) => null()
    integer, allocatable           :: tide_type(:,:,:)
    type(grid_spec_t), allocatable :: gs(:,:)
    integer                        :: l
    integer                        :: m
    integer                        :: k
    type(grid_t)                   :: gr
    type(nad_bvp_t)                :: bp
    real(WP)                       :: Upsilon
    complex(WP)                    :: w_i(3)
    complex(WP)                    :: w_o(3)
    type(c_state_t)                :: st
    type(wave_t)                   :: wv

    integer :: c_beg, c_end, c_rate
    integer :: c_ups, c_grid, c_bvp, c_solve, c_proc

    c_ups = 0
    c_grid = 0
    c_bvp = 0
    c_solve = 0
    c_proc = 0

    ! Calculate the orbital frequency and tidal strength

    Omega_orb = SQRT((1._WP + td_p%q)*td_p%R_a**3)

    eps_tide = (td_p%R_a)**3*td_p%q

    ! Set up the inertial-frame forcing frequencies array

    k_max = td_p%k_max

    allocate(omega(0:k_max))

    do k = 0, k_max
       omega(k) = k*Omega_orb
    end do

    ! Extract the model grid (used for checking co-rotating frequencies)

    ml_gr = ml%grid()

    ! Set up contexts and tide types

    l_max = td_p%l_max

    allocate(md_p(2:l_max,-l_max:l_max))
    allocate(cx(2:l_max,-l_max:l_max))
    allocate(tide_type(2:l_max,-l_max:l_max,0:k_max))
    allocate(gs(2:l_max,-l_max:l_max))

    tide_type = NO_TIDE

    cx_l_loop : do l = 2, l_max
       cx_m_loop : do m = -l, l

          ! Create the mode_par_t

          md_p(l,m) = mode_par_t(0, l=l, m=m, &
                                 n_pg_min=-HUGE(0), n_pg_max=HUGE(0), &
                                 rossby=.FALSE., tag='')

          ! Set up the context_t

          cx(l,m) = context_t(ml, gr_p, md_p(l,m), os_p)

          ! Classify the tide for eack k

          classify_loop : do k = 0, k_max
             tide_type(l,m,k) = classify_tide_(ml, ml_gr, cx(l,m), omega(k), td_p%omega_static)
             print *,'tide type:',l,m,k,tide_type(l,m,k)
          end do classify_loop

          ! Create the grid_spec_t

          gs(l,m) = grid_spec_t(cx(l,m), PACK(omega, MASK=(tide_type(l,m,:) == DYNAMIC_TIDE)))
             
       end do cx_m_loop
    end do cx_l_loop

    ! Use the grid_specs to create the grid

    call system_clock(c_beg, COUNT_RATE=c_rate)

    gr = grid_t(PACK(gs, MASK=ANY(tide_type == DYNAMIC_TIDE, DIM=3)), gr_p)

    call system_clock(c_end)

    c_grid = c_grid + (c_end - c_beg)

    ! Loop over l, m and k, calculating the tidal contribution

    l_loop : do l = 2, td_p%l_max

       m_loop : do m = -l, l

          ! Create the bvp_t

          call system_clock(c_beg)

          bp = nad_bvp_t(cx(l,m), gr, md_p(l,m), nm_p, os_p)

          call system_clock(c_end)
          
          c_bvp = c_bvp + (c_end - c_beg)

          k_loop : do k = 0, td_p%k_max

             select case (tide_type(l,m,k))

             case (DYNAMIC_TIDE)

                ! Dynamic tide: calculate the tidal potential coefficient

                call system_clock(c_beg)

                Upsilon = Upsilon_lmk(td_p, l, m, k)

                call system_clock(c_end)

                c_ups = c_ups + (c_end - c_beg)

                if (Upsilon /= 0._WP) then

                   ! Set up the inhomogeneous boundary conditions

                   w_i = 0._WP
         
                   w_o = 0._WP
                   w_o(2) = (2*l+1)*eps_tide*Upsilon
         
                   ! Solve for the wave function

                   call system_clock(c_beg)

                   st = c_state_t(CMPLX(omega(k), KIND=WP), omega(k))

                   wv = wave_t(bp, st, w_i, w_o)

                   call system_clock(c_end)

                   c_solve = c_solve + (c_end - c_beg)

                   ! Process it

                   call system_clock(c_beg)

                   call process_wave(wv)

                   call system_clock(c_end)

                   c_proc = c_proc + (c_end - c_beg)

                endif

             case (STATIC_TIDE)

                ! Ignore static tides

                if (check_log_level('INFO')) then
                   write(output_unit, 110) 'Info: ignoring static tide for l,m,k=', l, m, k
110                format(A,1X,I0,1X,I0,1X,I0)
                end if

             case (MIXED_TIDE)

                ! Ignore mixed (dynamic/static) tides

                if (check_log_level('WARN')) then
                   write(OUTPUT_UNIT, 110) 'Warning: ignoring mixed tide for l,m,k=', l, m, k
                endif

             case default

                $ABORT(Invalid tide_type)

             end select

          end do k_loop

       end do m_loop

    end do l_loop

    print *,'Time for ups:', REAL(c_ups)/c_rate
    print *,'Time for grids:', REAL(c_grid)/c_rate
    print *,'Time for bpv:', REAL(c_bvp)/c_rate
    print *,'Time for solve:', REAL(c_solve)/c_rate
    print *,'Time for proc:', REAL(c_proc)/c_rate

    ! Clean up

    deallocate(cx)

    ! Finish

    return

  contains

    function classify_tide_ (ml, gr, cx, omega, omega_static) result (tide_type)

      class(model_t), intent(in)  :: ml
      type(grid_t), intent(in)    :: gr
      type(context_t), intent(in) :: cx
      real(WP), intent(in)        :: omega
      real(WP), intent(in)        :: omega_static
      integer                     :: tide_type

      type(r_state_t) :: st
      integer         :: k
      real(WP)        :: Omega_rot
      real(WP)        :: omega_c(gr%n_k)

      ! Check co-rotating frequencies

      st = r_state_t(omega)

      !$OMP PARALLEL DO PRIVATE (Omega_rot)
      do k = 1, gr%n_k
         Omega_rot = ml%coeff(I_OMEGA_ROT, gr%pt(k))
         omega_c(k) = cx%omega_c(Omega_rot, st)
      end do

      if (ALL(ABS(omega_c) > omega_static)) then
         tide_type = DYNAMIC_TIDE
      elseif (ALL(ABS(omega_c) <= omega_static)) then
         tide_type = STATIC_TIDE
      else
         tide_type = MIXED_TIDE
      endif

      ! Finish

      return

    end function classify_tide_

  end subroutine eval_tide

end module gyre_tide
