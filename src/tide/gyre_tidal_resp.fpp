! Program  : gyre_tidal_resp
! Purpose  : routines to evaluate tidal response
!
! Copyright 2018-2022 Rich Townsend & The GYRE Team
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

module gyre_tidal_resp

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_constants
  use gyre_context
  use gyre_ext
  use gyre_freq
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_grid_spec
  use gyre_math
  use gyre_mode_par
  use gyre_model
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_orbit_par
  use gyre_osc_par
  use gyre_point
  use gyre_resp
  use gyre_rot_par
  use gyre_sad_bvp
  use gyre_state
  use gyre_tide_par
  use gyre_tidal_coeff
  use gyre_tidal_context
  use gyre_tnad_bvp
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

  ! Module variables

  integer, save :: id_m = 0

  ! Access specifiers

  private

  public :: eval_resp

  ! Procedures

contains

  subroutine eval_resp (ml, process_resp, gr_p, nm_p, or_p, os_p, rt_p, td_p)

    class(model_t), pointer, intent(in) :: ml
    interface
       subroutine process_resp (rs)
         use gyre_resp
         type(resp_t), intent(in) :: rs
       end subroutine process_resp
    end interface
    type(grid_par_t), intent(in)        :: gr_p
    type(num_par_t), intent(in)         :: nm_p
    type(orbit_par_t), intent(in)       :: or_p
    type(osc_par_t), intent(in)         :: os_p
    type(rot_par_t), intent(in)         :: rt_p
    type(tide_par_t), intent(in)        :: td_p

    real(WP), allocatable         :: omega(:)
    type(grid_t)                  :: gr
    type(mode_par_t), allocatable :: md_p(:,:)
    type(context_t), pointer      :: cx(:,:)
    real(WP), allocatable         :: Phi_T(:,:,:)
    integer, allocatable          :: tide_type(:,:,:)
    type(point_t)                 :: pt_o
    integer                       :: l
    integer                       :: m
    integer                       :: k
    class(c_bvp_t), allocatable   :: bp_nad
    type(sad_bvp_t)               :: bp_sad
    complex(WP)                   :: w_i_nad(3)
    complex(WP)                   :: w_o_nad(3)
    real(WP)                      :: w_i_sad(1)
    real(WP)                      :: w_o_sad(1)
    type(c_state_t)               :: st_nad
    type(r_state_t)               :: st_sad
    type(r_state_t)               :: st_null
    type(wave_t)                  :: wv
    type(resp_t)                  :: rs

    integer(I8) :: c_beg, c_end, c_rate
    integer(I8) :: c_bvp, c_solve, c_proc

    c_bvp = 0
    c_solve = 0
    c_proc = 0

    call system_clock(count_rate=c_rate)

    ! Set up tidal params, contexts, grid, etc

    call setup_tides_(ml, gr_p, or_p, os_p, rt_p, td_p, omega, gr, cx, md_p, Phi_T, tide_type)

    pt_o = gr%pt_o()

    ! Use the grid_specs to create the grid, filtering out non-dynamic/mixed tides

    ! Loop over l, m and k, solving for the tidal response

    id_m = id_m + 1

    l_loop : do l = td_p%l_min, td_p%l_max
       m_loop : do m = MAX(-l, td_p%m_min), MIN(l, td_p%m_max)

          ! Create the bvp_t's

          call system_clock(c_beg)

          if (os_p%alpha_trb > 0._WP) then
             allocate(bp_nad, SOURCE=tnad_bvp_t(cx(l,m), gr, md_p(l,m), nm_p, os_p))
          else
             allocate(bp_nad, SOURCE=nad_bvp_t(cx(l,m), gr, md_p(l,m), nm_p, os_p))
          endif

          bp_sad = sad_bvp_t(cx(l,m), gr, md_p(l,m), nm_p, os_p)

          call system_clock(c_end)
          c_bvp = c_bvp + (c_end - c_beg)

          k_loop : do k = td_p%k_min, td_p%k_max

             ! Evaluate the response

             call system_clock(c_beg)

             select case(tide_type(l,m,k))
             case (DYNAMIC_TIDE,MIXED_TIDE)

                if (ABS(Phi_T(l,m,k)) > td_p%Phi_T_thresh) then

                   w_i_nad = 0._WP

                   w_o_nad = 0._WP
                   w_o_nad(2) = (2*l+1)*(ml%coeff(I_C_1, pt_o)/pt_o%x**2)*Phi_T(l,m,k)

                   st_nad = c_state_t(CMPLX(omega(k), KIND=WP), omega(k))

                   select type (bp_nad)
                   type is (nad_bvp_t)
                      wv = wave_t(bp_nad, st_nad, w_i_nad, w_o_nad, id_m)
                   type is (tnad_bvp_t)
                      wv = wave_t(bp_nad, st_nad, w_i_nad, w_o_nad, id_m)
                   end select

                else

                   st_null = r_state_t(omega(k))

                   wv = null_wave_t_(st_null, cx(l,m), gr, md_p(l,m), nm_p, os_p, id_m, .FALSE.)

                end if

             case (STATIC_TIDE)

                if (ABS(Phi_T(l,m,k)) > td_p%Phi_T_thresh) then

                   w_i_sad = 0._WP

                   w_o_sad = 0._WP
                   w_o_sad(1) = (2*l+1)*(ml%coeff(I_C_1, pt_o)/pt_o%x**2)*Phi_T(l,m,k)

                   st_sad = r_state_t(omega(k))

                   wv = wave_t(bp_sad, st_sad, w_i_sad, w_o_sad, id_m)

                else

                   st_null = r_state_t(omega(k))

                   wv = null_wave_t_(st_null, cx(l,m), gr, md_p(l,m), nm_p, os_p, id_m, .TRUE.)

                end if

             case default

                $ABORT(invalid tide_type)

             end select

             rs = resp_t(wv, or_p, td_p, k)

             call system_clock(c_end)
             c_solve = c_solve + (c_end - c_beg)

             ! Process the response

             call system_clock(c_beg)

             call process_resp(rs)

             call system_clock(c_end)
             c_proc = c_proc + (c_end - c_beg)

          end do k_loop

          ! Clean up

          deallocate(bp_nad)

       end do m_loop
    end do l_loop

    ! Report timing

    if (check_log_level('DEBUG')) then
       write(OUTPUT_UNIT, *) 'Time for bpv:', REAL(c_bvp)/c_rate
       write(OUTPUT_UNIT, *) 'Time for solve:', REAL(c_solve)/c_rate
       write(OUTPUT_UNIT, *) 'Time for proc:', REAL(c_proc)/c_rate
    end if

    ! Clean up

    deallocate(cx)

    ! Finish

    return

  end subroutine eval_resp

  !****

  subroutine setup_tides_ (ml, gr_p, or_p, os_p, rt_p, td_p, &
                           omega, gr, cx, md_p, Phi_T, tide_type)

    class(model_t), pointer, intent(in)        :: ml
    type(grid_par_t), intent(in)               :: gr_p
    type(orbit_par_t), intent(in)              :: or_p
    type(osc_par_t), intent(in)                :: os_p
    type(rot_par_t), intent(in)                :: rt_p
    type(tide_par_t), intent(in)               :: td_p
    real(WP), allocatable, intent(out)         :: omega(:)
    type(grid_t)                               :: gr
    type(context_t), pointer, intent(out)      :: cx(:,:)
    type(mode_par_t), allocatable, intent(out) :: md_p(:,:)
    real(WP), allocatable, intent(out)         :: Phi_T(:,:,:)
    integer, allocatable, intent(out)          :: tide_type(:,:,:)

    type(point_t)                  :: pt_o
    integer                        :: l_min
    integer                        :: l_max
    integer                        :: m_min
    integer                        :: m_max
    integer                        :: k_min
    integer                        :: k_max
    real(WP)                       :: Omega_orb
    type(grid_spec_t), allocatable :: gs(:,:,:)
    logical, allocatable           :: gs_mask(:,:,:)
    integer                        :: l
    integer                        :: m
    integer                        :: k

    ! Set up the scaffold grid (used for tide classification and bc
    ! setting)

    gr = grid_t(ml%grid(), gr_p%x_i, gr_p%x_o)

    pt_o = gr%pt_o()

    ! Initialize arrays

    l_min = td_p%l_min
    l_max = td_p%l_max

    m_min = MAX(-l_max, td_p%m_min)
    m_max = MIN(l_max, td_p%m_max)

    k_min = td_p%k_min
    k_max = td_p%k_max

    allocate(omega(k_min:k_max))

    allocate(md_p(l_min:l_max,m_min:m_max))
    allocate(cx(l_min:l_max,m_min:m_max))

    allocate(Phi_T(l_min:l_max,m_min:m_max,k_min:k_max))
    allocate(tide_type(l_min:l_max,m_min:m_max,k_min:k_max))

    allocate(gs(l_min:l_max,m_min:m_max,k_min:k_max))
    allocate(gs_mask(l_min:l_max,m_min:m_max,k_min:k_max))

    ! Set up the forcing frequency

    Omega_orb = tidal_Omega_orb(ml, or_p)

    omega_loop : do k = k_min, k_max
       omega(k) = k*Omega_orb
    end do omega_loop

    ! Set up other tide params

    gs_mask = .FALSE.

    l_loop : do l = l_min, l_max
       m_loop : do m = MAX(-l, td_p%m_min), MIN(l, td_p%m_max)

          md_p(l,m) = mode_par_t(l=l, m=m)

          cx(l,m) = context_t(ml, gr_p, md_p(l,m), or_p, os_p, rt_p)

          k_loop : do k = k_min, k_max

             Phi_T(l,m,k) = tidal_Phi_T(ml, or_p, pt_o%x, l, m, k)
             tide_type(l,m,k) = classify_tide_(cx(l,m), gr, td_p, omega(k))

             gs_mask(l,m,k) = (tide_type(l,m,k) == DYNAMIC_TIDE .OR. &
                               tide_type(l,m,k) == MIXED_TIDE) .AND. &
                              ABS(Phi_T(l,m,k)) > td_p%Phi_T_thresh

             if (gs_mask(l,m,k)) then
                gs(l,m,k) = grid_spec_t(cx(l,m), [omega(k)])
             end if

             if (check_log_level('DEBUG')) then
                write(OUTPUT_UNIT, *) 'tide type:', l, m, k, tide_type(l,m,k), Phi_T(l,m,k)
             endif

          end do k_loop

       end do m_loop
    end do l_loop

    ! Set up the grid

    if (COUNT(gs_mask) > 0) then
       gr = grid_t(PACK(gs, MASK=gs_mask), gr_p, os_p)
    else
       gr = ml%grid()
    end if

    ! Finish

    return

  end subroutine setup_tides_

  !****

  function classify_tide_ (cx, gr, td_p, omega) result (tide_type)

    type(context_t), intent(in)  :: cx
    type(grid_t), intent(in)     :: gr
    type(tide_par_t), intent(in) :: td_p
    real(WP), intent(in)         :: omega
    integer                      :: tide_type

    type(r_state_t) :: st
    integer         :: j
    real(WP)        :: Omega_rot
    real(WP)        :: omega_c(gr%n)

    ! Check co-rotating frequencies

    st = r_state_t(omega)

    !$OMP PARALLEL DO PRIVATE (Omega_rot)
    do j = 1, gr%n
       Omega_rot = cx%Omega_rot(gr%pt(j))
       omega_c(j) = cx%omega_c(Omega_rot, st)
    end do

    if (ALL(abs(omega_c) > td_p%omega_c_thresh)) then
       tide_type = DYNAMIC_TIDE
    elseif (ALL(abs(omega_c) <= td_p%omega_c_thresh)) then
       tide_type = STATIC_TIDE
    else
       tide_type = MIXED_TIDE
    endif

    ! Finish

    return

  end function classify_tide_

  !****

  function null_wave_t_ (st, cx, gr, md_p, nm_p, os_p, id, static) result (wv)

    type(r_state_t), intent(in)          :: st
    type(context_t), pointer, intent(in) :: cx
    type(grid_t), intent(in)             :: gr
    type(mode_par_t), intent(in)         :: md_p
    type(num_par_t), intent(in)          :: nm_p
    type(osc_par_t), intent(in)          :: os_p
    integer, intent(in)                  :: id
    logical, intent(in)                  :: static
    type(wave_t)                         :: wv

    complex(WP)     :: y_c(6,gr%n)
    type(c_state_t) :: st_c
    type(c_ext_t)   :: discrim

    ! Create a null (zero-amplitudes) wave_t

    ! Set up complex eigenfunctions

    st_c = c_state_t(CMPLX(st%omega, KIND=WP), st%omega)

    y_c = 0._WP

    ! Construct the wave_t

    discrim = c_ext_t(0._WP)

    wv = wave_t(st_c, y_c, discrim, cx, gr, md_p, nm_p, os_p, id, static)

    ! Finish

    return

  end function null_wave_t_

end module gyre_tidal_resp
