! Program  : gyre_tide
! Purpose  : tidal response evaluation
!
! Copyright 2018-2020 Rich Townsend & The GYRE Team
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
  use gyre_rot_par
  use gyre_sad_bvp
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

  subroutine eval_tide (ml, gr, process_wave, os_p, rt_p, nm_p, or_p, td_p)

    class(model_t), pointer, intent(in) :: ml
    interface
       subroutine process_wave (wv, k)
         use gyre_wave
         type(wave_t), intent(in) :: wv
         integer, intent(in)      :: k
       end subroutine process_wave
    end interface
    type(osc_par_t), intent(in)         :: os_p
    type(rot_par_t), intent(in)         :: rt_p
    type(num_par_t), intent(in)         :: nm_p
    type(orbit_par_t), intent(in)       :: or_p(:)
    type(tide_par_t), intent(in)        :: td_p

    type(grid_t)                   :: ml_gr
    real(WP)                       :: Omega_orb
    real(WP)                       :: R_a
    real(WP)                       :: eps_tide
    type(rot_par_t)                :: rt_p_
    real(WP)                       :: Omega_sync
    integer                        :: k_min
    integer                        :: k_max
    real(WP), allocatable          :: omega(:)
    integer                        :: l_max
    type(mode_par_t), allocatable  :: md_p(:,:)
    type(context_t), pointer       :: cx(:,:) => null()
    integer, allocatable           :: tide_type(:,:,:)
    type(grid_spec_t), allocatable :: gs(:,:)
    integer                        :: l
    integer                        :: m
    integer                        :: k
    type(grid_t)                   :: gr
    type(nad_bvp_t)                :: bp_nad
    type(sad_bvp_t)                :: bp_sad
    real(WP)                       :: c
    complex(WP)                    :: w_i_nad(3)
    complex(WP)                    :: w_o_nad(3)
    real(WP)                       :: w_i_sad(1)
    real(WP)                       :: w_o_sad(1)
    type(c_state_t)                :: st_nad
    type(r_state_t)                :: st_sad
    type(wave_t)                   :: wv

    integer :: c_beg, c_end, c_rate
    integer :: c_c, c_grid, c_bvp, c_solve, c_proc

    c_c = 0
    c_grid = 0
    c_bvp = 0
    c_solve = 0
    c_proc = 0

    ! Set up index ranges

    l_max = td_p%l_max

    m_min = MAX(-l_max, td_p%m_min)
    m_max = MIN( l_max, td_p%m_max)

    if (td_p%reduce_order) then
       k_min = 0
    else
       k_min = -td_p%k_max
    endif

    k_max = td_p%k_max

    ! Allocate arrays for index-dependent quantities

    allocate(c(2:l_max,m_min:m_max,k_min:k_max))

    allocate(cx(2:l_max,-l_max:l_max))

    allocate(bp_sad(2:l_max,l_max:l_max))
    allocate(bp_nad(2:l_max,l_max:l_max))

    ! Initialize contexts and bvp's

    init_l_loop : do l = 2, l_max
       init_m_loop : do m = MAX(-l, m_min), MIN(l, m_max)

          ! Create the mode_par_t

          md_p = mode_par_t(l=l, m=m)

          ! Set up the context_t

          cx(l,m) = context_t(ml, gr_p, md_p, os_p, rt_p_)

          ! Create the bvp_t

          bp_nad(l,m) = nad_bvp_t(cx(l,m), gr, md_p, nm_p, os_p)
          bp_sad(l,m) = sad_bvp_t(cx(l,m), gr, md_p, nm_p, os_p)

       end do init_m_loop
    end do init_l_loop

          ! Classify the tide for eack k

          classify_loop : do k = k_min, k_max
             tide_type(l,m,k) = classify_tide_(cx(l,m), ml_gr, omega(k), td_p%omega_static)
             if (check_log_level('DEBUG')) then
                write(OUTPUT_UNIT, *) 'tide type:',l,m,k,tide_type(l,m,k),tidal_c(R_a, or_p%e, l, m, k)

             endif
          end do classify_loop

          ! Create the grid_spec_t

          gs(l,m) = grid_spec_t(cx(l,m), PACK(omega, MASK=(tide_type(l,m,:) == DYNAMIC_TIDE)))
             
          
    

    ! Loop over 

    ! Extract the model grid (used for tide classification)

    ml_gr = ml%grid()

    ! Calculate the dimensionless orbital frequency

    Omega_orb = or_p%Omega_orb/freq_scale(or_p%Omega_orb_units, ml)

    ! Calculate the orbital separation and tidal strength

    R_a = (Omega_orb**2/(1._WP + or_p%q))**(1._WP/3._WP)

    eps_tide = (R_a)**3*or_p%q

    ! If necessary, set up the synchronous rotation rate (based on
    ! eqn. 42 of Hut 1981)

    rt_p_ = rt_p

    if (or_p%sync_rot) then

       associate (e => or_p%e)
         Omega_sync = Omega_orb*(1 + 15*e**2/2 + 45*e**4/8 + 5*e**6/16)/ &
                                ((1 + 3*e**2 + 3*e**4/8)*sqrt(1 - e**2)**3)
       end associate

       rt_p_%Omega_rot_source = 'UNIFORM'
       rt_p_%Omega_rot = or_p%sync_fraction*Omega_sync
       rt_p_%Omega_rot_units = 'NONE'

    endif

    ! Set up the inertial-frame forcing frequencies array

    k_max = td_p%k_max

    if (td_p%combine_k) then
       k_min = 0
    else
       k_min = -k_max
    endif
    
    allocate(omega(k_min:k_max))

    do k = k_min, k_max
       omega(k) = -k*Omega_orb
    end do

    ! Set up contexts and tide types

    l_max = td_p%l_max

    allocate(md_p(2:l_max,-l_max:l_max))
    allocate(cx(2:l_max,-l_max:l_max))
    allocate(tide_type(2:l_max,-l_max:l_max,k_min:k_max))
    allocate(gs(2:l_max,-l_max:l_max))

    tide_type = NO_TIDE

    cx_l_loop : do l = 2, l_max
       cx_m_loop : do m = MAX(-l, td_p%m_min), MIN(l, td_p%m_max)

          ! Create the mode_par_t

          md_p(l,m) = mode_par_t(l=l, m=m)

          ! Set up the context_t

          cx(l,m) = context_t(ml, gr_p, md_p(l,m), os_p, rt_p_)

          ! Classify the tide for eack k

          classify_loop : do k = k_min, k_max
             tide_type(l,m,k) = classify_tide_(cx(l,m), ml_gr, omega(k), td_p%omega_static)
             if (check_log_level('DEBUG')) then
                write(OUTPUT_UNIT, *) 'tide type:',l,m,k,tide_type(l,m,k),tidal_c(R_a, or_p%e, l, m, k)

             endif
          end do classify_loop

          ! Create the grid_spec_t

          gs(l,m) = grid_spec_t(cx(l,m), PACK(omega, MASK=(tide_type(l,m,:) == DYNAMIC_TIDE)))
             
       end do cx_m_loop
    end do cx_l_loop

    ! Use the grid_specs to create the grid. Filter out static and
    ! mixed tides

    call system_clock(c_beg, COUNT_RATE=c_rate)

    gr = grid_t(PACK(gs, MASK=ANY(tide_type == DYNAMIC_TIDE, DIM=3)), gr_p)

    call system_clock(c_end)

    c_grid = c_grid + (c_end - c_beg)

    ! Loop over l, m and k, calculating the tidal contribution

    l_loop : do l = 2, td_p%l_max

       m_loop : do m = MAX(-l, td_p%m_min), MIN(l, td_p%m_max)

          ! Create the bvp_t's

          call system_clock(c_beg)

          bp_nad = nad_bvp_t(cx(l,m), gr, md_p(l,m), nm_p, os_p)
          bp_sad = sad_bvp_t(cx(l,m), gr, md_p(l,m), nm_p, os_p)

          call system_clock(c_end)
          
          c_bvp = c_bvp + (c_end - c_beg)

          k_loop : do k = k_min, k_max

             ! Calculate the tidal potential coefficient

             call system_clock(c_beg)
                
             c = tidal_c(R_a, or_p%e, l, m, k)

             call system_clock(c_end)

             c_c = c_c + (c_end - c_beg)

             if (c /= 0._WP) then

                select case (tide_type(l,m,k))

                case (DYNAMIC_TIDE)

                   ! Set up the inhomogeneous boundary conditions

                   w_i_nad = 0._WP
         
                   w_o_nad = 0._WP
                   w_o_nad(2) = -(2*l+1)*eps_tide/sqrt(4._WP*PI)*c
         
                   ! Solve for the wave function

                   call system_clock(c_beg)

                   st_nad = c_state_t(CMPLX(omega(k), KIND=WP), omega(k))

                   wv = wave_t(bp_nad, st_nad, w_i_nad, w_o_nad, 0)

                   call system_clock(c_end)

                   c_solve = c_solve + (c_end - c_beg)

                   ! Process it

                   call system_clock(c_beg)

                   call process_wave(wv, k)

                   call system_clock(c_end)

                   c_proc = c_proc + (c_end - c_beg)

                case (STATIC_TIDE)

                   ! Set up the inhomogeneous boundary conditions

                   w_i_sad = 0._WP
         
                   w_o_sad = 0._WP
                   w_o_sad(1) = -(2*l+1)*eps_tide/sqrt(4._WP*PI)*c
         
                   ! Solve for the wave function

                   call system_clock(c_beg)

                   st_sad = r_state_t(omega(k))

                   wv = wave_t(bp_sad, st_sad, w_i_sad, w_o_sad, 0)

                   call system_clock(c_end)

                   c_solve = c_solve + (c_end - c_beg)

                   ! Process it

                   call system_clock(c_beg)

                   call process_wave(wv, k)

                   call system_clock(c_end)

                   c_proc = c_proc + (c_end - c_beg)

                case (MIXED_TIDE)

                   ! Ignore mixed (dynamic/static) tides

                   if (check_log_level('WARN')) then
                      write(OUTPUT_UNIT, 110) 'Ignoring mixed tide for l,m,k=', l, m, k
110                   format(A,1X,I0,1X,I0,1X,I0)
                   endif

                case default

                   $ABORT(Invalid tide_type)

                end select

             end if

          end do k_loop

       end do m_loop

    end do l_loop

    ! Report timing

    if (check_log_level('DEBUG')) then
       write(OUTPUT_UNIT, *) 'Time for c:', REAL(c_c)/c_rate
       write(OUTPUT_UNIT, *) 'Time for grids:', REAL(c_grid)/c_rate
       write(OUTPUT_UNIT, *) 'Time for bpv:', REAL(c_bvp)/c_rate
       write(OUTPUT_UNIT, *) 'Time for solve:', REAL(c_solve)/c_rate
       write(OUTPUT_UNIT, *) 'Time for proc:', REAL(c_proc)/c_rate
    end if

    ! Clean up

    deallocate(cx)

    ! Finish

    return

  contains

    function classify_tide_ (cx, gr, omega, omega_static) result (tide_type)

      type(context_t), intent(in) :: cx
      type(grid_t), intent(in)    :: gr
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
         Omega_rot = cx%Omega_rot(gr%pt(k))
         omega_c(k) = cx%omega_c(Omega_rot, st)
      end do

      if (ALL(abs(omega_c) > omega_static)) then
         tide_type = DYNAMIC_TIDE
      elseif (ALL(abs(omega_c) <= omega_static)) then
         tide_type = STATIC_TIDE
      else
         tide_type = MIXED_TIDE
      endif

      ! Finish

      return

    end function classify_tide_

  end subroutine eval_tide

end module gyre_tide
