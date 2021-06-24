! Program  : gyre_force
! Purpose  : forced oscillation code
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

program gyre_force

  ! Uses

  use core_kinds, only : WP
  use core_parallel
  use core_system

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_constants
  use gyre_detail
  use gyre_evol_model
  use gyre_context
  use gyre_force_par
  use gyre_freq_context
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_math
  use gyre_mode_par
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_orbit_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_resp
  use gyre_rot_par
  use gyre_sad_bvp
  use gyre_scan
  use gyre_scan_par
  use gyre_state
  use gyre_summary
  use gyre_tide_util
  use gyre_tnad_bvp
  use gyre_util
  use gyre_version
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable      :: filename
  integer                        :: unit
  type(model_par_t)              :: ml_p
  type(mode_par_t), allocatable  :: md_p(:)
  type(osc_par_t), allocatable   :: os_p(:)
  type(rot_par_t), allocatable   :: rt_p(:)
  type(num_par_t), allocatable   :: nm_p(:)
  type(grid_par_t), allocatable  :: gr_p(:)
  type(scan_par_t), allocatable  :: sc_p(:)
  type(force_par_t), allocatable :: fr_p(:)
  type(orbit_par_t), allocatable :: or_p(:)
  type(out_par_t)                :: ot_p_ad
  type(out_par_t)                :: ot_p_nad
  class(model_t), pointer        :: ml => null()
  integer                        :: i
  type(osc_par_t)                :: os_p_sel
  type(rot_par_t)                :: rt_p_sel
  type(num_par_t)                :: nm_p_sel
  type(grid_par_t)               :: gr_p_sel
  type(scan_par_t), allocatable  :: sc_p_sel(:)
  type(force_par_t)              :: fr_p_sel
  type(orbit_par_t)              :: or_p_sel
  type(context_t), pointer       :: cx => null()
  type(summary_t)                :: sm_ad
  type(summary_t)                :: sm_nad
  type(detail_t)                 :: dt_ad
  type(detail_t)                 :: dt_nad
  real(WP), allocatable          :: omega(:)
  real(WP), allocatable          :: Omega_orb(:)
  real(WP), allocatable          :: Omega_rot(:)
  type(grid_t)                   :: gr
  class(r_bvp_t), allocatable    :: bp_ad
  class(c_bvp_t), allocatable    :: bp_nad
  integer                        :: j_ad
  integer                        :: j_nad

  ! Read command-line arguments

  $ASSERT(n_arg() == 1,Syntax: gyre_force <filename>)

  call get_arg(1, filename)

  ! Initialize

  call init_parallel()
  call init_math()

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_force ['//VERSION//']', '=')
100  format(A)

     write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
     write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110  format(A,1X,A)

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
120  format(A,1X,I0)
     
     write(OUTPUT_UNIT, 110) 'Input filename   :', filename

  endif

  ! Read the namelist file

  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  call read_constants(unit)

  call read_model_par(unit, ml_p)
  call read_mode_par(unit, md_p)
  call read_osc_par(unit, os_p)
  call read_rot_par(unit, rt_p)
  call read_num_par(unit, nm_p)
  call read_grid_par(unit, gr_p)
  call read_scan_par(unit, sc_p)
  call read_force_par(unit, fr_p)
  call read_orbit_par(unit, or_p)
  call read_out_par(unit, 'ad', ot_p_ad)
  call read_out_par(unit, 'nad', ot_p_nad)

  ! Check that GYRE_DIR is set

  $ASSERT(GYRE_DIR /= '',The GYRE_DIR environment variable is not set)

  ! Initialize the model

  if (check_log_level('INFO')) then
     write(OUTPUT_UNIT, 100) form_header('Model Init', '=')
  endif

  ml => model_t(ml_p)

  ! Allocate the context (will be initialized later on)

  allocate(cx)

  ! Initialize the summary and detail outputters

  sm_ad = summary_t(ot_p_ad)
  sm_nad = summary_t(ot_p_nad)

  dt_ad = detail_t(ot_p_ad)
  dt_nad = detail_t(ot_p_nad)

  ! Loop through md_p

  j_ad = 0
  j_nad = 0

  md_p_loop : do i = 1, SIZE(md_p)

     if (check_log_level('INFO')) then

        write(OUTPUT_UNIT, 100) form_header('Mode Search', '=')

        write(OUTPUT_UNIT, 100) 'Mode parameters'

        write(OUTPUT_UNIT, 130) 'l :', md_p(i)%l
        write(OUTPUT_UNIT, 130) 'm :', md_p(i)%m
130     format(3X,A,1X,I0)

        write(OUTPUT_UNIT, *)

     endif

     ! Select parameters according to tags

     call select_par(os_p, md_p(i)%tag, os_p_sel)
     call select_par(rt_p, md_p(i)%tag, rt_p_sel)
     call select_par(nm_p, md_p(i)%tag, nm_p_sel)
     call select_par(gr_p, md_p(i)%tag, gr_p_sel)
     call select_par(sc_p, md_p(i)%tag, sc_p_sel)
     call select_par(fr_p, md_p(i)%tag, fr_p_sel)
     call select_par(or_p, md_p(i)%tag, or_p_sel)

     ! Set up the context

     cx = context_t(ml, gr_p_sel, md_p(i), os_p_sel, rt_p_sel)

     ! Set up the frequency arrays

     call set_freqs(omega, Omega_orb, Omega_rot)

     if (SIZE(omega) == 0) then

        if (check_log_level('INFO')) then
           write(OUTPUT_UNIT, 100) 'Scan is empty, skipping mode...'
        endif
           
        cycle md_p_loop

     endif

     ! Set up the grid

     if (md_p(i)%static) then

        gr = ml%grid()

     else

        if (gr_p_sel%file /= '') then
           gr = grid_from_file(gr_p_sel%file)
        else
           gr = grid_t(cx, omega, gr_p_sel, os_p_sel)
        endif

     end if

     ! Calculate wavefunctions

     if (os_p_sel%adiabatic) then

        if (md_p(i)%static) then
           allocate(bp_ad, SOURCE=sad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))
        else
           allocate(bp_ad, SOURCE=ad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))
        endif
           
        call eval_wave_ad(omega, Omega_orb, Omega_rot)

        deallocate(bp_ad)

     endif

     if (os_p_sel%nonadiabatic) then

        if (md_p(i)%static) then
           $ABORT(Static nonadiabatic modes not currently implemented)
        elseif (os_p_sel%alpha_trb > 0._WP) then
           allocate(bp_nad, SOURCE=tnad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))
        else
           allocate(bp_nad, SOURCE=nad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))
        endif
           
        call eval_wave_nad(omega, Omega_orb, Omega_rot)

        deallocate(bp_nad)

     end if

  end do md_p_loop

  ! Write the summaries

  call sm_ad%write()
  call sm_nad%write()

  ! Clean up

  deallocate(ml)

  ! Finish

  close(unit)

  call final_parallel()

contains

  subroutine set_freqs(omega, Omega_orb, Omega_rot)

     real(WP), allocatable, intent(out) :: omega(:)
     real(WP), allocatable, intent(out) :: Omega_orb(:)
     real(WP), allocatable, intent(out) :: Omega_rot(:)

     ! Set up the frequency arrays -- omega (the forcing frequency),
     ! Omega_orb (the orbital frequency) and Omega_rot (the rotation
     ! frequency)

     select case (fr_p_sel%force_type)

     case ('FIXED')

        call build_scan(cx, md_p(i), os_p_sel, sc_p_sel, 'REAL', omega)

        ! Omega_orb and Omega_rot not set

     case ('BINARY')

        select case (fr_p_sel%scan_var)

        case ('OMEGA')

           $ASSERT(fr_p_sel%k /= 0,Cannot scan over omega when k == 0)

           call build_scan(cx, md_p(i), os_p_sel, sc_p_sel, 'REAL', omega)

           Omega_orb = omega/fr_p_sel%k

           ! Omega_rot not set

        case ('OMEGA_ORB')

           call build_scan(cx, md_p(i), os_p_sel, sc_p_sel, 'REAL', Omega_orb)

           omega = fr_p_sel%k*Omega_orb

           ! Omega_rot not set

        case ('OMEGA_ROT')

           call build_scan(cx, md_p(i), os_p_sel, sc_p_sel, 'REAL', Omega_rot)

           allocate(Omega_orb(SIZE(Omega_rot)))

           Omega_orb = or_p_sel%Omega_orb/freq_scale(or_p_sel%Omega_orb_units, ml)
           omega = fr_p_sel%k*Omega_orb

        case default

           $ABORT(Invalid scan_var)

        end select

     case default

        $ABORT(Invalid force_type)

     end select

     ! Finish

     return

  end subroutine set_freqs

  !****

  subroutine eval_wave_ad (omega, Omega_orb, Omega_rot)

    real(WP), intent(in)              :: omega(:)
    real(WP), intent(in), allocatable :: Omega_orb(:)
    real(WP), intent(in), allocatable :: Omega_rot(:)

    integer          :: j
    real(WP)         :: v_i(bp_ad%n_i)
    real(WP)         :: v_o(bp_ad%n_o)
    type(r_state_t)  :: st
    type(wave_t)     :: wv
    type(resp_t)     :: rs

    ! Scan over frequencies

    omega_loop : do j = 1, SIZE(omega)

       j_ad = j_ad + 1

       ! If necessary, set the rotation rate

       if (ALLOCATED(Omega_rot)) then
          call cx%set_Omega_rot(Omega_rot(j))
       end if
          
       ! Set up the inhomogeneous boundary terms

       v_i = 0._WP
       v_o = 0._WP
       
       select type (bp_ad)
       type is (sad_bvp_t)
          v_o(1) = Phi_force(j, Omega_orb)
       type is (ad_bvp_t)
          v_o(2) = Phi_force(j, Omega_orb)
       class default
          $ABORT(Invalid bp_ad class)
       end select

       ! Solve for the wave function and response

       st = r_state_t(omega(j))
       
       select type (bp_ad)
       type is (sad_bvp_t)
          wv = wave_t(bp_ad, st, v_i, v_o, j_ad)
       type is (ad_bvp_t)
          wv = wave_t(bp_ad, st, v_i, v_o, j_ad)
       class default
          $ABORT(Invalid bp_ad class)
       end select

       if (ALLOCATED(Omega_orb)) then
          rs = resp_t(wv, or_p_sel, Omega_orb(j), fr_p_sel%k)
       else
          rs = resp_t(wv, or_p_sel, 0._WP, fr_p_sel%k)
       end if

       ! Cache/write the response
    
       call sm_ad%cache(rs)
       call dt_ad%write(rs)

    end do omega_loop

    ! Finish

    return

  end subroutine eval_wave_ad

  !****

  subroutine eval_wave_nad (omega, Omega_orb, Omega_rot)

    real(WP), intent(in)              :: omega(:)
    real(WP), allocatable, intent(in) :: Omega_orb(:)
    real(WP), allocatable, intent(in) :: Omega_rot(:)

    integer         :: j
    complex(WP)     :: v_i(bp_nad%n_i)
    complex(WP)     :: v_o(bp_nad%n_o)
    type(c_state_t) :: st
    type(wave_t)    :: wv
    type(resp_t)    :: rs

    ! Scan over frequencies

    omega_loop : do j = 1, SIZE(omega)

       j_nad = j_nad + 1

       ! If necessary, set the rotation rate

       if (ALLOCATED(Omega_rot)) then
          call cx%set_Omega_rot(Omega_rot(j))
       end if
          
       ! Set up the inhomogeneous boundary terms

       v_i = 0._WP
         
       v_o = 0._WP
       v_o(2) = Phi_force(j, Omega_orb)
         
       ! Solve for the wave function and response

       st = c_state_t(CMPLX(omega(j), KIND=WP))
       
       select type (bp_nad)
       type is (nad_bvp_t)
          wv = wave_t(bp_nad, st, v_i, v_o, j_nad)
       type is (tnad_bvp_t)
          wv = wave_t(bp_nad, st, v_i, v_o, j_nad)
       class default
          $ABORT(Invalid bp_nad class)
       end select

       if (ALLOCATED(Omega_orb)) then
          rs = resp_t(wv, or_p_sel, Omega_orb(j), fr_p_sel%k)
       else
          rs = resp_t(wv, or_p_sel, 0._WP, fr_p_sel%k)
       endif

       ! Cache/write the response
    
       call sm_nad%cache(rs)
       call dt_nad%write(rs)

    end do omega_loop

    ! Finish

    return

  end subroutine eval_wave_nad

  !****

  function Phi_force (j, Omega_orb)

    integer, intent(in)               :: j
    real(WP), allocatable, intent(in) :: Omega_orb(:)
    real(WP)                          :: Phi_force

    real(WP) :: R_a
    real(WP) :: eps_tide

    ! Evaluate the surface forcing potential

    select case (fr_p_sel%force_type)

    case ('FIXED')

       Phi_force = -(2*md_p(i)%l+1)*fr_p_sel%Phi

    case ('BINARY')

       R_a = tidal_R_a(Omega_orb(j), or_p_sel%q)

       eps_tide = R_a**3*or_p_sel%q

       Phi_force = -(2*md_p(i)%l+1)*eps_tide*tidal_c(R_a, or_p_sel%e, md_p(i)%l, md_p(i)%m, fr_p_sel%k)

    case default

       $ABORT(Invalid force_type)

    end select

    ! Finish

    return

  end function Phi_force

  !****

  function grid_from_file (file) result (gr)

     use core_hgroup

     character(*), intent(in) :: file
     type(grid_t)             :: gr

     type(hgroup_t) :: hg
     real(WP), allocatable :: x(:)

     ! Create the grid from the file

     hg = hgroup_t(file, OPEN_FILE)

     call read_dset_alloc(hg, 'x', x)

     call hg%final()

     gr = grid_t(x)

     ! Finish

     return

  end function grid_from_file
  
end program gyre_force
