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
  use gyre_osc_par
  use gyre_out_par
  use gyre_rot_par
  use gyre_sad_bvp
  use gyre_scan
  use gyre_scan_par
  use gyre_state
  use gyre_summary
  use gyre_tide_util
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
  type(context_t), pointer       :: cx => null()
  type(summary_t)                :: sm_ad
  type(summary_t)                :: sm_nad
  type(detail_t)                 :: dt_ad
  type(detail_t)                 :: dt_nad
  real(WP), allocatable          :: omega(:)
  real(WP), allocatable          :: Omega_orb(:)
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
  call read_out_par(unit, 'ad', ot_p_ad)
  call read_out_par(unit, 'nad', ot_p_nad)

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

     ! Set up the context

     cx = context_t(ml, gr_p_sel, md_p(i), os_p_sel, rt_p_sel)

     ! Set up the frequency arrays

     select case (fr_p_sel%force_type)
     case ('FIXED')
        call build_scan(cx, md_p(i), os_p_sel, sc_p_sel, omega)
     case ('BINARY')
        call build_scan(cx, md_p(i), os_p_sel, sc_p_sel, Omega_orb)
        omega = Omega_orb*fr_p_sel%k
     case default
        $ABORT(Invalid force_type)
     end select

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

        gr = grid_t(cx, omega, gr_p_sel, os_p_sel)

     end if

     ! Calculate wavefunctions

     if (os_p_sel%adiabatic) then

        if (md_p(i)%static) then
           allocate(bp_ad, SOURCE=sad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))
        else
           allocate(bp_ad, SOURCE=ad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))
        endif
           
        call eval_wave_ad(omega, Omega_orb)

        deallocate(bp_ad)

     endif

     if (os_p_sel%nonadiabatic) then

        if (md_p(i)%static) then
           $ABORT(Static nonadiabatic modes not currently implemented)
        else
           allocate(bp_nad, SOURCE=nad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))
        endif
           
        call eval_wave_nad(omega, Omega_orb)

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

  subroutine eval_wave_ad (omega, Omega_orb)

    real(WP), intent(in) :: omega(:)
    real(WP), intent(in) :: Omega_orb(:)

    integer         :: j
    real(WP)        :: v_i(bp_ad%n_i)
    real(WP)        :: v_o(bp_ad%n_o)
    type(r_state_t) :: st
    type(wave_t)    :: wv

    ! Scan over frequencies

    omega_loop : do j = 1, SIZE(omega)

       j_ad = j_ad + 1

       ! Set up the inhomogeneous boundary terms

       v_i = 0._WP
       v_o = 0._WP
       
       select type (bp_ad)
       type is (sad_bvp_t)
          v_o(1) = Phi_force(omega(j), Omega_orb(j))
       type is (ad_bvp_t)
          v_o(2) = Phi_force(omega(j), Omega_orb(j))
       class default
          $ABORT(Invalid bp_ad class)
       end select
         
       ! Solve for the wave function

       st = r_state_t(omega(j))
       
       select type (bp_ad)
       type is (sad_bvp_t)
          wv = wave_t(bp_ad, st, v_i, v_o, j_ad)
       type is (ad_bvp_t)
          wv = wave_t(bp_ad, st, v_i, v_o, j_ad)
       class default
          $ABORT(Invalid bp_ad class)
       end select

       ! Cache/write the mode
    
       call sm_ad%cache(wv)
       call dt_ad%write(wv)

    end do omega_loop

    ! Finish

    return

  end subroutine eval_wave_ad

  !****

  subroutine eval_wave_nad (omega, Omega_orb)

    real(WP), intent(in) :: omega(:)
    real(WP), intent(in) :: Omega_orb(:)

    integer         :: j
    complex(WP)     :: v_i(bp_nad%n_i)
    complex(WP)     :: v_o(bp_nad%n_o)
    type(c_state_t) :: st
    type(wave_t)    :: wv

    ! Scan over frequencies

    omega_loop : do j = 1, SIZE(omega)

       j_nad = j_nad + 1

       ! Set up the inhomogeneous boundary terms

       v_i = 0._WP
         
       v_o = 0._WP
       v_o(2) = Phi_force(omega(j), Omega_orb(j))
         
       ! Solve for the wave function

       st = c_state_t(CMPLX(omega(j), KIND=WP))
       
       select type (bp_nad)
       type is (nad_bvp_t)
          wv = wave_t(bp_nad, st, v_i, v_o, j_nad)
       class default
          $ABORT(Invalid bp_nad class)
       end select

       ! Cache/write the mode
    
       call sm_nad%cache(wv)
       call dt_nad%write(wv)

    end do omega_loop

    ! Finish

    return

  end subroutine eval_wave_nad

  !****

  function Phi_force (omega, Omega_orb)

    real(WP), intent(in) :: omega
    real(WP), intent(in) :: Omega_orb
    real(WP)             :: Phi_force

    real(WP) :: M_pri
    real(WP) :: M_sec
    real(WP) :: R_pri
    real(WP) :: P
    real(WP) :: a
    real(WP) :: eps_tide

    ! Evaluate the surface forcing potential at forcing frequency
    ! omega (in units of GM/R)

    select case (fr_p_sel%force_type)

    case ('FIXED')

       Phi_force = -(2*md_p(i)%l+1)*fr_p_sel%Phi

    case ('BINARY')

       select type (ml)
       class is (evol_model_t)
          M_pri = ml%M_star
          M_sec = fr_p_sel%q*M_pri
          R_pri = ml%R_star
       class default
          $ABORT(Invalid model class)
       end select

       P = 1._WP/(Omega_orb*freq_scale('HZ', cx, md_p(i), os_p_sel))

       a = (G_GRAVITY*(M_pri + M_sec)*P**2/(4.*PI**2))**(1._WP/3._WP)

       eps_tide = (R_pri/a)**3*(M_sec/M_pri)

       Phi_force = -(2*md_p(i)%l+1)*eps_tide/sqrt(4._WP*PI)*tidal_c(R_pri/a, fr_p_sel%e, md_p(i)%l, md_p(i)%m, fr_p_sel%k)

    case default

       $ABORT(Invalid force_type)

    end select

    ! Finish

    return

  end function Phi_force

end program gyre_force
