! Program  : gyre
! Purpose  : oscillation code
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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

program gyre

  ! Uses

  use core_kinds, only : WP
  use core_memory
  use core_parallel
  use core_system

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_constants
  use gyre_context
  use gyre_detail
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_math
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_rad_bvp
  use gyre_rot_par
  use gyre_scan
  use gyre_scan_par
  use gyre_search
  use gyre_summary
  use gyre_util
  use gyre_version

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable            :: filename
  integer                              :: unit
  type(model_par_t)                    :: ml_p
  type(mode_par_t), allocatable        :: md_p(:)
  type(osc_par_t), allocatable         :: os_p(:)
  type(rot_par_t), allocatable         :: rt_p(:)
  type(num_par_t), allocatable         :: nm_p(:)
  type(grid_par_t), allocatable        :: gr_p(:)
  type(scan_par_t), allocatable        :: sc_p(:)
  type(out_par_t)                      :: ot_p_ad
  type(out_par_t)                      :: ot_p_nad
  class(model_t), pointer              :: ml => null()
  integer                              :: i
  type(osc_par_t)                      :: os_p_sel
  type(rot_par_t)                      :: rt_p_sel
  type(num_par_t)                      :: nm_p_sel
  type(grid_par_t)                     :: gr_p_sel
  type(scan_par_t), allocatable        :: sc_p_sel(:)
  type(context_t), pointer             :: cx => null()
  type(summary_t)                      :: sm_ad
  type(summary_t)                      :: sm_nad
  type(detail_t)                       :: dt_ad
  type(detail_t)                       :: dt_nad
  real(WP), allocatable                :: omega(:)
  real(WP)                             :: omega_min
  real(WP)                             :: omega_max
  type(grid_t)                         :: gr
  class(r_bvp_t), allocatable          :: bp_ad
  class(c_bvp_t), allocatable          :: bp_nad
  integer                              :: n_ad
  integer                              :: d_ad
  complex(WP), allocatable             :: omega_ad(:)
  integer, allocatable                 :: j_ad(:)

  ! Read command-line arguments

  $ASSERT(n_arg() == 1,Syntax: gyre <filename>)

  call get_arg(1, filename)

  ! Initialize

  call init_parallel()
  call init_math()

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre ['//VERSION//']', '-')
100  format(A)

     if (check_log_level('DEBUG')) then
        write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
        write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110     format(A,1X,A)
     endif

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
120  format(A,1X,I0)

     write(OUTPUT_UNIT, 110) 'Input filename   :', filename

     write(OUTPUT_UNIT, *)

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
  call read_out_par(unit, 'ad', ot_p_ad)
  call read_out_par(unit, 'nad', ot_p_nad)

  ! Initialize the model

  if (check_log_level('INFO')) then
     write(OUTPUT_UNIT, 100) form_header('Model Init', '-')
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

  md_p_loop : do i = 1, SIZE(md_p)

     if (check_log_level('INFO')) then

        write(OUTPUT_UNIT, 100) form_header('Mode Search', '-')

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

     ! Set up the context

     cx = context_t(ml, gr_p_sel, md_p(i), os_p_sel, rt_p_sel)

     ! Set up the frequency array

     call build_scan(cx, md_p(i), os_p_sel, sc_p_sel, omega)

     if (SIZE(omega) < 2) then

        if (check_log_level('INFO')) then
           write(OUTPUT_UNIT, 100) 'Scan is empty, skipping mode...'
        endif

        cycle md_p_loop

     endif

     ! Create the grid

     gr = grid_t(cx, omega, gr_p_sel, os_p_sel)

     ! Set frequency bounds and perform checks

     if (nm_p_sel%restrict_roots) then
        omega_min = MINVAL(omega)
        omega_max = MAXVAL(omega)
     else
        omega_min = -HUGE(0._WP)
        omega_max = HUGE(0._WP)
     endif

     call check_scan(cx, gr, omega, md_p(i), os_p_sel)

     ! Find adiabatic modes

     if (os_p_sel%adiabatic) then
        
        d_ad = 128
        n_ad = 0

        allocate(omega_ad(d_ad))
        allocate(j_ad(d_ad))

        if (md_p(i)%l == 0 .AND. os_p_sel%reduce_order) then
           allocate(bp_ad, SOURCE=rad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))
        else
           allocate(bp_ad, SOURCE=ad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))
        endif

        if (check_log_level('INFO')) then
           write(OUTPUT_UNIT, 100) 'Starting search (adiabatic)'
           write(OUTPUT_UNIT, *)
        endif

        select case (md_p(i)%ad_search)
        case ('SCAN')
           call scan_search(bp_ad, omega, omega_min, omega_max, process_mode_ad, nm_p_sel)
        case default
           $ABORT(Invalid ad_search)
        end select

        deallocate(bp_ad)

     endif

     ! Find non-adiabatic modes

     if (os_p_sel%nonadiabatic) then

        allocate(bp_nad, SOURCE=nad_bvp_t(cx, gr, md_p(i), nm_p_sel, os_p_sel))

        if (check_log_level('INFO')) then
           write(OUTPUT_UNIT, 100) 'Starting search (non-adiabatic)'
           write(OUTPUT_UNIT, *)
        endif

        select case (md_p(i)%nad_search)
        case ('AD')
           $ASSERT(os_p_sel%adiabatic,No adiabatic modes to start from)
           call mode_search(bp_nad, omega_ad(:n_ad), j_ad(:n_ad), omega_min, omega_max, process_mode_nad, nm_p_sel)
        case ('SCAN')
           call scan_search(bp_nad, omega, omega_min, omega_max, process_mode_nad, nm_p_sel)
        case default
           $ABORT(Invalid nad_start)
        end select

        deallocate(bp_nad)

     endif

     ! Deallocate adiabatic data

     if (os_p_sel%adiabatic) then
        deallocate(omega_ad)
        deallocate(j_ad)
     endif

  end do md_p_loop

  ! Write the summaries

  call sm_ad%write()
  call sm_nad%write()

  ! Clean up

  deallocate(cx)

  deallocate(ml)

  ! Finish

  close(unit)

  call final_parallel()

contains

  subroutine process_mode_ad (md, n_iter, chi)

    type(mode_t), intent(in)  :: md
    integer, intent(in)       :: n_iter
    type(r_ext_t), intent(in) :: chi

    ! Process the adiabatic mode

    if (md%n_pg < md_p(i)%n_pg_min .OR. md%n_pg > md_p(i)%n_pg_max) return

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) md%l, md%m, md%n_pg, md%n_p, md%n_g, &
            md%omega, real(chi), n_iter
100    format(1X,I3,1X,I4,1X,I7,1X,I6,1X,I6,1X,E15.8,1X,E15.8,1X,E10.4,1X,I6)
    endif

    ! Store omega and j

    n_ad = n_ad + 1

    if (n_ad > d_ad) then
       d_ad = 2*d_ad
       call reallocate(omega_ad, [d_ad])
       call reallocate(j_ad, [d_ad])
    endif

    omega_ad(n_ad) = md%omega
    j_ad(n_ad) = md%j

    ! Cache/write the mode

    call sm_ad%cache(md)
    call dt_ad%write(md)

    ! Finish

    return

  end subroutine process_mode_ad

  !****

  subroutine process_mode_nad (md, n_iter, chi)

    type(mode_t), intent(in)  :: md
    integer, intent(in)       :: n_iter
    type(r_ext_t), intent(in) :: chi

    ! Process the non-adiabatic mode

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 100) md%l, md%m, md%n_pg, md%n_p, md%n_g, &
            md%omega, real(chi), n_iter
100    format(1X,I3,1X,I4,1X,I7,1X,I6,1X,I6,1X,E15.8,1X,E15.8,1X,E10.4,1X,I6)
    endif

    ! Cache/write the mode

    call sm_nad%cache(md)
    call dt_nad%write(md)

    ! Finish

    return

  end subroutine process_mode_nad

end program gyre
