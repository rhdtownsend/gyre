! Program  : gyre
! Purpose  : oscillation code
!
! Copyright 2013-2017 Rich Townsend
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
  use core_parallel
  use core_system

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_constants
  use gyre_context
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_output
  use gyre_rad_bvp
  use gyre_scan_par
  use gyre_search
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
  type(num_par_t), allocatable         :: nm_p(:)
  type(grid_par_t), allocatable        :: gr_p(:)
  type(scan_par_t), allocatable        :: sc_p(:)
  type(out_par_t)                      :: ot_p_ad
  type(out_par_t)                      :: ot_p_nad
  class(model_t), pointer              :: ml => null()
  integer                              :: i
  type(osc_par_t)                      :: os_p_sel
  type(num_par_t)                      :: nm_p_sel
  type(grid_par_t)                     :: gr_p_sel
  type(scan_par_t), allocatable        :: sc_p_sel(:)
  type(context_t), pointer             :: cx(:) => null()
  real(WP), allocatable                :: omega(:)
  real(WP)                             :: omega_min
  real(WP)                             :: omega_max
  type(grid_t)                         :: gr
  class(r_bvp_t), allocatable          :: bp_ad
  class(c_bvp_t), allocatable          :: bp_nad
  integer                              :: n_md_ad
  integer                              :: d_md_ad
  type(mode_t), allocatable            :: md_ad(:)
  integer                              :: n_md_nad
  integer                              :: d_md_nad
  integer                              :: i_ad_a
  integer                              :: i_ad_b
  type(mode_t), allocatable            :: md_nad(:)

  ! Read command-line arguments

  $ASSERT(n_arg() == 1,Syntax: gyre <filename>)

  call get_arg(1, filename)

  ! Initialize

  call init_parallel()

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

  ! Allocate the contexts array (will be initialized later on)

  allocate(cx(SIZE(md_p)))

  ! Loop through md_p

  d_md_ad = 128
  n_md_ad = 0

  allocate(md_ad(d_md_ad))

  d_md_nad = 128
  n_md_nad = 0

  allocate(md_nad(d_md_nad))

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
     call select_par(nm_p, md_p(i)%tag, nm_p_sel)
     call select_par(gr_p, md_p(i)%tag, gr_p_sel)
     call select_par(sc_p, md_p(i)%tag, sc_p_sel)

     ! Set up the context

     cx(i) = context_t(ml, gr_p_sel, md_p(i), os_p_sel)

     ! Set up the frequency array

     call build_scan(cx(i), md_p(i), os_p_sel, sc_p_sel, omega)

     if (SIZE(omega) < 2) then

        if (check_log_level('INFO')) then
           write(OUTPUT_UNIT, 100) 'Scan is empty, skipping mode...'
        endif

        cycle md_p_loop

     endif

     ! Create the grid

     gr = grid_t(cx(i), omega, gr_p_sel)

     ! Set frequency bounds and perform checks

     if (nm_p_sel%restrict_roots) then
        omega_min = MINVAL(omega)
        omega_max = MAXVAL(omega)
     else
        omega_min = -HUGE(0._WP)
        omega_max = HUGE(0._WP)
     endif

     call check_scan(ml, gr, omega, md_p(i), os_p_sel)

     ! Find adiabatic modes

     if (md_p(i)%l == 0 .AND. os_p_sel%reduce_order) then
        allocate(bp_ad, SOURCE=rad_bvp_t(cx(i), gr, md_p(i), nm_p_sel, os_p_sel))
     else
        allocate(bp_ad, SOURCE=ad_bvp_t(cx(i), gr, md_p(i), nm_p_sel, os_p_sel))
     endif

     i_ad_a = n_md_ad + 1

     if (check_log_level('INFO')) then
        write(OUTPUT_UNIT, 100) 'Starting search (adiabatic)'
        write(OUTPUT_UNIT, *)
     endif

     call scan_search(bp_ad, omega, omega_min, omega_max, process_mode_ad, nm_p_sel)

     deallocate(bp_ad)

     ! Find non-adiabatic modes

     if (os_p_sel%nonadiabatic) then

        allocate(bp_nad, SOURCE=nad_bvp_t(cx(i), gr, md_p(i), nm_p_sel, os_p_sel))

        i_ad_b = n_md_ad

        if (check_log_level('INFO')) then
           write(OUTPUT_UNIT, 100) 'Starting search (non-adiabatic)'
           write(OUTPUT_UNIT, *)
        endif

        call prox_search(bp_nad, md_ad(i_ad_a:i_ad_b), omega_min, omega_max, process_mode_nad, md_p(i), nm_p_sel, os_p_sel)

        deallocate(bp_nad)

     endif

  end do md_p_loop

  ! Write summary files

  call write_summary(md_ad(:n_md_ad), ot_p_ad)
  call write_summary(md_nad(:n_md_nad), ot_p_nad)

  ! Clean up

  deallocate(md_ad)
  deallocate(md_nad)

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

    ! Store it

    n_md_ad = n_md_ad + 1

    if (n_md_ad > d_md_ad) then
       d_md_ad = 2*d_md_ad
       call reallocate(md_ad, [d_md_ad])
    endif

    md_ad(n_md_ad) = md

    ! Write it

    call write_mode(md_ad(n_md_ad), ot_p_ad)

    ! If necessary, prune it

    if (ot_p_ad%prune_modes) call md_ad(n_md_ad)%prune()

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

    ! Store it

    n_md_nad = n_md_nad + 1

    if (n_md_nad > d_md_nad) then
       d_md_nad = 2*d_md_nad
       call reallocate(md_nad, [d_md_nad])
    endif

    md_nad(n_md_nad) = md

    ! Write it

    call write_mode(md_nad(n_md_nad), ot_p_nad)

    ! If necessary, prune it

    if (ot_p_nad%prune_modes) call md_nad(n_md_nad)%prune()

    ! Finish

    return

  end subroutine process_mode_nad

end program gyre
