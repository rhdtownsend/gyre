! Program  : gyre_nad
! Purpose  : nonadiabatic oscillation code
!
! Copyright 2013-2016 Rich Townsend
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

program gyre_nad

  ! Uses

  use core_kinds, SP_ => SP
  use gyre_constants
  use core_parallel

  use gyre_ad_bep
  use gyre_bep
  use gyre_ext
  use gyre_input
  use gyre_grid_par
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_model_par
  use gyre_nad_bep
  use gyre_num_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_output
  use gyre_rad_bep
  use gyre_scan_par
  use gyre_search
  use gyre_trad
  use gyre_util
  use gyre_version

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable     :: filename
  integer                       :: unit
  type(model_par_t)             :: ml_p
  type(mode_par_t), allocatable :: md_p(:)
  type(osc_par_t), allocatable  :: os_p(:)
  type(num_par_t), allocatable  :: nm_p(:)
  type(grid_par_t), allocatable :: gr_p(:)
  type(scan_par_t), allocatable :: sc_p(:)
  type(out_par_t)               :: ot_p
  class(model_t), pointer       :: ml => null()
  integer                       :: i
  type(osc_par_t), allocatable  :: os_p_sel(:)
  type(num_par_t), allocatable  :: nm_p_sel(:)
  type(grid_par_t), allocatable :: gr_p_sel(:)
  type(scan_par_t), allocatable :: sc_p_sel(:)
  real(WP), allocatable         :: omega(:)
  class(r_bep_t), allocatable   :: bp_ad
  class(c_bep_t), allocatable   :: bp_nad
  integer                       :: n_md_ad
  integer                       :: d_md_ad
  type(mode_t), allocatable     :: md_ad(:)
  integer                       :: n_md_nad
  integer                       :: d_md_nad
  type(mode_t), allocatable     :: md_nad(:)

  ! Initialize

  call init_parallel()
  call init_system(filename)

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_nad ['//TRIM(version)//']', '=')
100  format(A)

     write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
     write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110  format(A,1X,A)

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
120  format(A,1X,I0)
     
     write(OUTPUT_UNIT, 110) 'Input filename   :', filename
     write(OUTPUT_UNIT, 110) 'GYRE_DIR         :', gyre_dir

     write(OUTPUT_UNIT, 100) form_header('Initialization', '=')

  endif

  ! Process arguments

  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  call read_model_par(unit, ml_p)
  call read_constants(unit)
  call read_mode_par(unit, md_p)
  call read_osc_par(unit, os_p)
  call read_num_par(unit, nm_p)
  call read_grid_par(unit, gr_p)
  call read_scan_par(unit, sc_p)
  call read_out_par(unit, ot_p)

  ! Read the model

  call read_model(ml_p, ml)

  ! Loop through md_p

  d_md_ad = 128
  n_md_ad = 0

  allocate(md_ad(d_md_ad))

  d_md_nad = 128
  n_md_nad = 0

  allocate(md_nad(d_md_nad))

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

     call select_par(os_p, md_p(i)%tag, os_p_sel, last=.TRUE.)
     call select_par(nm_p, md_p(i)%tag, nm_p_sel, last=.TRUE.)
     call select_par(gr_p, md_p(i)%tag, gr_p_sel, last=.TRUE.)
     call select_par(sc_p, md_p(i)%tag, sc_p_sel)

     $ASSERT(SIZE(os_p_sel) == 1,No matching osc parameters)
     $ASSERT(SIZE(nm_p_sel) == 1,No matching num parameters)
     $ASSERT(SIZE(gr_p_sel) == 1,No matching grid parameters)
     $ASSERT(SIZE(sc_p_sel) >= 1,No matching scan parameters)

     ! Set up the frequency array

     call build_scan(ml, md_p(i), os_p_sel(1), sc_p_sel, omega)

     ! Set up bp_ad and bp_nad

     if (md_p(i)%l == 0 .AND. os_p_sel(1)%reduce_order) then
        allocate(bp_ad, SOURCE=rad_bep_t(ml, omega, gr_p_sel(1), md_p(i), nm_p_sel(1), os_p_sel(1)))
     else
        allocate(bp_ad, SOURCE=ad_bep_t(ml, omega, gr_p_sel(1), md_p(i), nm_p_sel(1), os_p_sel(1)))
     endif

     allocate(bp_nad, SOURCE=nad_bep_t(ml, omega, gr_p_sel(1), md_p(i), nm_p_sel(1), os_p_sel(1)))

     ! Find modes

     n_md_ad = 0

     call scan_search(bp_ad, omega, process_root_ad, nm_p_sel(1))

     call prox_search(bp_nad, md_ad(:n_md_ad), process_root_nad, md_p(i), nm_p_sel(1), os_p_sel(1))

     ! Clean up

     deallocate(bp_ad)
     deallocate(bp_nad)

  end do md_p_loop

  ! Write the summary file
 
  call write_summary(md_nad(:n_md_nad), ot_p)

  ! Finish

  close(unit)

  call final_parallel()

contains

  subroutine process_root_ad (omega, n_iter, discrim_ref)

    real(WP), intent(in)      :: omega
    integer, intent(in)       :: n_iter
    type(r_ext_t), intent(in) :: discrim_ref

    type(sol_t)   :: sl
    type(mode_t)  :: md_new
    type(r_ext_t) :: chi

    ! Create the sol_t

    select type (bp_ad)
    type is (ad_bep_t)
       sl = sol_t(bp_ad, omega)
    type is (rad_bep_t)
       sl = sol_t(bp_ad, omega)
    class default
       $ABORT(Invalid bp class)
    end select

    ! Construct the new mode

    md_new = mode_t(ml, sl, md_p(i), os_p_sel(1))

    if (md_new%n_pg < md_p(i)%n_pg_min .OR. md_new%n_pg > md_p(i)%n_pg_max) return

    chi = ABS(sl%discrim)/ABS(discrim_ref)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) md_new%l, md_new%n_pg, md_new%n_p, md_new%n_g, &
            md_new%omega, real(chi), n_iter
120    format(4(2X,I8),3(2X,E24.16),2X,I6)
    endif

    ! Store it

    n_md_ad = n_md_ad + 1

    if (n_md_ad > d_md_ad) then
       d_md_ad = 2*d_md_ad
       call reallocate(md_ad, [d_md_ad])
    endif

    md_ad(n_md_ad) = md_new

    ! If necessary, prune it

    if (ot_p%prune_modes) call md_ad(n_md_ad)%prune()

    ! Finish

    return

  end subroutine process_root_ad

  !****

  subroutine process_root_nad (omega, n_iter, discrim_ref)

    complex(WP), intent(in)   :: omega
    integer, intent(in)       :: n_iter
    type(r_ext_t), intent(in) :: discrim_ref

    type(sol_t)   :: sl
    type(mode_t)  :: md_new
    type(r_ext_t) :: chi

    ! Create the sol_t

    select type (bp_nad)
    type is (nad_bep_t)
       sl = sol_t(bp_nad, omega)
    class default
       $ABORT(Invalid bp_nad class)
    end select

    ! Construct the new mode

    md_new = mode_t(ml, sl, md_p(i), os_p_sel(1))

    chi = ABS(sl%discrim)/ABS(discrim_ref)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) md_new%l, md_new%n_pg, md_new%n_p, md_new%n_g, &
            md_new%omega, real(chi), n_iter, md_new%n_k
120    format(4(2X,I8),3(2X,E24.16),2X,I6)
    endif

    ! Store it

    n_md_nad = n_md_nad + 1

    if (n_md_nad > d_md_nad) then
       d_md_nad = 2*d_md_nad
       call reallocate(md_nad, [d_md_nad])
    endif

    md_nad(n_md_nad) = md_new

    ! Write it

    call write_mode(md_nad(n_md_nad), n_md_nad, ot_p)

    ! If necessary, prune it

    if (ot_p%prune_modes) call md_nad(n_md_nad)%prune()

    ! Finish

    return

    ! Finish

    return

  end subroutine process_root_nad

end program gyre_nad
