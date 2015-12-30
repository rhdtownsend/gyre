! Program  : gyre_ad
! Purpose  : adiabatic oscillation code
!
! Copyright 2013-2015 Rich Townsend
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

program gyre_ad

  ! Uses

  use core_kinds, SP_ => SP
  use gyre_constants
  use core_parallel

  use gyre_ad_bep
  use gyre_bep
  use gyre_ext
  use gyre_grid_par
  use gyre_input
  use gyre_model
  use gyre_model_par
  use gyre_mode
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_output
!  use gyre_rad_bvp
  use gyre_search
  use gyre_scan_par
  use gyre_mode
  use gyre_trad
  use gyre_util
  use gyre_version

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable      :: filename
  integer                        :: unit
  type(model_par_t)              :: ml_p
  type(mode_par_t), allocatable  :: md_p(:)
  type(osc_par_t), allocatable   :: os_p(:)
  type(num_par_t), allocatable   :: nm_p(:)
  type(grid_par_t), allocatable  :: gr_p(:)
  type(scan_par_t), allocatable  :: sc_p(:)
  type(out_par_t)                :: ot_p
  class(model_t), pointer        :: ml => null()
  integer                        :: i
  type(osc_par_t), allocatable   :: os_p_sel(:)
  type(num_par_t), allocatable   :: nm_p_sel(:)
  type(grid_par_t), allocatable  :: gr_p_sel(:)
  type(scan_par_t), allocatable  :: sc_p_sel(:)
  real(WP), allocatable          :: omega(:)
!  class(r_bep_t), allocatable    :: bp
  class(ad_bep_t), allocatable   :: bp
  integer                        :: n_md
  integer                        :: d_md
  type(mode_t), allocatable      :: md(:)

  ! Initialize

  call init_parallel()
  call init_system(filename)

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_ad ['//TRIM(version)//']', '=')
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

  ! Loop through mp

  d_md = 128
  n_md = 0

  allocate(md(d_md))

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
     call select_par(gr_p, md_p(i)%tag, gr_p_sel)
     call select_par(sc_p, md_p(i)%tag, sc_p_sel)

     $ASSERT(SIZE(os_p_sel) == 1,No matching osc parameters)
     $ASSERT(SIZE(nm_p_sel) == 1,No matching num parameters)
     $ASSERT(SIZE(gr_p_sel) >= 1,No matching grid parameters)
     $ASSERT(SIZE(sc_p_sel) >= 1,No matching scan parameters)

     ! Set up the frequency array

     call build_scan(ml, md_p(i), os_p_sel(1), sc_p_sel, omega)

     ! Set up bp

     !if (mp(i)%l == 0 .AND. os_p_sel(1)%reduce_order) then
     !   allocate(bp, SOURCE=rad_bvp_t(x_sh, ml, mp(i), os_p_sel(1), np_sel(1), omega_min, omega_max))
     !else
     allocate(bp, SOURCE=ad_bep_t(ml, omega, gr_p_sel, md_p(i), nm_p_sel(1), os_p_sel(1)))
     !endif

     ! Find roots

     call scan_search(bp, omega, process_root, nm_p_sel(1))

     ! Clean up

     deallocate(bp)

  end do md_p_loop

  ! Write the summary file
 
  call write_summary(md(:n_md), ot_p)

  ! Finish

  close(unit)

  call final_parallel()

contains

  subroutine process_root (omega, n_iter, discrim_ref)

    real(WP), intent(in)      :: omega
    integer, intent(in)       :: n_iter
    type(r_ext_t), intent(in) :: discrim_ref

    type(mode_t) :: md_new

    ! Create the mode

    md_new = mode_t(bp, omega)

    if (md_new%n_pg < md_p(i)%n_pg_min .OR. md_new%n_pg > md_p(i)%n_pg_max) return

    md_new%n_iter = n_iter
    md_new%chi = ABS(md_new%discrim)/ABS(discrim_ref)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) md_new%l, md_new%n_pg, md_new%n_p, md_new%n_g, &
            md_new%omega, real(md_new%chi), md_new%n_iter, md_new%n_k
120    format(4(2X,I8),3(2X,E24.16),2X,I6,2X,I7)
    endif

    ! Store it

    n_md = n_md + 1

    if (n_md > d_md) then
       d_md = 2*d_md
       call reallocate(md, [d_md])
    endif

    md(n_md) = md_new

    ! Write it

    call write_mode(md(n_md), n_md, ot_p)

    ! If necessary, prune it

    if (ot_p%prune_modes) call md(n_md)%prune()

    ! Finish

    return

  end subroutine process_root

end program gyre_ad
