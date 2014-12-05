! Program  : gyre_ad
! Purpose  : adiabatic oscillation code
!
! Copyright 2013 Rich Townsend
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

  use gyre_version
  use gyre_ext
  use gyre_model
  use gyre_modepar
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_scanpar
  use gyre_outpar
  use gyre_bvp
  use gyre_ad_bvp
  use gyre_rad_bvp
  use gyre_search
  use gyre_mode
  use gyre_input
  use gyre_output
  use gyre_util
  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(:), allocatable    :: filename
  integer                      :: unit
  real(WP), allocatable        :: x_ml(:)
  class(model_t), pointer      :: ml => null()
  type(modepar_t), allocatable :: mp(:)
  type(oscpar_t), allocatable  :: op(:)
  type(numpar_t), allocatable  :: np(:)
  type(gridpar_t), allocatable :: shoot_gp(:)
  type(gridpar_t), allocatable :: recon_gp(:)
  type(scanpar_t), allocatable :: sp(:)
  type(outpar_t)               :: up
  integer                      :: i
  type(oscpar_t), allocatable  :: op_sel(:)
  type(numpar_t), allocatable  :: np_sel(:)
  type(gridpar_t), allocatable :: shoot_gp_sel(:)
  type(gridpar_t), allocatable :: recon_gp_sel(:)
  type(scanpar_t), allocatable :: sp_sel(:)
  integer                      :: n_op_sel
  integer                      :: n_np_sel
  real(WP), allocatable        :: omega(:)
  real(WP), allocatable        :: x_sh(:)
  class(r_bvp_t), allocatable  :: bp
  integer                      :: n_md
  integer                      :: d_md
  type(mode_t), allocatable    :: md(:)

  ! Initialize

  call init_parallel()

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_ad ['//TRIM(version)//']', '=')
100  format(A)

     write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
     write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110  format(A,1X,A)

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
120  format(A,1X,I0)

     write(OUTPUT_UNIT, 100) form_header('Initialization', '=')

  endif

  ! Process arguments

  call parse_args(filename)
     
  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  call read_model(unit, x_ml, ml)
  call read_constants(unit)
  call read_modepar(unit, mp)
  call read_oscpar(unit, op)
  call read_numpar(unit, np)
  call read_shoot_gridpar(unit, shoot_gp)
  call read_recon_gridpar(unit, recon_gp)
  call read_scanpar(unit, sp)
  call read_outpar(unit, up)

  ! Loop through modepars

  d_md = 128
  n_md = 0

  allocate(md(d_md))

  mp_loop : do i = 1, SIZE(mp)

     if (check_log_level('INFO')) then

        write(OUTPUT_UNIT, 100) form_header('Mode Search', '=')

        write(OUTPUT_UNIT, 100) 'Mode parameters'

        write(OUTPUT_UNIT, 130) 'l :', mp(i)%l
        write(OUTPUT_UNIT, 130) 'm :', mp(i)%m
130     format(3X,A,1X,I0)

        write(OUTPUT_UNIT, *)

     endif

     ! Select parameters according to tags

     call select_par(op, mp(i)%tag, op_sel, last=.TRUE.)
     call select_par(np, mp(i)%tag, np_sel, last=.TRUE.)
     call select_par(shoot_gp, mp(i)%tag, shoot_gp_sel)
     call select_par(recon_gp, mp(i)%tag, recon_gp_sel)
     call select_par(sp, mp(i)%tag, sp_sel)

     n_op_sel = SIZE(op_sel)
     n_np_sel = SIZE(np_sel)

     $ASSERT(n_op_sel >= 1,No matching osc parameters)
     $ASSERT(n_np_sel >= 1,No matching num parameters)
     $ASSERT(SIZE(shoot_gp_sel) >= 1,No matching shoot_grid parameters)
     $ASSERT(SIZE(recon_gp_sel) >= 1,No matching recon_grid parameters)
     $ASSERT(SIZE(sp_sel) >= 1,No matching scan parameters)

     if (n_op_sel > 1 .AND. check_log_level('WARN')) then
        write(OUTPUT_UNIT, 140) 'Warning: multiple matching osc namelists, using final match'
140     format('!!',1X,A)
     endif
        
     if (n_np_sel > 1 .AND. check_log_level('WARN')) then
        write(OUTPUT_UNIT, 140) 'Warning: multiple matching num namelists, using final match'
     endif
        
     ! Set up the frequency array

     call build_scan(sp_sel, ml, mp(i), op_sel(n_op_sel), shoot_gp_sel, x_ml, omega)

     ! Set up the shooting grid

     call build_grid(shoot_gp_sel, ml, mp(i), MINVAL(omega), MAXVAL(omega), x_ml, x_sh, verbose=.TRUE.)

     ! Set up bp

     if (mp(i)%l == 0 .AND. op_sel(n_op_sel)%reduce_order) then
        allocate(bp, SOURCE=rad_bvp_t(x_sh, ml, mp(i), op_sel(n_op_sel), np_sel(n_np_sel)))
     else
        allocate(bp, SOURCE=ad_bvp_t(x_sh, ml, mp(i), op_sel(n_op_sel), np_sel(n_np_sel)))
     endif

     ! Find roots

     call scan_search(bp, np_sel(n_np_sel), omega, process_root)

     ! Clean up

     deallocate(bp)

  end do mp_loop

  ! Write the summary file
 
  call write_summary(up, md(:n_md))

  ! Finish

  close(unit)

  call final_parallel()

contains

  subroutine process_root (omega, n_iter, discrim_ref)

    real(WP), intent(in)      :: omega
    integer, intent(in)       :: n_iter
    type(r_ext_t), intent(in) :: discrim_ref

    real(WP), allocatable :: x_rc(:)
    integer               :: n
    real(WP)              :: x_ref
    real(WP), allocatable :: y(:,:)
    real(WP)              :: y_ref(6)
    type(r_ext_t)         :: discrim
    type(mode_t)          :: md_new

    ! Build the reconstruction grid

    call build_grid(recon_gp_sel, ml, mp(i), omega, omega, x_sh, x_rc, verbose=.FALSE.)

    ! Reconstruct the solution

    x_ref = MIN(MAX(op_sel(n_op_sel)%x_ref, x_sh(1)), x_sh(SIZE(x_sh)))

    n = SIZE(x_rc)

    allocate(y(6,n))

    call bp%recon(omega, x_rc, x_ref, y, y_ref, discrim)

    ! Create the mode

    md_new = mode_t(ml, mp(i), op_sel(n_op_sel), CMPLX(omega, KIND=WP), c_ext_t(discrim), &
                    x_rc, CMPLX(y, KIND=WP), x_ref, CMPLX(y_ref, KIND=WP))

    if (md_new%n_pg < mp(i)%n_pg_min .OR. md_new%n_pg > mp(i)%n_pg_max) return

    md_new%n_iter = n_iter
    md_new%chi = ABS(discrim)/ABS(discrim_ref)

    if (check_log_level('INFO')) then
       write(OUTPUT_UNIT, 120) md_new%mp%l, md_new%n_pg, md_new%n_p, md_new%n_g, &
            md_new%omega, real(md_new%chi), md_new%n_iter, md_new%n
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

    call write_mode(up, md(n_md), n_md)

    ! If necessary, prune it

    if (up%prune_modes) call md(n_md)%prune()

    ! Finish

    return

  end subroutine process_root

end program gyre_ad
