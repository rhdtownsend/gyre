! Program  : gyre_nad
! Purpose  : nonadiabatic oscillation code
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

program gyre_nad

  ! Uses

  use core_kinds, SP_ => SP
  use gyre_constants
  use core_parallel

  use gyre_version
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
  use gyre_nad_bvp
  use gyre_search
  use gyre_mode
  use gyre_input
  use gyre_output
  use gyre_util

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
  integer                       :: n_op_sel
  integer                       :: n_np_sel
  real(WP), allocatable        :: omega(:)
  class(bvp_t), allocatable    :: ad_bp
  class(bvp_t), allocatable    :: nad_bp
  integer                      :: n_md_ad
  integer                      :: d_md_ad
  type(mode_t), allocatable    :: md_ad(:)
  integer                      :: n_md_nad
  integer                      :: d_md_nad
  type(mode_t), allocatable    :: md_nad(:)

  ! Initialize

  call init_parallel()

  call set_log_level($str($LOG_LEVEL))

  if(check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_nad ['//TRIM(version)//']', '=')
100  format(A)

     write(OUTPUT_UNIT, 110) 'Compiler         : ', COMPILER_VERSION()
     write(OUTPUT_UNIT, 110) 'Compiler options : ', COMPILER_OPTIONS()
110  format(2A)

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   : ', OMP_SIZE_MAX
120  format(A,I0)
     
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

  d_md_ad = 128
  n_md_ad = 0

  allocate(md_ad(d_md_ad))

  d_md_nad = 128
  n_md_nad = 0

  allocate(md_nad(d_md_nad))

  op_loop : do i = 1, SIZE(mp)

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

     ! Store the frequency range in shoot_gp_sel

     shoot_gp_sel%omega_a = MINVAL(omega)
     shoot_gp_sel%omega_b = MAXVAL(omega)

     ! Set up the bvp's

     if(ALLOCATED(ad_bp)) deallocate(ad_bp)

     if(mp(i)%l == 0 .AND. np_sel(n_np_sel)%reduce_order) then
        allocate(ad_bp, SOURCE=rad_bvp_t(ml, mp(i), op_sel(n_op_sel), np_sel(n_np_sel), shoot_gp_sel, recon_gp_sel, x_ml))
     else
        allocate(ad_bp, SOURCE=ad_bvp_t(ml, mp(i), op_sel(n_op_sel), np_sel(n_np_sel), shoot_gp_sel, recon_gp_sel, x_ml))
     endif
 
     allocate(nad_bp, SOURCE=nad_bvp_t(ml, mp(i), op_sel(n_op_sel), np_sel(n_np_sel), shoot_gp_sel, recon_gp_sel, x_ml))

     ! Find modes

     n_md_ad = 0

     call scan_search(ad_bp, np_sel(n_np_sel), omega, process_mode_ad)

     call prox_search(nad_bp, np_sel(n_np_sel), md_ad(:n_md_ad), process_mode_nad)

     ! Clean up

     deallocate(ad_bp)
     deallocate(nad_bp)

  end do op_loop

  ! Write the summary file
 
  call write_summary(up, md_nad(:n_md_nad))

  ! Finish

  close(unit)

  call final_parallel()

contains

  subroutine process_mode_ad (md_new)

    type(mode_t), intent(in) :: md_new

    ! Store the mode

    n_md_ad = n_md_ad + 1

    if (n_md_ad > d_md_ad) then
       d_md_ad = 2*d_md_ad
       call reallocate(md_ad, [d_md_ad])
    endif

    md_ad(n_md_ad) = md_new

    ! Prune it

    call md_ad(n_md_ad)%prune()

    ! Finish

    return

  end subroutine process_mode_ad

!****

  subroutine process_mode_nad (md_new)

    type(mode_t), intent(in) :: md_new

    ! Store the mode

    n_md_nad = n_md_nad + 1

    if (n_md_nad > d_md_nad) then
       d_md_nad = 2*d_md_nad
       call reallocate(md_nad, [d_md_nad])
    endif

    md_nad(n_md_nad) = md_new

    ! Write it

    call write_mode(up, md_nad(n_md_nad), n_md_nad)

    ! If necessary, prune it

    if (up%prune_modes) call md_nad(n_md_nad)%prune()

    ! Finish

    return

  end subroutine process_mode_nad

end program gyre_nad
