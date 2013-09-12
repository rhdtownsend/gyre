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
  use core_constants
  use core_parallel

  use gyre_version
  use gyre_coeffs
  $if($MPI)
  use gyre_coeffs_mpi
  $endif
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_scanpar
  use gyre_bvp
  use gyre_bvp_ad
  use gyre_bvp_rad
  use gyre_search
  use gyre_mode
  use gyre_input
  use gyre_output
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(LEN=:), allocatable :: filename
  integer                       :: unit
  real(WP), allocatable         :: x_cf(:)
  class(coeffs_t), allocatable  :: cf
  type(oscpar_t), allocatable   :: op(:)
  type(numpar_t), allocatable   :: np(:)
  type(gridpar_t), allocatable  :: shoot_gp(:)
  type(gridpar_t), allocatable  :: recon_gp(:)
  type(scanpar_t), allocatable  :: sp(:)
  integer                       :: i
  type(numpar_t), allocatable   :: np_sel(:)
  type(gridpar_t), allocatable  :: shoot_gp_sel(:)
  type(gridpar_t), allocatable  :: recon_gp_sel(:)
  type(scanpar_t), allocatable  :: sp_sel(:)
  real(WP), allocatable         :: omega(:)
  class(bvp_t), allocatable     :: bp
  type(mode_t), allocatable     :: md(:)
  type(mode_t), allocatable     :: md_all(:)
  type(mode_t), allocatable     :: md_tmp(:)

  ! Initialize

  call init_parallel()

  call set_log_level($str($LOG_LEVEL))

  if(check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_ad ['//TRIM(version)//']', '=')
100  format(A)

     write(OUTPUT_UNIT, 110) 'Compiler         : ', COMPILER_VERSION()
     write(OUTPUT_UNIT, 110) 'Compiler options : ', COMPILER_OPTIONS()
110  format(2A)

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   : ', OMP_SIZE_MAX
     write(OUTPUT_UNIT, 120) 'MPI Processors   : ', MPI_SIZE
120  format(A,I0)

     write(OUTPUT_UNIT, 100) form_header('Initialization', '=')

  endif

  ! Process arguments

  if(MPI_RANK == 0) then

     call parse_args(filename)
     
     open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

     call read_coeffs(unit, x_cf, cf)
     call read_oscpar(unit, op)
     call read_numpar(unit, np)
     call read_shoot_gridpar(unit, shoot_gp)
     call read_recon_gridpar(unit, recon_gp)
     call read_scanpar(unit, sp)

  end if

  $if($MPI)
  call bcast_alloc(x_cf, 0)
  call bcast_alloc(cf, 0)
  call bcast_alloc(op, 0)
  call bcast_alloc(np, 0)
  call bcast_alloc(shoot_gp, 0)
  call bcast_alloc(recon_gp, 0)
  call bcast_alloc(sp, 0)
  $endif

  ! Loop through oscpars

  allocate(md_all(0))

  op_loop : do i = 1, SIZE(op)

     ! Select parameters according to tags

     call select_par(np, op(i)%tag, np_sel, last=.TRUE.)
     call select_par(shoot_gp, op(i)%tag, shoot_gp_sel)
     call select_par(recon_gp, op(i)%tag, recon_gp_sel)
     call select_par(sp, op(i)%tag, sp_sel)

     $ASSERT(SIZE(np_sel) == 1,No matching num parameters)
     $ASSERT(SIZE(shoot_gp_sel) >= 1,No matching shoot_grid parameters)
     $ASSERT(SIZE(recon_gp_sel) >= 1,No matching recon_grid parameters)
     $ASSERT(SIZE(sp_sel) >= 1,No matching scan parameters)

     ! Set up the frequency array

     call build_scan(sp_sel, cf, op(i), shoot_gp_sel, x_cf, omega)

     ! Store the frequency range in shoot_gp

     shoot_gp_sel%omega_a = MINVAL(omega)
     shoot_gp_sel%omega_b = MAXVAL(omega)

     ! Set up bp

     if(ALLOCATED(bp)) deallocate(bp)

     if(op(i)%l == 0 .AND. np_sel(1)%reduce_order) then
        allocate(bvp_rad_t::bp)
     else
        allocate(bvp_ad_t::bp)
     endif

     call bp%init(cf, op(i), np_sel(1), shoot_gp_sel, recon_gp_sel, x_cf)

     ! Find modes

     call scan_search(bp, omega, md)

     ! (The following could be simpler, but this is a workaround for a
     ! gfortran issue, likely PR 57839)

     md_tmp = [md_all, md]
     call MOVE_ALLOC(md_tmp, md_all)

  end do op_loop

  ! Write output
 
  if(MPI_RANK == 0) then
     call write_data(unit, md_all)
  endif

  ! Finish

  close(unit)

  call final_parallel()

end program gyre_ad
