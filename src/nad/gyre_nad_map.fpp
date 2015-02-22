! Program  : gyre_nad_map
! Purpose  : nonadiabatic discriminant mapping code
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

program gyre_nad_map

  ! Uses

  use core_kinds, SP_ => SP
  use gyre_constants
  use core_hgroup
  use core_order
  use core_parallel

  use gyre_bvp
  use gyre_ext
  use gyre_grid
  use gyre_grid_par
  use gyre_input
  use gyre_mode
  use gyre_model
  $if($MPI)
  use gyre_model_mpi
  $endif
  use gyre_mode_par
  use gyre_nad_bvp
  use gyre_num_par
  use gyre_osc_par
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
  character(:), allocatable     :: gyre_dir
  integer                       :: unit
  real(WP), allocatable         :: x_ml(:)
  class(model_t), pointer       :: ml => null()
  type(mode_par_t), allocatable :: mp(:)
  type(osc_par_t), allocatable  :: op(:)
  type(num_par_t), allocatable  :: np(:)
  type(grid_par_t), allocatable :: shoot_gp(:)
  type(grid_par_t), allocatable :: recon_gp(:)
  type(scan_par_t), allocatable :: sp(:)
  type(scan_par_t), allocatable :: sp_re(:)
  type(scan_par_t), allocatable :: sp_im(:)
  real(WP)                      :: x_i
  real(WP)                      :: x_o
  real(WP), allocatable         :: omega_re(:)
  real(WP), allocatable         :: omega_im(:)
  real(WP), allocatable         :: x_sh(:)
  class(c_bvp_t), allocatable   :: nad_bp
  integer                       :: n_omega_re
  integer                       :: n_omega_im
  complex(WP), allocatable      :: discrim_map_f(:,:)
  integer, allocatable          :: discrim_map_e(:,:)
  integer, allocatable          :: k_part(:)
  integer                       :: k
  integer                       :: i(2)
  $if($MPI)
  integer                       :: p
  $endif
  complex(WP)                   :: omega
  type(c_ext_t)                 :: discrim
  character(FILENAME_LEN)       :: map_filename
  integer                       :: i_percent
  integer                       :: n_percent
  type(hgroup_t)                :: hg

  ! Initialize

  call init_parallel()
  call init_system(filename, gyre_dir)

  call set_log_level($str($LOG_LEVEL))

  if(check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_nad_map ['//TRIM(version)//']', '=')
100  format(A)

     write(OUTPUT_UNIT, 110) 'Compiler         : ', COMPILER_VERSION()
     write(OUTPUT_UNIT, 110) 'Compiler options : ', COMPILER_OPTIONS()
110  format(2A)

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   : ', OMP_SIZE_MAX
     write(OUTPUT_UNIT, 120) 'MPI Processors   : ', MPI_SIZE
120  format(A,I0)
     
     write(OUTPUT_UNIT, 110) 'Input filename   :', filename
     write(OUTPUT_UNIT, 110) 'GYRE_DIR         :', gyre_dir

     write(OUTPUT_UNIT, 100) form_header('Initialization', '=')

  endif

  call init_trad(gyre_dir)

  ! Process arguments

  if(MPI_RANK == 0) then

     open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

     call read_constants(unit)
     call read_model(unit, x_ml, ml)
     call read_mode_par(unit, mp)
     call read_osc_par(unit, op)
     call read_num_par(unit, np)
     call read_shoot_grid_par(unit, shoot_gp)
     call read_recon_grid_par(unit, recon_gp)
     call read_scan_par(unit, sp)
     call read_out_par(unit, map_filename)

  endif

  $if($MPI)
  call bcast_constants(0)
  call bcast_alloc(x_ml, 0)
  call bcast_alloc(ml, 0)
  call bcast_alloc(mp, 0)
  call bcast_alloc(op, 0)
  call bcast_alloc(np, 0)
  call bcast_alloc(shoot_gp, 0)
  call bcast_alloc(recon_gp, 0)
  call bcast_alloc(sp, 0)
  $endif

  ! Select parameters according to tags

  call select_par(sp, 'REAL', sp_re)
  call select_par(sp, 'IMAG', sp_im)

  $ASSERT(SIZE(op) == 1,No matching osc parameters)
  $ASSERT(SIZE(np) == 1,No matching num parameters)
  $ASSERT(SIZE(shoot_gp) >= 1,No matching shoot_grid parameters)
  $ASSERT(SIZE(recon_gp) >= 1,No matching recon_grid parameters)
  $ASSERT(SIZE(sp_im) == 1,No matching scan parameters)
  $ASSERT(SIZE(sp_re) == 1,No matching scan parameters)

  ! Set up the frequency arrays

  if (allocated(x_ml)) then
     x_i = x_ml(1)
     x_o = x_ml(SIZE(x_ml))
  else
     x_i = 0._WP
     x_o = 1._WP
  endif

  call build_scan(sp_re, ml, mp(1), op(1), x_i, x_o, omega_re)
  call build_scan(sp_im, ml, mp(1), op(1), x_i, x_o, omega_im)

  ! Set up the shooting grid

  call build_grid(shoot_gp, ml, mp(1), op(1), omega_re, x_ml, x_sh, verbose=.TRUE.)

  ! Set up the bvp

  allocate(nad_bp, SOURCE=nad_bvp_t(x_sh, ml, mp(1), op(1), np(1)))

  ! Map the discriminant

  n_omega_re = SIZE(omega_re)
  n_omega_im = SIZE(omega_im)

  allocate(discrim_map_f(n_omega_re,n_omega_im))
  allocate(discrim_map_e(n_omega_re,n_omega_im))

  allocate(k_part(MPI_SIZE+1))

  call partition_tasks(n_omega_re*n_omega_im, 1, k_part)

  n_percent = 0

  do k = k_part(MPI_RANK+1), k_part(MPI_RANK+2)-1

     i = index_nd(k, [n_omega_re,n_omega_im])

     omega = CMPLX(omega_re(i(1)), omega_im(i(2)), KIND=WP)

     if (MPI_RANK == 0) then
        i_percent = FLOOR(100._WP*REAL(k-k_part(MPI_RANK+1))/REAL(k_part(MPI_RANK+2)-k_part(MPI_RANK+1)-1))
        if (i_percent > n_percent) then
           print *,'Percent complete: ', i_percent
           n_percent = i_percent
        end if
     endif
     
     discrim = nad_bp%discrim(omega)

     discrim_map_f(i(1),i(2)) = FRACTION(discrim)
     discrim_map_e(i(1),i(2)) = EXPONENT(discrim)

  end do

  $if($MPI)

  do p = 1,MPI_SIZE
     call bcast_seq(discrim_map_f, k_part(p), k_part(p+1)-1, p-1)
     call bcast_seq(discrim_map_e, k_part(p), k_part(p+1)-1, p-1)
  end do

  $endif

  ! Write out the map

  if (MPI_RANK == 0) then

     hg = hgroup_t(map_filename, CREATE_FILE)

     call write_dset(hg, 'omega_re', omega_re)
     call write_dset(hg, 'omega_im', omega_im)

     call write_dset(hg, 'discrim_map_f', discrim_map_f)
     call write_dset(hg, 'discrim_map_e', discrim_map_e)

     call hg%final()

  end if

  ! Finish

  call final_parallel()

contains

  subroutine read_out_par (unit, map_file)

    integer, intent(in)       :: unit
    character(*), intent(out) :: map_file

    namelist /output/ map_file

    ! Read output parameters

    rewind(unit)

    read(unit, NML=output)

    ! Finish
    
    return

  end subroutine read_out_par
          
end program gyre_nad_map
