! Program  : gyre_nad_map
! Purpose  : nonadiabatic discriminant mapping code
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

program gyre_nad_map

  ! Uses

  use core_kinds
  use gyre_constants
  use core_parallel
  use core_order
  use core_hgroup

  use gyre_version
  use gyre_base_coeffs
  use gyre_therm_coeffs
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_bvp_nad
  use gyre_mode
  use gyre_input
  use gyre_output
  use gyre_util
  use gyre_ext_arith
  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  character(LEN=:), allocatable      :: filename
  integer                            :: unit
  real(WP), allocatable              :: x_bc(:)
  class(base_coeffs_t), allocatable  :: bc
  class(therm_coeffs_t), allocatable :: tc
  type(oscpar_t)                     :: op
  type(numpar_t)                     :: np
  real(WP), allocatable              :: omega_re(:)
  real(WP), allocatable              :: omega_im(:)
  type(gridpar_t), allocatable       :: shoot_gp(:)
  type(gridpar_t), allocatable       :: recon_gp(:)
  type(bvp_nad_t)                    :: nad_bp
  character(LEN=FILENAME_LEN)        :: map_file
  
  integer                   :: n_omega_re
  integer                   :: n_omega_im
  complex(WP), allocatable  :: discrim_map_f(:,:)
  integer, allocatable      :: discrim_map_e(:,:)
  integer, allocatable      :: k_part(:)
  integer                   :: k
  integer                   :: i(2)
  $if($MPI)
  integer                   :: p
  $endif
  complex(WP)               :: omega
  type(ext_complex_t)       :: discrim
  type(hgroup_t)            :: hg

  ! Initialize

  call init_parallel()

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
     
     write(OUTPUT_UNIT, 100) form_header('Initialization', '=')

  endif

  ! Process arguments

  if(MPI_RANK == 0) then

     call parse_args(filename)
     
     open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

     call read_coeffs(unit, x_bc, bc, tc)

     $ASSERT(ALLOCATED(tc),No therm_coeffs data)

     call read_oscpar(unit, op)
     call read_numpar(unit, np)
     call read_shoot_gridpar(unit, shoot_gp)
     call read_recon_gridpar(unit, recon_gp)
     call read_mappar(unit, bc, op, shoot_gp, x_bc, omega_re, omega_im)
     call read_outpar(unit, map_file)

     call nad_bp%init(bc, op, np, shoot_gp, recon_gp, x_bc, tc)

  endif

  $if($MPI)
  call bcast_constants(0)
  call bcast_alloc(omega_re, 0)
  call bcast_alloc(omega_im, 0)
  call bcast(nad_bp, 0)
  $endif

  ! Map the discriminant

  n_omega_re = SIZE(omega_re)
  n_omega_im = SIZE(omega_im)

  allocate(discrim_map_f(n_omega_re,n_omega_im))
  allocate(discrim_map_e(n_omega_re,n_omega_im))

  allocate(k_part(MPI_SIZE+1))

  call partition_tasks(n_omega_re*n_omega_im, 1, k_part)

  do k = k_part(MPI_RANK+1), k_part(MPI_RANK+2)-1

     i = index_nd(k, [n_omega_re,n_omega_im])

     omega = CMPLX(omega_re(i(1)), omega_im(i(2)), KIND=WP)

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

  call hg%init(map_file, CREATE_FILE)

  call write_dset(hg, 'omega_re', omega_re)
  call write_dset(hg, 'omega_im', omega_im)

  call write_dset(hg, 'discrim_map_f', discrim_map_f)
  call write_dset(hg, 'discrim_map_e', discrim_map_e)

  call hg%final()

  ! Finish

  call final_parallel()

contains

  subroutine read_mappar (unit, bc, op, gp, x_in, omega_re, omega_im)

    integer, intent(in)                :: unit
    class(base_coeffs_t), intent(in)   :: bc
    type(oscpar_t), intent(in)         :: op
    type(gridpar_t), intent(inout)     :: gp(:)
    real(WP), allocatable, intent(in)  :: x_in(:)
    real(WP), allocatable, intent(out) :: omega_re(:)
    real(WP), allocatable, intent(out) :: omega_im(:)

    character(LEN=256) :: grid_type_re
    character(LEN=256) :: grid_type_im
    real(WP)           :: freq_re_min
    real(WP)           :: freq_re_max
    integer            :: n_freq_re
    real(WP)           :: freq_im_min
    real(WP)           :: freq_im_max
    integer            :: n_freq_im
    character(LEN=256) :: freq_units
    real(WP)           :: x_i
    real(WP)           :: x_o
    real(WP)           :: omega_min
    real(WP)           :: omega_max
    integer            :: i

    namelist /map/ grid_type_re, grid_type_im, &
          freq_re_min, freq_re_max, n_freq_re, &
          freq_im_min, freq_im_max, n_freq_im, &
          freq_units

    ! Determine the grid range

    call grid_range(gp, bc, op, x_in, x_i, x_o)

    ! Read map parameters

    rewind(unit)

    grid_type_re = 'LINEAR'
    grid_type_im = 'LINEAR'

    freq_re_min = 1._WP
    freq_re_max = 11._WP
    n_freq_re = 10
          
    freq_im_min = -5._WP
    freq_im_max = 5._WP
    n_freq_im = 10
          
    freq_units = 'NONE'

    read(unit, NML=map)
          
    ! Set up the frequency grids

    omega_min = freq_re_min/freq_scale(bc, op, x_o, freq_units)
    omega_max = freq_re_max/freq_scale(bc, op, x_o, freq_units)
       
    select case(grid_type_re)
    case('LINEAR')
       omega_re = [(((n_freq_re-i)*omega_min + (i-1)*omega_max)/(n_freq_re-1), i=1,n_freq_re)]
    case('INVERSE')
       omega_re = [((n_freq_re-1)/((n_freq_re-i)/omega_min + (i-1)/omega_max), i=1,n_freq_re)]
    case default
       $ABORT(Invalid grid_type_re)
    end select

    omega_min = freq_im_min/freq_scale(bc, op, x_o, freq_units)
    omega_max = freq_im_max/freq_scale(bc, op, x_o, freq_units)
       
    select case(grid_type_im)
    case('LINEAR')
       omega_im = [(((n_freq_im-i)*omega_min + (i-1)*omega_max)/(n_freq_im-1), i=1,n_freq_im)]
    case('INVERSE')
       omega_im = [((n_freq_im-1)/((n_freq_im-i)/omega_min + (i-1)/omega_max), i=1,n_freq_im)]
    case default
       $ABORT(Invalid grid_type_im)
    end select

    ! Store the frequency range in gp

    gp%omega_a = MINVAL(omega_re)
    gp%omega_b = MAXVAL(omega_re)

    ! Finish

    return

  end subroutine read_mappar

!****

  subroutine read_outpar (unit, file)

    integer, intent(in)           :: unit
    character(LEN=*), intent(out) :: file

    namelist /output/ file

    ! Read output parameters

    rewind(unit)

    read(unit, NML=output)

    ! Finish
    
    return

  end subroutine read_outpar
          
end program gyre_nad_map
