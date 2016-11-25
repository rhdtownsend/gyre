! Program  : gyre_map
! Purpose  : discriminant mapping code
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

program gyre_map

  ! Uses

  use core_kinds, only : WP
  use core_hgroup
  use core_order
  use core_parallel
  use core_system

  use gyre_bep
  use gyre_constants
  use gyre_ext
  use gyre_grid
  use gyre_grid_factory
  use gyre_grid_par
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_model_factory
  use gyre_model_par
  use gyre_nad_bep
  use gyre_num_par
  use gyre_osc_par
  use gyre_scan_par, only : scan_par_t
  use gyre_search
  use gyre_status
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
  type(scan_par_t), allocatable :: sc_p_re(:)
  type(scan_par_t), allocatable :: sc_p_im(:)
  character(FILENAME_LEN)       :: map_filename
  class(model_t), pointer       :: ml => null()
  type(osc_par_t)               :: os_p_sel
  type(num_par_t)               :: nm_p_sel
  type(grid_par_t)              :: gr_p_sel
  type(scan_par_t), allocatable :: sc_p_re_sel(:)
  type(scan_par_t), allocatable :: sc_p_im_sel(:)
  type(grid_t)                  :: gr
  real(WP), allocatable         :: omega_re(:)
  real(WP), allocatable         :: omega_im(:)
  type(nad_bep_t)               :: bp
  integer                       :: n_omega_re
  integer                       :: n_omega_im
  complex(WP), allocatable      :: discrim_map_f(:,:)
  integer, allocatable          :: discrim_map_e(:,:)
  integer, allocatable          :: k_part(:)
  integer                       :: k
  integer                       :: i(2)
  $if ($MPI)
  integer                       :: p
  $endif
  complex(WP)                   :: omega
  type(c_ext_t)                 :: discrim
  integer                       :: status
  integer                       :: i_percent
  integer                       :: n_percent
  type(hgroup_t)                :: hg

  ! Read command-line arguments

  $ASSERT(n_arg() == 1,Syntax: gyre_map <filename>)

  call get_arg(1, filename)

  ! Initialize

  call init_parallel()

  call set_log_level($str($LOG_LEVEL))

  if (check_log_level('INFO')) then

     write(OUTPUT_UNIT, 100) form_header('gyre_map ['//VERSION//']', '-')
100  format(A)

     if (check_log_level('DEBUG')) then
        write(OUTPUT_UNIT, 110) 'Compiler         :', COMPILER_VERSION()
        write(OUTPUT_UNIT, 110) 'Compiler options :', COMPILER_OPTIONS()
110     format(A,1X,A)
     endif

     write(OUTPUT_UNIT, 120) 'OpenMP Threads   :', OMP_SIZE_MAX
     write(OUTPUT_UNIT, 120) 'MPI Processors   :', MPI_SIZE
120  format(A,1X,I0)
     
     write(OUTPUT_UNIT, 110) 'Input filename   :', filename
     write(OUTPUT_UNIT, 110) 'GYRE_DIR         :', gyre_dir

     write(OUTPUT_UNIT, 100) form_header('Initialization', '-')

  endif

  ! Process arguments

  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  call read_constants(unit)

  call read_model_par(unit, ml_p)
  call read_mode_par(unit, md_p)
  call read_osc_par(unit, os_p)
  call read_num_par(unit, nm_p)
  call read_grid_par(unit, gr_p)
  call read_scan_par(unit, 're', sc_p_re)
  call read_scan_par(unit, 'im', sc_p_im)
  call read_out_par(unit, map_filename)

  $ASSERT(SIZE(md_p) == 1,Must be exactly one mode parameter)

  ! Construct the model

  ml => model_t(ml_p)

  ! Select parameters according to tags

  call select_par(os_p, md_p(1)%tag, os_p_sel)
  call select_par(nm_p, md_p(1)%tag, nm_p_sel)
  call select_par(gr_p, md_p(1)%tag, gr_p_sel)
  call select_par(sc_p_re, md_p(1)%tag, sc_p_re_sel)
  call select_par(sc_p_im, md_p(1)%tag, sc_p_im_sel)
  
  ! Create the scaffold grid (used in setting up the frequency arrays)

  gr = grid_t(ml%grid(), gr_p_sel%x_i, gr_p_sel%x_o)

  ! Set up the frequency arrays

  call build_scan(ml, gr, md_p(1), os_p_sel, sc_p_re_sel, omega_re)
  call build_scan(ml, gr, md_p(1), os_p_sel, sc_p_im_sel, omega_im)

  ! Create the full grid

  gr = grid_t(ml, omega_re, gr_p_sel, md_p(1), os_p_sel)

  ! Set up the bep

  bp = nad_bep_t(ml, gr, omega_re, md_p(1), nm_p_sel, os_p_sel)

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
     
     call bp%eval_discrim(omega, discrim, status)
     $ASSERT(status == STATUS_OK,Invalid status)

     discrim_map_f(i(1),i(2)) = FRACTION(discrim)
     discrim_map_e(i(1),i(2)) = EXPONENT(discrim)

  end do

  $if ($MPI)

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

  subroutine read_scan_par (unit, axis, sc_p)

    integer, intent(in)                        :: unit
    character(*), intent(in)                   :: axis
    type(scan_par_t), allocatable, intent(out) :: sc_p(:)

    integer                         :: n_sc_p
    integer                         :: i
    real(WP)                        :: freq_min
    real(WP)                        :: freq_max
    integer                         :: n_freq
    character(LEN(sc_p%freq_units)) :: freq_units
    character(LEN(sc_p%freq_frame)) :: freq_frame
    character(LEN(sc_p%grid_type))  :: grid_type
    character(LEN(sc_p%grid_frame)) :: grid_frame
    character(LEN(sc_p%tag_list))   :: tag_list

    namelist /re_scan/ freq_min, freq_max, n_freq, freq_units, freq_frame, &
         grid_type, grid_frame, tag_list

    namelist /im_scan/ freq_min, freq_max, n_freq, freq_units, freq_frame, &
         grid_type, grid_frame, tag_list

    ! Count the number of scan namelists

    rewind(unit)

    n_sc_p = 0

    count_loop : do
       select case (axis)
       case ('re')
          read(unit, NML=re_scan, END=100)
       case ('im')
          read(unit, NML=im_scan, END=100)
       case default
          $ABORT(Invalid axis)
       end select
       n_sc_p = n_sc_p + 1
    end do count_loop

100 continue

    ! Read scan parameters

    rewind(unit)

    allocate(sc_p(n_sc_p))

    read_loop : do i = 1, n_sc_p

       freq_min = 1._WP
       freq_max = 10._WP
       n_freq = 10
          
       freq_units = 'NONE'
       freq_frame = 'INERTIAL'

       grid_type = 'LINEAR'
       grid_frame = 'INERTIAL'

       tag_list = ''

       select case (axis)
       case ('re')
          read(unit, NML=re_scan)
       case ('im')
          read(unit, NML=im_scan)
       case default
          $ABORT(Invalid axis)
       end select

       ! Initialize the scan_par

       sc_p(i) = scan_par_t(freq_min=freq_min, &
                            freq_max=freq_max, &
                            n_freq=n_freq, &
                            freq_units=freq_units, &
                            freq_frame=freq_frame, &
                            grid_type=grid_type, &
                            grid_frame=grid_frame, &
                            tag_list=tag_list)

    end do read_loop

    ! Finish

    return

  end subroutine read_scan_par

  !****

  subroutine read_out_par (unit, map_filename)

    integer, intent(in)       :: unit
    character(*), intent(out) :: map_filename

    namelist /output/ map_filename

    ! Read output parameters

    rewind(unit)

    read(unit, NML=output)

    ! Finish
    
    return

  end subroutine read_out_par
          
end program gyre_map
