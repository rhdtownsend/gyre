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

  use core_kinds
  use core_constants
  use core_parallel
  use core_hgroup
  use core_order

  use gyre_mech_coeffs
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_ad_bvp
  use gyre_ad_search
  use gyre_mode
  use gyre_frontend
  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  integer                           :: unit
  real(WP), allocatable             :: x_mc(:)
  class(mech_coeffs_t), allocatable :: mc
  type(oscpar_t)                    :: op
  type(numpar_t)                    :: np
  type(ad_bvp_t)                    :: bp 
  real(WP), allocatable             :: omega(:)
  type(mode_t), allocatable         :: md(:)

  ! Initialize

  call init_parallel()

  if(MPI_RANK == 0) then

     write(OUTPUT_UNIT, '(A)') 'gyre_ad [hg]'
     write(OUTPUT_UNIT, '(A,2X,I0)') 'OpenMP Threads :', OMP_SIZE_MAX
     write(OUTPUT_UNIT, '(A,2X,I0)') 'MPI Processors :', MPI_SIZE
     
     call open_input(unit)

  endif

  call write_header('Initialization', '=')

  if(MPI_RANK == 0) then

     call init_coeffs(unit, x_mc, mc)
     call init_oscpar(unit, op)
     call init_numpar(unit, np)
     call init_scan(unit, mc, omega)
     call init_bvp(unit, x_mc, mc, op, np, omega, bp)

  end if

  $if($MPI)
  call bcast_alloc(omega, 0)
  call bcast(bp, 0)
  $endif

  ! Find modes

  call ad_scan_search(bp, omega, md)

  ! Write output
 
  if(MPI_RANK == 0) then
     call write_eigdata(unit, bp, md)
  endif

  ! Finish

  call final_parallel()

contains

!****

  subroutine init_bvp (unit, x_mc, mc, op, np, omega, bp)

    integer, intent(in)                      :: unit
    real(WP), intent(in), allocatable        :: x_mc(:)
    class(mech_coeffs_t), intent(in), target :: mc
    type(oscpar_t), intent(in)               :: op
    type(numpar_t), intent(in)               :: np
    real(WP), intent(in)                     :: omega(:)
    type(ad_bvp_t), intent(out)              :: bp

    character(LEN=256)    :: grid_type
    real(WP)              :: alpha_osc
    real(WP)              :: alpha_exp
    integer               :: n_center
    integer               :: n_floor
    real(WP)              :: s
    integer               :: n_grid
    integer               :: dn(SIZE(x_mc)-1)
    integer               :: i
    real(WP), allocatable :: x_sh(:)
    type(gridpar_t)       :: gp

    namelist /shoot_grid/ grid_type, alpha_osc, alpha_exp, &
         n_center, n_floor, s, n_grid

    namelist /recon_grid/ alpha_osc, alpha_exp, n_center, n_floor

    ! Read shooting grid parameters

    grid_type = 'DISPERSION'

    alpha_osc = 0._WP
    alpha_exp = 0._WP
       
    n_center = 0
    n_floor = 0

    s = 100._WP
    n_grid = 100

    rewind(unit)
    read(unit, NML=shoot_grid, END=100)

100 continue

    ! Initialize the shooting grid

    select case(grid_type)
    case('GEOM')
       call build_geom_grid(s, n_grid, x_sh)
    case('LOG')
       call build_log_grid(s, n_grid, x_sh)
    case('INHERIT')
       $ASSERT(ALLOCATED(x_mc),No input grid)
       x_sh = x_mc
    case('DISPERSION')
       $ASSERT(ALLOCATED(x_mc),No input grid)
       dn = 0
       do i = 1,SIZE(omega)
          call plan_dispersion_grid(x_mc, mc, CMPLX(omega(i), KIND=WP), op, alpha_osc, alpha_exp, n_center, n_floor, dn)
       enddo
       call build_oversamp_grid(x_mc, dn, x_sh)
    case default
       $ABORT(Invalid grid_type)
    end select

    ! Read recon grid parameters

    alpha_osc = 0._WP
    alpha_exp = 0._WP
       
    n_center = 0
    n_floor = 0

    rewind(unit)
    read(unit, NML=recon_grid, END=200)

200 continue

    call gp%init(alpha_osc, alpha_exp, n_center, n_floor, 0._WP, 0, 'DISP')

    ! Initialize the bvp

    call bp%init(mc, op, gp, np, x_sh)

    ! Finish

    return

  end subroutine init_bvp

!****

  subroutine write_eigdata (unit, bp, md)

    integer, intent(in)        :: unit
    type(ad_bvp_t), intent(in) :: bp
    type(mode_t), intent(in)   :: md(:)

    character(LEN=256)               :: freq_units
    character(LEN=FILENAME_LEN)      :: eigval_file
    character(LEN=FILENAME_LEN)      :: eigfunc_prefix
    integer                          :: i
    complex(WP)                      :: freq(SIZE(md))
    type(hgroup_t)                   :: hg
    character(LEN=FILENAME_LEN)      :: eigfunc_file

    namelist /output/ freq_units, eigval_file, eigfunc_prefix

    ! Read output parameters

    freq_units = 'NONE'

    eigval_file = ''
    eigfunc_prefix = ''

    rewind(unit)
    read(unit, NML=output, END=100)

100 continue

    ! Write eigenvalues

    freq_loop : do i = 1,SIZE(md)
       freq(i) = bp%mc%conv_freq(md(i)%omega, 'NONE', freq_units)
    end do freq_loop

    if(eigval_file /= '') then

       call hg%init(eigval_file, CREATE_FILE)
          
       call write_attr(hg, 'n_md', SIZE(md))

       call write_attr(hg, 'n', bp%n)
       call write_attr(hg, 'n_e', bp%n_e)

       call write_attr(hg, 'l', bp%op%l)

       call write_dset(hg, 'n_p', md%n_p)
       call write_dset(hg, 'n_g', md%n_g)

       call write_dset(hg, 'freq', freq)
       call write_attr(hg, 'freq_units', freq_units)

       call write_dset(hg, 'E', md%E)

       call hg%final()

    end if

    ! Write eigenfunctions

    if(eigfunc_prefix /= '') then

       mode_loop : do i = 1,SIZE(md)

          write(eigfunc_file, 200) TRIM(eigfunc_prefix), i, '.h5'
200       format(A,I4.4,A)

          call md(i)%write(eigfunc_file, 'GYRE')

       end do mode_loop

    end if

    ! Finish

    return

  end subroutine write_eigdata

end program gyre_ad
