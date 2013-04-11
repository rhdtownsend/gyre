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

  use core_kinds
  use core_constants
  use core_parallel
  use core_hgroup

  use gyre_base_coeffs
  use gyre_therm_coeffs
  use gyre_oscpar
  use gyre_gridpar
  use gyre_numpar
  use gyre_ad_bvp
  use gyre_ad_search
  use gyre_nad_bvp
  use gyre_nad_search
  use gyre_eigfunc
  use gyre_frontend
  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  integer                            :: unit
  real(WP), allocatable              :: x_bc(:)
  class(base_coeffs_t), allocatable  :: bc
  class(therm_coeffs_t), allocatable :: tc
  type(oscpar_t)                     :: op
  type(numpar_t)                     :: np
  type(ad_bvp_t)                     :: ad_bp
  type(nad_bvp_t)                    :: nad_bp
  real(WP), allocatable              :: omega(:)
  type(eigfunc_t), allocatable       :: ad_ef(:)
  type(eigfunc_t), allocatable       :: nad_ef(:)

  ! Initialize

  call init_parallel()

  if(MPI_RANK == 0) then

     write(OUTPUT_UNIT, '(A)') 'gyre_nad [hg]'
     write(OUTPUT_UNIT, '(A,2X,I0)') 'OpenMP Threads :', OMP_SIZE_MAX
     write(OUTPUT_UNIT, '(A,2X,I0)') 'MPI Processors :', MPI_SIZE
     
     call open_input(unit)

  endif

  call write_header('Initialization', '=')

  if(MPI_RANK == 0) then

     call init_coeffs(unit, x_bc, bc, tc)

     $ASSERT(ALLOCATED(tc),No therm_coeffs data)

     call init_oscpar(unit, op)
     call init_numpar(unit, np)
     call init_scan(unit, bc, omega)
     call init_bvp(unit, x_bc, bc, tc, op, np, omega, ad_bp, nad_bp)

  endif

  $if($MPI)
  call bcast_alloc(omega, 0)
  call bcast(ad_bp, 0)
  call bcast(nad_bp, 0)
  $endif

  ! Search for modes

  call ad_scan_search(ad_bp, omega, ad_ef)
  call nad_prox_search(nad_bp, ad_ef, nad_ef)

  ! Write output
 
  if(MPI_RANK == 0) then
     call write_data(unit, nad_ef, bc)
  endif

  ! Finish

  call final_parallel()

contains

  subroutine init_bvp (unit, x_bc, bc, tc, op, np, omega, ad_bp, nad_bp)

    integer, intent(in)                            :: unit
    real(WP), intent(in), allocatable              :: x_bc(:)
    class(base_coeffs_t), intent(in)               :: bc
    class(therm_coeffs_t), allocatable, intent(in) :: tc
    type(oscpar_t), intent(in)                     :: op
    type(numpar_t), intent(in)                     :: np
    real(WP), intent(in)                           :: omega(:)
    type(ad_bvp_t), intent(out)                    :: ad_bp
    type(nad_bvp_t), intent(out)                   :: nad_bp

    character(LEN=256)    :: grid_type
    real(WP)              :: alpha_osc
    real(WP)              :: alpha_exp
    integer               :: n_center
    integer               :: n_floor
    real(WP)              :: s
    integer               :: n_grid
    integer               :: dn(SIZE(x_bc)-1)
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
    read(unit, NML=shoot_grid, END=900)

    ! Initialize the shooting grid

    select case(grid_type)
    case('GEOM')
       call build_geom_grid(s, n_grid, x_sh)
    case('LOG')
       call build_log_grid(s, n_grid, x_sh)
    case('INHERIT')
       $ASSERT(ALLOCATED(x_bc),No input grid)
       x_sh = x_bc
    case('DISPERSION')
       $ASSERT(ALLOCATED(x_bc),No input grid)
       dn = 0
       do i = 1,SIZE(omega)
          call plan_dispersion_grid(x_bc, bc, CMPLX(omega(i), KIND=WP), op, alpha_osc, alpha_exp, n_center, n_floor, dn)
       enddo
       call build_oversamp_grid(x_bc, dn, x_sh)
    case default
       $ABORT(Invalid grid_type)
    end select

    ! Read recon grid parameters

    alpha_osc = 0._WP
    alpha_exp = 0._WP
       
    n_center = 0
    n_floor = 0

    rewind(unit)
    read(unit, NML=recon_grid, END=910)
       
    call gp%init(alpha_osc, alpha_exp, n_center, n_floor, 0._WP, 0, 'DISP')

    ! Initialize the bvps

    call ad_bp%init(bc, tc, op, gp, np, x_sh)
    call nad_bp%init(bc, tc, op, gp, np, x_sh)

    ! Finish

    return

    ! Jump-in points for end-of-file

900 continue

    $ABORT(No &shoot_grid namelist in input file)

910 continue

    $ABORT(No &recon_grid namelist in input file)

  end subroutine init_bvp

end program gyre_nad
