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

  use gyre_version
  use gyre_base_coeffs
  use gyre_therm_coeffs
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_ad_bvp
  use gyre_nad_bvp
  use gyre_ad_search
  use gyre_nad_search
  use gyre_eigfunc
  use gyre_frontend

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
  real(WP), allocatable              :: omega(:)
  type(gridpar_t), allocatable       :: shoot_gp(:)
  type(gridpar_t), allocatable       :: recon_gp(:)
  type(ad_bvp_t)                     :: ad_bp
  type(nad_bvp_t)                    :: nad_bp
  type(eigfunc_t), allocatable       :: ad_ef(:)
  type(eigfunc_t), allocatable       :: nad_ef(:)

  ! Initialize

  call init_parallel()

  call write_header('gyre_nad ['//TRIM(version)//']', '=')

  if(MPI_RANK == 0) then

     write(OUTPUT_UNIT, 100) 'Compler         : ', COMPILER_VERSION()
     write(OUTPUT_UNIT, 100) 'Compler options : ', COMPILER_OPTIONS()
100  format(3A)

     write(OUTPUT_UNIT, 110) 'OpenMP Threads  : ', OMP_SIZE_MAX
     write(OUTPUT_UNIT, 110) 'MPI Processors  : ', MPI_SIZE
110  format(A,I0)
     
     call parse_args(filename)
     
     open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  endif

  call write_header('Initialization', '=')

  if(MPI_RANK == 0) then

     call init_coeffs(unit, x_bc, bc, tc)

     $ASSERT(ALLOCATED(tc),No therm_coeffs data)

     call init_oscpar(unit, op)
     call init_numpar(unit, np)
     call init_scan(unit, bc, omega)
     call init_shoot_grid(unit, MINVAL(omega), MAXVAL(omega), shoot_gp)
     call init_recon_grid(unit, 0._WP, 0._WP, recon_gp)

     call ad_bp%init(bc, tc, op, np, shoot_gp, recon_gp, x_bc)
     call nad_bp%init(bc, tc, op, np, shoot_gp, recon_gp, x_bc)

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

  close(unit)

  call final_parallel()

end program gyre_nad
