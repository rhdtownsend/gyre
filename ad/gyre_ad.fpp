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
  use core_memory

  use gyre_version
  use gyre_base_coeffs
  use gyre_therm_coeffs
  use gyre_oscpar
  use gyre_numpar
  use gyre_gridpar
  use gyre_ad_bvp
  use gyre_ad_search
  use gyre_mode
  use gyre_input
  use gyre_input_coeffs
  use gyre_output
  use gyre_util

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
  type(ad_bvp_t)                     :: bp
  type(mode_t), allocatable          :: md(:)

  ! Initialize

  call init_parallel()

  call write_header('gyre_ad ['//TRIM(version)//']', '=')

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

     call read_coeffs(unit, x_bc, bc, tc)
     call read_oscpar(unit, op)
     call read_numpar(unit, np)
     call read_shoot_gridpar(unit, shoot_gp)
     call read_recon_gridpar(unit, recon_gp)
     call read_scanpar(unit, bc, op, shoot_gp, x_bc, omega)

     call bp%init(bc, tc, op, np, shoot_gp, recon_gp, x_bc)

  end if

  $if($MPI)
  call bcast_alloc(omega, 0)
  call bcast(bp, 0)
  $endif

  ! Find modes

  call ad_scan_search(bp, omega, md)

  ! Write output
 
  if(MPI_RANK == 0) then
     call write_data(unit, md)
  endif

  ! Finish

  close(unit)

  call final_parallel()

end program gyre_ad
