! Program  : build_poly
! Purpose  : build a polytrope, possibly with disctontinuities
!
! Copyright 2015-2017 Rich Townsend
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

program build_poly

  ! Uses

  use core_kinds
  use gyre_constants
  use core_hgroup
  use core_system
  use core_memory

  use gyre_lane_emden

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameters

  integer, parameter :: D = 128

  ! Variables

  character(:), allocatable :: in_filename
  integer                   :: unit
  integer                   :: n_d
  real(WP), allocatable     :: n_poly(:)
  real(WP)                  :: Gamma_1
  real(WP), allocatable     :: xi_d(:)
  real(WP), allocatable     :: Delta_d(:)
  real(WP)                  :: dxi
  real(WP)                  :: toler
  character(FILENAME_LEN)   :: filename
  real(WP), allocatable     :: xi(:)
  real(WP), allocatable     :: Theta(:)
  real(WP), allocatable     :: dTheta(:)
  type(hgroup_t)            :: hg

  namelist /poly/ n_d, n_poly, Gamma_1, xi_d, Delta_d
  namelist /num/ dxi, toler
  namelist /out/ filename

  ! Read parameters

  $ASSERT(n_arg() == 1,Syntax: build_poly_disc in_filename)

  call get_arg(1, in_filename)

  open(NEWUNIT=unit, FILE=in_filename, STATUS='OLD')

  allocate(n_poly(D))
  allocate(xi_d(D))
  allocate(Delta_d(D))

  n_d = 0
  Gamma_1 = 5._WP/3._WP

  rewind(unit)
  read(unit, NML=poly)

  call reallocate(n_poly, [n_d+1])
  call reallocate(xi_d, [n_d])
  call reallocate(Delta_d, [n_d])

  rewind(unit)
  read(unit, NML=num)

  rewind(unit)
  read(unit, NML=out)

  close(unit)

  ! Solve the discontinuous Lane-Emden equation

  call solve_lane_emden(n_poly, xi_d, Delta_d, dxi, toler, xi, Theta, dTheta)

  ! Write the model

  hg = hgroup_t(filename, CREATE_FILE)

  call write_attr(hg, 'n', SIZE(xi))
  call write_attr(hg, 'n_d', n_d)
  call write_attr(hg, 'n_poly', n_poly)
  if (n_d > 0) then
     call write_attr(hg, 'Delta_d', Delta_d)
  endif
  call write_attr(hg, 'Gamma_1', Gamma_1)

  call write_dset(hg, 'xi', xi)
  call write_dset(hg, 'Theta', Theta)
  call write_dset(hg, 'dTheta', dTheta)

  call hg%final()

  ! Finish

end program build_poly
