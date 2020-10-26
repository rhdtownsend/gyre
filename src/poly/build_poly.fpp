! Program  : build_poly
! Purpose  : build a composite polytrope
!
! Copyright 2015-2020 Rich Townsend & The GYRE Team
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

  character(:), allocatable :: filename
  integer                   :: unit
  integer                   :: n_r
  real(WP), allocatable     :: n_poly(:)
  real(WP)                  :: Gamma_1
  real(WP), allocatable     :: z_b(:)
  real(WP), allocatable     :: Delta_b(:)
  real(WP)                  :: dz
  real(WP)                  :: toler
  character(FILENAME_LEN)   :: file
  real(WP), allocatable     :: z(:)
  real(WP), allocatable     :: theta(:)
  real(WP), allocatable     :: dtheta(:)
  type(hgroup_t)            :: hg

  namelist /poly/ n_r, n_poly, Gamma_1, z_b, Delta_b
  namelist /num/ dz, toler
  namelist /out/ file

  ! Read command-line arguments

  $ASSERT(n_arg() == 1,Syntax: build_poly <filename>)

  call get_arg(1, filename)

  ! Set defaults

  allocate(n_poly(D))
  allocate(z_b(D))
  allocate(Delta_b(D))

  n_r = 1
  Gamma_1 = 5._WP/3._WP

  n_poly(1) = 0._WP

  dz = 0.01
  toler = 1E-10

  ! Read parameters

  open(NEWUNIT=unit, FILE=filename, STATUS='OLD')

  rewind(unit)
  read(unit, NML=poly)

  call reallocate(n_poly, [n_r])
  call reallocate(z_b, [n_r-1])
  call reallocate(Delta_b, [n_r-1])

  rewind(unit)
  read(unit, NML=num)

  rewind(unit)
  read(unit, NML=out)

  close(unit)

  ! Solve the Lane-Emden equation

  call solve_lane_emden(n_poly, z_b, Delta_b, dz, toler, z, theta, dtheta)

  ! Write the model

  hg = hgroup_t(file, CREATE_FILE)

  call write_attr(hg, 'n', SIZE(z))
  call write_attr(hg, 'n_r', n_r)
  call write_attr(hg, 'n_poly', n_poly)
  if (n_r > 1) then
     call write_attr(hg, 'Delta_b', Delta_b)
  endif
  call write_attr(hg, 'Gamma_1', Gamma_1)

  call write_dset(hg, 'z', z)
  call write_dset(hg, 'theta', theta)
  call write_dset(hg, 'dtheta', dtheta)

  call hg%final()

  ! Finish

end program build_poly
