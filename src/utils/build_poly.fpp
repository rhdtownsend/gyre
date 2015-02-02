! Program  : build_poly
! Purpose  : build a polytrope
!
! Copyright 2013-2014 Rich Townsend
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

  use gyre_lane_emden

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  real(WP)                  :: n_poly
  real(WP)                  :: Gamma_1
  real(WP)                  :: dxi
  real(WP)                  :: toler
  character(:), allocatable :: filename
  real(WP), allocatable     :: xi(:)
  real(WP), allocatable     :: Theta(:)
  real(WP), allocatable     :: dTheta(:)
  integer                   :: n
  type(hgroup_t)            :: hg

  ! Read parameters

  $ASSERT(n_arg() == 5,Syntax: build_poly n_poly Gamma_1 dxi toler filename)

  call get_arg(1, n_poly)
  call get_arg(2, Gamma_1)
  call get_arg(3, dxi)
  call get_arg(4, toler)
  call get_arg(5, filename)

  ! Solve the Lane-Emden equation

  call solve_lane_emden(n_poly, dxi, toler, xi, Theta, dTheta)

  n = SIZE(xi)

  ! Write the model

  hg = hgroup_t(filename, CREATE_FILE)

  call write_attr(hg, 'n', n)

  call write_attr(hg, 'n_poly', n_poly)
  call write_attr(hg, 'Gamma_1', Gamma_1)

  call write_dset(hg, 'xi', xi)
  call write_dset(hg, 'Theta', Theta)
  call write_dset(hg, 'dTheta', dTheta)
  call write_dset(hg, 'Omega_rot', SPREAD(0._WP, DIM=1, NCOPIES=n))

  call hg%final()

  ! Finish

end program build_poly
