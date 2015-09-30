! Program  : build_poly
! Purpose  : build a polytrope with a single discontinuity
!
! Copyright 2015 Rich Townsend
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

program build_poly_disc

  ! Uses

  use core_kinds
  use gyre_constants
  use core_hgroup
  use core_system

  use gyre_lane_emden_disc

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  real(WP)                  :: n_poly
  real(WP)                  :: Gamma_1
  real(WP)                  :: Theta_d
  real(WP)                  :: Delta_d
  real(WP)                  :: dxi
  real(WP)                  :: toler
  character(:), allocatable :: filename
  real(WP), allocatable     :: xi(:)
  real(WP), allocatable     :: Theta(:)
  real(WP), allocatable     :: dTheta(:)
  real(WP)                  :: xi_d
  integer                   :: n
  type(hgroup_t)            :: hg

  ! Read parameters

  $ASSERT(n_arg() == 7,Syntax: build_poly_disc n_poly Gamma_1 Theta_d Delta_d dxi toler filename)

  call get_arg(1, n_poly)
  call get_arg(2, Gamma_1)
  call get_arg(3, Theta_d)
  call get_arg(4, Delta_d)
  call get_arg(5, dxi)
  call get_arg(6, toler)
  call get_arg(7, filename)

  ! Solve the discontinuous Lane-Emden equation

  call solve_lane_emden(n_poly, Theta_d, Delta_d, dxi, toler, xi, Theta, dTheta, xi_d)

  n = SIZE(xi)

  ! Write the model

  hg = hgroup_t(filename, CREATE_FILE)

  call write_attr(hg, 'n', n)

  call write_attr(hg, 'n_poly', n_poly)
  call write_attr(hg, 'Gamma_1', Gamma_1)
  call write_attr(hg, 'Theta_d', Theta_d)
  call write_attr(hg, 'Delta_d', Delta_d)
  call write_attr(hg, 'xi_d', xi_d)

  call write_dset(hg, 'xi', xi)
  call write_dset(hg, 'Theta', Theta)
  call write_dset(hg, 'dTheta', dTheta)
  call write_dset(hg, 'Omega_rot', SPREAD(0._WP, DIM=1, NCOPIES=n))

  call hg%final()

  ! Finish

end program build_poly_disc
