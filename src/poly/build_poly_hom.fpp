! Program  : build_poly_hom
! Purpose  : build homogeneous polytrope
!
! Copyright 2013-214 Rich Townsend
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

program build_poly_hom

  ! Uses

  use core_kinds
  use gyre_constants
  use core_hgroup
  use core_system

  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  real(WP)                   :: Gamma_1
  character(:), allocatable  :: grid_type
  real(WP)                   :: s
  integer                    :: n
  character(:), allocatable  :: filename
  real(WP), allocatable      :: x(:)
  real(WP), allocatable      :: xi(:)
  real(WP), allocatable      :: Theta(:)
  real(WP), allocatable      :: dTheta(:)
  type(hgroup_t)             :: hg

  ! Read parameters

  $ASSERT(n_arg() == 5,Syntax: build_poly_hom Gamma_1 grid_type s n filename)

  call get_arg(1, Gamma_1)
  call get_arg(2, grid_type)
  call get_arg(3, s)
  call get_arg(4, n)
  call get_arg(5, filename)

  ! Set up the grid

  select case (grid_type)
  case ('UNIFORM')
     call create_uniform(n, x)
  case ('GEOM')
     call create_geom(s, n, x)
  case ('LOG')
     call create_log(s, n, x)
  case default
     $ABORT(Invalid grid_type)
  end select

  ! Calculate structure variables

  allocate(xi(n))

  allocate(Theta(n))
  allocate(dTheta(n))

  xi = SQRT(6._WP)*x

  Theta = 1._WP - xi**2/6._WP
  dTheta = -xi/3._WP

  ! Write the model

  hg = hgroup_t(filename, CREATE_FILE)

  call write_attr(hg, 'n', n)

  call write_attr(hg, 'n_poly', 0._WP)
  call write_attr(hg, 'Gamma_1', Gamma_1)

  call write_dset(hg, 'xi', xi)
  call write_dset(hg, 'Theta', Theta)
  call write_dset(hg, 'dTheta', dTheta)
  call write_dset(hg, 'Omega_rot', SPREAD(0._WP, DIM=1, NCOPIES=n))

  call hg%final()

  ! Finish

end program build_poly_hom
