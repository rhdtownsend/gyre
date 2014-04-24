! Program  : build_poly_hom
! Purpose  : build homogeneous polytrope
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

program build_poly_hom

  ! Uses

  use core_kinds
  use gyre_constants
  use core_hgroup

  use gyre_grid

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  real(WP)                    :: Gamma_1
  real(WP)                    :: Omega_rot
  character(LEN=256)          :: grid_type
  real(WP)                    :: s
  integer                     :: n
  character(LEN=FILENAME_LEN) :: file
  real(WP), allocatable       :: x(:)
  real(WP), allocatable       :: xi(:)
  real(WP), allocatable       :: Theta(:)
  real(WP), allocatable       :: dTheta(:)
  type(hgroup_t)              :: hg

  namelist /poly/ Gamma_1, Omega_rot
  namelist /grid/ grid_type, s, n
  namelist /output/ file

  ! Read parameters

  Gamma_1 = 5._WP/3._WP
  Omega_rot = 0._WP

  read(INPUT_UNIT, NML=poly)
  
  grid_type = 'GEOM'

  s = 100._WP
  n = 100

  read(INPUT_UNIT, NML=grid)

  read(INPUT_UNIT, NML=output)

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

  hg = hgroup_t(file, CREATE_FILE)

  call write_attr(hg, 'n', n)

  call write_attr(hg, 'n_poly', 0._WP)
  call write_attr(hg, 'Gamma_1', Gamma_1)

  call write_dset(hg, 'xi', xi)
  call write_dset(hg, 'Theta', Theta)
  call write_dset(hg, 'dTheta', dTheta)
  call write_dset(hg, 'Omega_rot', SPREAD(Omega_rot, DIM=1, NCOPIES=n))

  call hg%final()

  ! Finish

end program build_poly_hom
