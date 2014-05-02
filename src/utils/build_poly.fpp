! Program  : build_poly
! Purpose  : build a polytrope
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

program build_poly

  ! Uses

  use core_kinds
  use gyre_constants
  use core_hgroup

  use gyre_lane_emden

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  real(WP)                    :: Gamma_1
  real(WP)                    :: Omega_rot
  real(WP)                    :: n_poly
  real(WP)                    :: dxi
  real(WP)                    :: toler
  character(LEN=FILENAME_LEN) :: file
  real(WP), allocatable       :: xi(:)
  real(WP), allocatable       :: Theta(:)
  real(WP), allocatable       :: dTheta(:)
  integer                     :: n
  type(hgroup_t)              :: hg

  namelist /poly/ n_poly, Gamma_1, Omega_rot, dxi, toler
  namelist /output/ file

  ! Read parameters

  Gamma_1 = 5._WP/3._WP
  Omega_rot = 0._WP

  dxi = 1.E-3_WP
  toler = EPSILON(0._WP)

  read(INPUT_UNIT, NML=poly)
  
  read(INPUT_UNIT, NML=output)

  ! Solve the Lane-Emden equation

  call solve_lane_emden(n_poly, dxi, toler, xi, Theta, dTheta)

  n = SIZE(xi)

  ! Write the model

  hg = hgroup_t(file, CREATE_FILE)

  call write_attr(hg, 'n', n)

  call write_attr(hg, 'n_poly', n_poly)
  call write_attr(hg, 'Gamma_1', Gamma_1)

  call write_dset(hg, 'xi', xi)
  call write_dset(hg, 'Theta', Theta)
  call write_dset(hg, 'dTheta', dTheta)
  call write_dset(hg, 'Omega_rot', SPREAD(Omega_rot, DIM=1, NCOPIES=n))

  call hg%final()

  ! Finish

end program build_poly
