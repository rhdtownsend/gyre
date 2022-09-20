! Program  : poly_to_parfait
! Purpose  : convert a polytrope to PARFAIT format
!
! Copyright 2022 Rich Townsend & The GYRE Team
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

program poly_to_parfait

  ! Uses

  use core_kinds
  use core_hgroup
  use core_system

  use gyre_constants
  use gyre_grid
  use gyre_math
  use gyre_model
  use gyre_model_par
  use gyre_point
  use gyre_poly_file
  use gyre_poly_model

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameters

  integer, parameter :: IVERS = 1300
  integer, parameter :: ICONST = 15
  integer, parameter :: IVAR = 40

  ! Variables

  character(:), allocatable :: in_filename
  character(:), allocatable :: out_filename
  integer                   :: n_zone

  type(model_par_t)       :: ml_p
  class(model_t), pointer :: ml
  type(grid_t)            :: gr
  type(point_t)           :: pt_i
  type(point_t)           :: pt_o
  integer                 :: n
  integer                 :: k
  real(WP), allocatable   :: x(:)
  real(WP), allocatable   :: m(:)
  real(WP), allocatable   :: Gamma_1(:)
  type(hgroup_t)          :: hg

  ! Read parameters

  $ASSERT(n_arg() == 3,Syntax: poly_to_parfait in_filename out_filename n_zone)

  call get_arg(1, in_filename)
  call get_arg(2, out_filename)
  call get_arg(3, n_zone)
  
  ! Initialize

  call init_math()

  ! Read the POLY file

  ml_p%file = TRIM(in_filename)

  call read_poly_model(ml_p, ml)

  ! Extract the grid

  gr = ml%grid()

  pt_i = gr%pt_i()
  pt_o = gr%pt_o()

  ! Set up dimensionless radii and masses

  n = n_zone + 1

  x = [((pt_i%x*(n-k) + pt_o%x*(k-1))/(n-1), k=1,n)]

  allocate(m(n))

  do k = 1, n
     m(k) = x(k)**3/ml%coeff(I_C_1, gr%pt_x(x(k)))
  end do

  ! Set up Gamma_1 data. Note that currently this uses the Gamma_1 of
  ! the inner-most segment. It should calculate an appropriately
  ! weighted average instead

  allocate(Gamma_1(n_zone))

  Gamma_1 = ml%coeff(I_GAMMA_1, pt_i)

  ! Write out the PARFAIT file

  hg = hgroup_t(out_filename, CREATE_FILE)

  call write_dset(hg, 'x', x)
  call write_dset(hg, 'm', m)

  call write_dset(hg, 'Gamma_1', Gamma_1)

  call hg%final()

  ! Finish

end program poly_to_parfait
