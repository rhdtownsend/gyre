! Program  : build_tar_fit
! Purpose  : build tar_fit_t types
!
! Copyright 2016 Rich Townsend
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

program build_tar_fit

  ! Uses

  use core_kinds
  use gyre_constants
  use core_hgroup
  use core_system

  use gyre_tar_fit

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  integer                   :: m
  integer                   :: k
  real(WP)                  :: cheb_tol
  character(:), allocatable :: filename

  type(tar_fit_t) :: tf
  type(hgroup_t)  :: hg

  ! Read parameters

  $ASSERT(n_arg() == 4,Syntax: build_tar_fit m k cheb_tol filename)

  call get_arg(1, m)
  call get_arg(2, k)
  call get_arg(3, cheb_tol)
  call get_arg(4, filename)

  ! Construct the tar_fit_t

  tf = tar_fit_t(m, k, cheb_tol)

  ! Write it out

  hg = hgroup_t(filename, CREATE_FILE)
  call write(hg, tf)
  call hg%final()

  ! Finish

end program build_tar_fit
