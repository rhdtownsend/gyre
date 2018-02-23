! Module   : gyre_grid_spec
! Purpose  : specification for grid construction
!
! Copyright 2018 Rich Townsend
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
$include 'core_parallel.inc'

module gyre_grid_spec

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_context

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: grid_spec_t
     type(context_t), pointer :: cx
     real(WP), allocatable    :: omega(:)
  end type grid_spec_t

  ! Access specifiers

  private

  public :: grid_spec_t

end module gyre_grid_spec
