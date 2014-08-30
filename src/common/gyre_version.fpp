! Program  : gyre_version
! Purpose  : versioning info
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

module gyre_version

  ! Uses

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameters

  character(LEN=*), parameter :: version = '3.2.1'

  ! Access specifiers

  private

  public :: version

end module gyre_version
