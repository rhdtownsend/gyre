! Module   : gyre_ivp_colloc
! Purpose  : solve initial-value problems (interface, for collocation schemes)
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

module gyre_ivp_colloc

  ! Uses

  use core_kinds
  use core_linalg

  use gyre_jacobian
  use gyre_ivp

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract, extends (ivp_t) :: ivp_colloc_t
  end type ivp_colloc_t

  ! Access specifiers

  private

  public :: ivp_colloc_t

end module gyre_ivp_colloc
