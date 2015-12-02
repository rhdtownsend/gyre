! Module   : gyre_model
! Purpose  : stellar model
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

module gyre_model

  ! Uses

  use core_kinds

  use gyre_model_seg

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: model_t
     type(model_seg_t), allocatable :: ms(:)
     integer                        :: n_s
  end type model_t

  ! Access specifiers

  private

  public :: model_t

end module gyre_model
