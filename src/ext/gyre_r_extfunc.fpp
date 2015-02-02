! Module   : gyre_r_extfunc
! Purpose  : monovariate functions with extended-range arithmetic (real)
!
! Copyright 2013-2015 Rich Townsend
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

module gyre_r_extfunc

  ! Uses

  use core_kinds

  use gyre_ext

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: r_extfunc_t
   contains
     private
     procedure(eval_), deferred, public :: eval
  end type r_extfunc_t

  ! Interfaces

  abstract interface
    function eval_ (this, rx) result (f_rx)
      use gyre_r_ext
      import r_extfunc_t
      class(r_extfunc_t), intent(inout) :: this
      type(r_ext_t), intent(in)         :: rx
      type(r_ext_t)                     :: f_rx
    end function eval_
  end interface

  ! Access specifiers

  private

  public :: r_extfunc_t

end module gyre_r_extfunc
