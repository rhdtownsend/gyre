! Module   : gyre_c_extfunc
! Purpose  : monovariate functions with extended-range arithmetic (complex)
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

module gyre_c_extfunc

  ! Uses

  use core_kinds

  use gyre_ext

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: c_extfunc_t
   contains
     private
     procedure(eval_), deferred, public :: eval
  end type c_extfunc_t

  ! Interfaces

  abstract interface
    function eval_ (this, cx) result (f_cx)
      use gyre_c_ext
      import c_extfunc_t
      class(c_extfunc_t), intent(inout) :: this
      type(c_ext_t), intent(in)         :: cx
      type(c_ext_t)                     :: f_cx
    end function eval_
  end interface

  ! Access specifiers

  private

  public :: c_extfunc_t

end module gyre_c_extfunc
