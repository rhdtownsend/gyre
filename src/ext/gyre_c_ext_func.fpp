! Module   : gyre_c_ext_func
! Purpose  : monovariate functions with extended-range arithmetic (complex)
!
! Copyright 2013-2017 Rich Townsend
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

module gyre_c_ext_func

  ! Uses

  use core_kinds

  use gyre_ext

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: c_ext_func_t
   contains
     procedure(eval), deferred :: eval
  end type c_ext_func_t

  ! Interfaces

  abstract interface
     subroutine eval (this, cx, f_cx, status)
      use gyre_c_ext
      import c_ext_func_t
      class(c_ext_func_t), intent(inout) :: this
      type(c_ext_t), intent(in)          :: cx
      type(c_ext_t), intent(out)         :: f_cx
      integer, intent(out)               :: status
    end subroutine eval
  end interface

  ! Access specifiers

  private

  public :: c_ext_func_t

end module gyre_c_ext_func
