! Module   : gyre_seg
! Purpose  : segment (interface)
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

module gyre_seg

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: seg_t
   contains
     private
     procedure(x_min_), deferred, public :: x_min
     procedure(x_max_), deferred, public :: x_max
  end type seg_t

  ! Interfaces

  abstract interface

     function x_min_ (this) result (x_min)
       use gyre_kinds
       import seg_t
       class(seg_t), intent(in) :: this
       real(WP), intent(out)    :: x_min
     end function x_min_

     function x_max_ (this) result (x_max)
       use gyre_kinds
       import seg_t
       class(seg_t), intent(in) :: this
       real(WP), intent(out)    :: x_max
     end function x_max_

  end interface

  ! Access specifiers

  private

  public :: seg_t

end module gyre_seg
