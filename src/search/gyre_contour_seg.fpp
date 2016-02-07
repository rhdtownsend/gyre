! Module   : gyre_contour_seg
! Purpose  : contour segments
!
! Copyright 2015-2016 Rich Townsend
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

module gyre_contour_seg

  ! Uses

  use core_kinds

  use gyre_ext

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: contour_seg_t
     type(r_ext_t) :: x(2)
     type(r_ext_t) :: y(2)
   contains
     private
  end type contour_seg_t

  ! Interfaces

  interface contour_seg_t
     module procedure contour_seg_t_
  end interface contour_seg_t

  ! Access specifiers

  private

  public :: contour_seg_t

  ! Procedures

contains

  function contour_seg_t_ (x, y) result (cs)

    type(r_ext_t), intent(in) :: x(:)
    type(r_ext_t), intent(in) :: y(:)
    type(contour_seg_t)       :: cs

    $CHECK_BOUNDS(SIZE(x),2)
    $CHECK_BOUNDS(SIZE(y),2)

    ! Construct the contour_seg_t

    cs%x = x
    cs%y = y

    ! Finish

    return

  end function contour_seg_t_

end module gyre_contour_seg
