! Module   : gyre_point
! Purpose  : segmented grid point
!
! Copyright 2013-2020 Rich Townsend & The GYRE Team
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

module gyre_point

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: point_t
     integer  :: s
     real(WP) :: x
   contains
     private
     procedure       :: op_eq_
     generic, public :: operator(==) => op_eq_
  end type point_t

  ! Access specifiers

  private

  public :: point_t

  ! Procedures

contains

  elemental function op_eq_ (this, that) result (eq)

    class(point_t), intent(in) :: this
    class(point_t), intent(in) :: that
    logical                    :: eq

    ! Evaluate the equality operator

    eq = this%s == that%s .AND. this%x == that%x

    ! Finish

    return

  end function op_eq_

end module gyre_point
