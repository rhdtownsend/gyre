! Module   : gyre_point
! Purpose  : segmented grid point
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_point

  ! Uses

  use core_kinds
  use core_parallel

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

  ! Interfaces
     
  $if ($MPI)
  interface bcast
     module procedure bcast_0_
     module procedure bcast_1_
  end interface bcast
  interface bcast_alloc
     module procedure bcast_alloc_0_
     module procedure bcast_alloc_1_
  end interface bcast_alloc
  $endif

  ! Access specifiers

  private

  public :: point_t
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

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

  !****

  $if ($MPI)

  subroutine bcast_0_ (pt, root_rank)

    type(point_t), intent(inout) :: pt
    integer, intent(in)          :: root_rank

    ! Broadcast the point_t

    call bcast(pt%s, root_rank)
    call bcast(pt%x, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  !****

  subroutine bcast_1_ (pt, root_rank)

    type(point_t), intent(inout) :: pt(:)
    integer, intent(in)          :: root_rank

    ! Broadcast the point_t

    call bcast(pt%s, root_rank)
    call bcast(pt%x, root_rank)

    ! Finish

    return

  end subroutine bcast_1_

  $BCAST_ALLOC(type(point_t),0)
  $BCAST_ALLOC(type(point_t),1)

  $endif

end module gyre_point
