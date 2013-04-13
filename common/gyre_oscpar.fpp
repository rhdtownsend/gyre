! Module   : gyre_oscpar
! Purpose  : oscillation parameters
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

module gyre_oscpar

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: oscpar_t
     private
     integer, public           :: l = 0
     character(LEN=64), public :: outer_bound_type = ''
  end type oscpar_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_op
  end interface bcast

  $endif

 ! Access specifiers

  private

  public :: oscpar_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  $if($MPI)

  subroutine bcast_op (op, root_rank)

    type(oscpar_t), intent(inout) :: op
    integer, intent(in)           :: root_rank

    ! Broadcast the oscpar

    call bcast(op%l, root_rank)
    call bcast(op%outer_bound_type, root_rank)

    ! Finish

    return

  end subroutine bcast_op

  $endif

end module gyre_oscpar
