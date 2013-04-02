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
     real(WP), public          :: lambda_0
     integer, public           :: l
     character(LEN=64), public :: outer_bound_type
   contains
     private
     procedure, public :: init
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

  subroutine init (this, l, outer_bound_type)

    class(oscpar_t), intent(out) :: this
    integer, intent(in)          :: l
    character(LEN=*), intent(in) :: outer_bound_type

    ! Initialize the oscpar

    this%lambda_0 = l - 2._WP
    this%l = l

    this%outer_bound_type = outer_bound_type

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_op (op, root_rank)

    type(oscpar_t), intent(inout) :: op
    integer, intent(in)           :: root_rank

    integer                                 :: l
    character(LEN=LEN(op%outer_bound_type)) :: outer_bound_type

    ! Broadcast the oscpar

    if(MPI_RANK == root_rank) then

       call bcast(op%l, root_rank)
       call bcast(op%outer_bound_type, root_rank)

    else

       call bcast(l, root_rank)
       call bcast(outer_bound_type, root_rank)

       call op%init(l, outer_bound_type)

    endif

    ! Finish

    return

  end subroutine bcast_op

  $endif

end module gyre_oscpar
