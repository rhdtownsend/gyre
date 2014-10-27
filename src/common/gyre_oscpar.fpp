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
$include 'core_parallel.inc'

module gyre_oscpar

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: oscpar_t
     real(WP)        :: x_ref
     character(64)   :: variables_type
     character(64)   :: outer_bound_type
     character(64)   :: inertia_norm_type
     character(2048) :: tag_list
     logical         :: reduce_order
  end type oscpar_t

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

  public :: oscpar_t
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  $if ($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_${RANK}_ (op, root_rank)

    type(oscpar_t), intent(inout) :: op$ARRAY_SPEC($RANK)
    integer, intent(in)           :: root_rank

    ! Broadcast the oscpar_t

    call bcast(op%x_ref, root_rank)

    call bcast(op%variables_type, root_rank)
    call bcast(op%outer_bound_type, root_rank)
    call bcast(op%inertia_norm_type, root_rank)
    call bcast(op%tag_list, root_rank)

    call bcast(np%reduce_order, root_rank)

    ! Finish

    return

  end subroutine bcast_${RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $BCAST_ALLOC(type(oscpar_t),0)
  $BCAST_ALLOC(type(oscpar_t),1)

  $endif

end module gyre_oscpar
