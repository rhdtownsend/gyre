! Module   : gyre_gridpar
! Purpose  : grid parameters
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

module gyre_gridpar

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: gridpar_t
     private
     real(WP), public          :: alpha_osc = 0._WP
     real(WP), public          :: alpha_exp = 0._WP
     real(WP), public          :: omega_a = 0._WP
     real(WP), public          :: omega_b = 0._WP
     real(WP), public          :: s = 0
     integer, public           :: n = 0
     character(LEN=64), public :: op_type = ''
  end type gridpar_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_gp_0
     module procedure bcast_gp_1
  end interface bcast

  interface bcast_alloc
     module procedure bcast_alloc_gp_0
     module procedure bcast_alloc_gp_1
  end interface bcast_alloc

  $endif

  ! Access specifiers

  private

  public :: gridpar_t
  $if($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  $if($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_gp_$RANK (gp, root_rank)

    type(gridpar_t), intent(inout) :: gp$ARRAY_SPEC($RANK)
    integer, intent(in)            :: root_rank

    ! Broadcast the gridpar

    call bcast(gp%alpha_osc, root_rank)
    call bcast(gp%alpha_exp, root_rank)

    call bcast(gp%omega_a, root_rank)
    call bcast(gp%omega_b, root_rank)

    call bcast(gp%s, root_rank)

    call bcast(gp%n, root_rank)

    call bcast(gp%op_type, root_rank)

    ! Finish

    return

  end subroutine bcast_gp_$RANK

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $define $BCAST_ALLOC $sub

  $local $RANK $1

  subroutine bcast_alloc_gp_$RANK (gp, root_rank)

    type(gridpar_t), allocatable, intent(inout) :: gp$ARRAY_SPEC($RANK)
    integer, intent(in)                         :: root_rank

    logical :: alloc
    $if($RANK > 0)
    integer :: s(SIZE(SHAPE(gp)))
    $endif

    ! Deallocate the gridpar on non-root processors

    if(MPI_RANK /= root_rank .AND. ALLOCATED(gp)) then
       deallocate(gp)
    endif

    ! Check if the gridpar is allocated on the root processor

    if(MPI_RANK == root_rank) alloc = ALLOCATED(gp)
    call bcast(alloc, root_rank)

    if(alloc) then

       ! Broadcast the shape

       $if($RANK > 0)

       if(MPI_RANK == root_rank) s = SHAPE(gp)

       call bcast(s, root_rank)

       $endif

       ! Allocate the buffer

       $if($RANK > 0)

       if(MPI_RANK /= root_rank) then
          if(ALLOCATED(gp)) deallocate(gp)
          allocate(gp($ARRAY_EXPAND(s,$RANK)))
       endif

       $else

       if(MPI_RANK /= root_rank) then
          if(ALLOCATED(gp)) deallocate(gp)
          allocate(gp)
       endif
       
       $endif

       ! Broadcast the gridpar

       call bcast(gp, root_rank)

    endif

    ! Finish

    return

  end subroutine bcast_alloc_gp_$RANK

  $endsub

  $BCAST_ALLOC(0)
  $BCAST_ALLOC(1)

  $endif

end module gyre_gridpar
