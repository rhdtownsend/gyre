! Module   : gyre_numpar
! Purpose  : numerics parameters
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

module gyre_numpar

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: numpar_t
     integer             :: n_iter_max
     real(WP)            :: theta_ad
     logical             :: reduce_order
     logical             :: use_banded
     logical             :: use_trad_approx
     character(LEN=64)   :: ivp_solver_type
     character(LEN=2048) :: tag_list
  end type numpar_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_np_0
     module procedure bcast_np_1
  end interface bcast

  interface bcast_alloc
     module procedure bcast_alloc_np_0
     module procedure bcast_alloc_np_1
  end interface bcast_alloc

  $endif

  ! Access specifiers

  private

  public :: numpar_t
  $if($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  $if($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_np_$RANK (np, root_rank)

    type(numpar_t), intent(inout) :: np$ARRAY_SPEC($RANK)
    integer, intent(in)           :: root_rank

    ! Broadcast the numpar

    call bcast(np%n_iter_max, root_rank)
    call bcast(np%theta_ad, root_rank)

    call bcast(np%reduce_order, root_rank)
    call bcast(np%use_banded, root_rank)
    call bcast(np%use_trad_approx, root_rank)

    call bcast(np%ivp_solver_type, root_rank)
    call bcast(np%tag_list, root_rank)

    ! Finish

    return

  end subroutine bcast_np_$RANK

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $BCAST_ALLOC(np,type(numpar_t),0)
  $BCAST_ALLOC(np,type(numpar_t),1)

  $endif

end module gyre_numpar
