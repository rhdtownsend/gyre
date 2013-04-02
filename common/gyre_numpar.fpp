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

module gyre_numpar

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: numpar_t
     private
     integer, public           :: n_iter_max
     character(LEN=64), public :: ivp_solver_type
   contains
     private
     procedure, public :: init
  end type numpar_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_np
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: numpar_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, n_iter_max, ivp_solver_type)

    class(numpar_t), intent(out) :: this
    integer, intent(in)          :: n_iter_max
    character(LEN=*), intent(in) :: ivp_solver_type

    ! Initialize the numpar

    this%n_iter_max = n_iter_max

    this%ivp_solver_type = ivp_solver_type
    
    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_np (np, root_rank)

    type(numpar_t), intent(inout) :: np
    integer, intent(in)           :: root_rank

    ! Broadcast the numpar

    call bcast(np%n_iter_max, root_rank)
    call bcast(np%ivp_solver_type, root_rank)

    ! Finish

    return

  end subroutine bcast_np

  $endif

end module gyre_numpar
