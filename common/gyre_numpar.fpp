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
     integer           :: n_iter_max
     real(WP)          :: theta_ad
     character(LEN=64) :: ivp_solver_type
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

  $if($MPI)

  subroutine bcast_np (np, root_rank)

    type(numpar_t), intent(inout) :: np
    integer, intent(in)           :: root_rank

    ! Broadcast the numpar

    call bcast(np%n_iter_max, root_rank)
    call bcast(np%theta_ad, root_rank)

    call bcast(np%ivp_solver_type, root_rank)

    ! Finish

    return

  end subroutine bcast_np

  $endif

end module gyre_numpar
