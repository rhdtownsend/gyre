! Module   : gyre_therm_coeffs_mpi
! Purpose  : MPI support for gyre_therm_coeffs
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

module gyre_therm_coeffs_mpi

  ! Uses

  use core_parallel

  use gyre_therm_coeffs
  use gyre_evol_therm_coeffs

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  $if($MPI)
  interface bcast_alloc
     module procedure bcast_alloc_tc
  end interface bcast_alloc
  $endif

  ! Access specifiers

  private

  $if($MPI)
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  $if($MPI)

  subroutine bcast_alloc_tc (tc, root_rank)

    class(therm_coeffs_t), allocatable, intent(inout) :: tc
    integer, intent(in)                               :: root_rank

    integer, parameter :: EVOL_TYPE = 1

    logical :: alloc
    integer :: type

    ! Deallocate the therm_coeffs on non-root processors

    if(MPI_RANK /= root_rank .AND. ALLOCATED(tc)) then
       deallocate(tc)
    endif

    ! Check if the therm_coeffs is allocated on the root processor

    if(MPI_RANK == root_rank) alloc = ALLOCATED(tc)
    call bcast(alloc, root_rank)

    if(alloc) then

       ! Broadcast the dynamic type

       if(MPI_RANK == root_rank) then

          select type (tc)
          type is (evol_therm_coeffs_t)
             type = EVOL_TYPE
          class default
             $ABORT(Unsupported type)
          end select
          
       end if

       call bcast(type, root_rank)

       ! Allocate the therm_coeffs

       if(MPI_RANK /= root_rank) then
          select case (type)
          case (EVOL_TYPE)
             allocate(evol_therm_coeffs_t::tc)
          case default
             $ABORT(Unsupported type)
          end select
       endif

       ! Broadcast the therm_coeffs

       call tc%bcast(root_rank)

    endif

    ! Finish

  end subroutine bcast_alloc_tc

  $endif

end module gyre_therm_coeffs_mpi
