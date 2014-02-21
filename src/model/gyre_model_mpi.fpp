! Module   : gyre_model_mpi
! Purpose  : MPI support for gyre_model
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

module gyre_model_mpi

  ! Uses

  use core_parallel

  use gyre_model
  use gyre_model_evol
  use gyre_model_poly
  use gyre_model_hom

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  $if ($MPI)
  interface bcast_alloc
     module procedure bcast_alloc_
  end interface bcast_alloc
  $endif

  ! Access specifiers

  private

  $if ($MPI)
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  $if ($MPI)

  subroutine bcast_alloc_ (cf, root_rank)

    class(model_t), allocatable, intent(inout) :: cf
    integer, intent(in)                        :: root_rank

    integer, parameter :: EVOL_TYPE = 1
    integer, parameter :: POLY_TYPE = 2
    integer, parameter :: HOM_TYPE = 3

    logical :: alloc
    integer :: type

    ! Deallocate the model on non-root processors

    if(MPI_RANK /= root_rank .AND. ALLOCATED(cf)) then
       deallocate(cf)
    endif

    ! Check if the model is allocated on the root processor

    if(MPI_RANK == root_rank) alloc = ALLOCATED(cf)
    call bcast(alloc, root_rank)

    if(alloc) then

       ! Broadcast the dynamic type

       if(MPI_RANK == root_rank) then

          select type (cf)
          type is (model_evol_t)
             type = EVOL_TYPE
          type is (model_poly_t)
             type = POLY_TYPE
          type is (model_hom_t)
             type = HOM_TYPE
          class default
             $ABORT(Unsupported type)
          end select
          
       end if

       call bcast(type, root_rank)

       ! Allocate the model

       if(MPI_RANK /= root_rank) then
          select case (type)
          case (EVOL_TYPE)
             allocate(model_evol_t::cf)
          case (POLY_TYPE)
             allocate(model_poly_t::cf)
          case(HOM_TYPE)
             allocate(model_hom_t::cf)
          case default
             $ABORT(Unsupported type)
          end select
       endif

       ! Broadcast the model

       select type (cf)
       type is (model_evol_t)
          call bcast(cf, root_rank)
       type is (model_poly_t)
          call bcast(cf, root_rank)
       type is (model_hom_t)
          call bcast(cf, root_rank)
       class default
          $ABORT(Unsupported type)
       end select
          
    endif

    ! Finish

  end subroutine bcast_alloc_

  $endif

end module gyre_model_mpi
