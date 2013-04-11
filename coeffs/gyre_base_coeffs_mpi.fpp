! Module   : gyre_base_coeffs_mpi
! Purpose  : MPI support for gyre_base_coeffs
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

module gyre_base_coeffs_mpi

  ! Uses

  use core_parallel

  use gyre_base_coeffs
  use gyre_evol_base_coeffs
  use gyre_poly_base_coeffs
  use gyre_hom_base_coeffs

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  $if($MPI)
  interface bcast_alloc
     module procedure bcast_alloc_bc
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

  subroutine bcast_alloc_bc (bc, root_rank)

    class(base_coeffs_t), allocatable, intent(inout) :: bc
    integer, intent(in)                              :: root_rank

    integer, parameter :: EVOL_TYPE = 1
    integer, parameter :: POLY_TYPE = 2
    integer, parameter :: HOM_TYPE = 3

    logical :: alloc
    integer :: type

    ! Deallocate the base_coeffs on non-root processors

    if(MPI_RANK /= root_rank .AND. ALLOCATED(bc)) then
       deallocate(bc)
    endif

    ! Check if the base_coeffs is allocated on the root processor

    if(MPI_RANK == root_rank) alloc = ALLOCATED(bc)
    call bcast(alloc, root_rank)

    if(alloc) then

       ! Broadcast the dynamic type

       if(MPI_RANK == root_rank) then

          select type (bc)
          type is (evol_base_coeffs_t)
             type = EVOL_TYPE
          type is (poly_base_coeffs_t)
             type = POLY_TYPE
          type is (hom_base_coeffs_t)
             type = HOM_TYPE
          class default
             $ABORT(Unsupported type)
          end select
          
       end if

       call bcast(type, root_rank)

       ! Allocate the base_coeffs

       if(MPI_RANK /= root_rank) then
          select case (type)
          case (EVOL_TYPE)
             allocate(evol_base_coeffs_t::bc)
          case (POLY_TYPE)
             allocate(poly_base_coeffs_t::bc)
          case(HOM_TYPE)
             allocate(hom_base_coeffs_t::bc)
          case default
             $ABORT(Unsupported type)
          end select
       endif

       ! Broadcast the base_coeffs

       select type (bc)
       type is (evol_base_coeffs_t)
          call bcast(bc, root_rank)
       type is (poly_base_coeffs_t)
          call bcast(bc, root_rank)
       type is (hom_base_coeffs_t)
          call bcast(bc, root_rank)
       class default
          $ABORT(Unsupported type)
       end select
          
    endif

    ! Finish

  end subroutine bcast_alloc_bc

  $endif

end module gyre_base_coeffs_mpi
