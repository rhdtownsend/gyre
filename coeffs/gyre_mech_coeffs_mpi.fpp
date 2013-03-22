! Module   : gyre_mech_coeffs_mpi
! Purpose  : MPI support for gyre_mech_coeffs
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

module gyre_mech_coeffs_mpi

  ! Uses

  use core_parallel

  use gyre_mech_coeffs
  use gyre_mech_coeffs_evol
  use gyre_mech_coeffs_poly
  use gyre_mech_coeffs_hom

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  $if($MPI)
  interface bcast_alloc
     module procedure bcast_alloc_mc
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

  subroutine bcast_alloc_mc (mc, root_rank)

    class(mech_coeffs_t), allocatable, intent(inout) :: mc
    integer, intent(in)                              :: root_rank

    integer, parameter :: EVOL_TYPE = 1
    integer, parameter :: POLY_TYPE = 2
    integer, parameter :: HOM_TYPE = 3

    logical :: alloc
    integer :: type

    ! Deallocate the mech_coeffs on non-root processors

    if(MPI_RANK /= root_rank .AND. ALLOCATED(mc)) then
       deallocate(mc)
    endif

    ! Check if the mech_coeffs is allocated on the root processor

    if(MPI_RANK == root_rank) alloc = ALLOCATED(mc)
    call bcast(alloc, root_rank)

    if(alloc) then

       ! Broadcast the dynamic type

       if(MPI_RANK == root_rank) then

          select type (mc)
          type is (mech_coeffs_evol_t)
             type = EVOL_TYPE
          type is (mech_coeffs_poly_t)
             type = POLY_TYPE
          type is (mech_coeffs_hom_t)
             type = HOM_TYPE
          class default
             $ABORT(Unsupported type)
          end select
          
       end if

       call bcast(type, root_rank)

       ! Allocate the mech_coeffs

       if(MPI_RANK /= root_rank) then
          select case (type)
          case (EVOL_TYPE)
             allocate(mech_coeffs_evol_t::mc)
          case (POLY_TYPE)
             allocate(mech_coeffs_poly_t::mc)
          case(HOM_TYPE)
             allocate(mech_coeffs_hom_t::mc)
          case default
             $ABORT(Unsupported type)
          end select
       endif

       ! Broadcast the mech_coeffs

       call mc%bcast(root_rank)

    endif

    ! Finish

  end subroutine bcast_alloc_mc

  $endif

end module gyre_mech_coeffs_mpi
