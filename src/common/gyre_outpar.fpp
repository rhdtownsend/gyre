! Module   : gyre_outpar
! Purpose  : output parameters
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

module gyre_outpar

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_constants

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: outpar_t
     character(256)          :: freq_units
     character(FILENAME_LEN) :: summary_file
     character(256)          :: summary_file_format
     character(2048)         :: summary_item_list
     character(FILENAME_LEN) :: mode_prefix
     character(FILENAME_LEN) :: mode_template
     character(256)          :: mode_file_format
     character(2048)         :: mode_item_list
     logical                 :: prune_modes
  end type outpar_t

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

  public :: outpar_t
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  $if ($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_${RANK}_ (up, root_rank)

    type(outpar_t), intent(inout) :: up$ARRAY_SPEC($RANK)
    integer, intent(in)           :: root_rank

    ! Broadcast the outpar_t

    call bcast(up%freq_units, root_rank)
    call bcast(up%summary_file, root_rank)
    call bcast(up%summary_file_format, root_rank)
    call bcast(up%summary_item_list, root_rank)
    call bcast(up%mode_prefix, root_rank)
    call bcast(up%mode_template, root_rank)
    call bcast(up%mode_file_format, root_rank)
    call bcast(up%mode_item_list, root_rank)
    
    call bcast(up%prune_modes, root_rank)

    ! Finish

    return

  end subroutine bcast_${RANK}_

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $BCAST_ALLOC(type(outpar_t),0)
  $BCAST_ALLOC(type(outpar_t),1)

  $endif

end module gyre_outpar
