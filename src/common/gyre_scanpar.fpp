! Module   : gyre_scanpar
! Purpose  : frequency scan parameters
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

module gyre_scanpar

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: scanpar_t
     real(WP)            :: freq_min
     real(WP)            :: freq_max
     integer             :: n_freq
     character(LEN=64)   :: grid_type
     character(LEN=64)   :: freq_units
     character(LEN=2048) :: tag_list
  end type scanpar_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_sp_0
     module procedure bcast_sp_1
  end interface bcast

  interface bcast_alloc
     module procedure bcast_alloc_sp_0
     module procedure bcast_alloc_sp_1
  end interface bcast_alloc

  $endif

  ! Access specifiers

  private

  public :: scanpar_t
  $if($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  $if($MPI)

  $define $BCAST $sub

  $local $RANK $1

  subroutine bcast_sp_$RANK (sp, root_rank)

    type(scanpar_t), intent(inout) :: sp$ARRAY_SPEC($RANK)
    integer, intent(in)            :: root_rank

    ! Broadcast the scanpar

    call bcast(sp%freq_min, root_rank)
    call bcast(sp%freq_max, root_rank)
    call bcast(sp%n_freq, root_rank)

    call bcast(sp%grid_type, root_rank)
    call bcast(sp%freq_units, root_rank)
    call bcast(sp%tag_list, root_rank)

    ! Finish

    return

  end subroutine bcast_sp_$RANK

  $endsub

  $BCAST(0)
  $BCAST(1)

!****

  $BCAST_ALLOC(sp,type(scanpar_t),0)
  $BCAST_ALLOC(sp,type(scanpar_t),1)

  $endif

end module gyre_scanpar
