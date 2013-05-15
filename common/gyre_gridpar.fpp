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
$include 'core_parallel.inc'

module gyre_gridpar

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: gridpar_t
     real(WP)          :: x_i = 0._WP
     real(WP)          :: x_o = 1._WP
     real(WP)          :: alpha_osc = 0._WP
     real(WP)          :: alpha_exp = 0._WP
     real(WP)          :: omega_a = 0._WP
     real(WP)          :: omega_b = 0._WP
     real(WP)          :: s = 0
     integer           :: n = 0
     character(LEN=64) :: op_type = ''
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

    call bcast(gp%x_i, root_rank)
    call bcast(gp%x_o, root_rank)

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

  $BCAST_ALLOC(gp,type(gridpar_t),0)
  $BCAST_ALLOC(gp,type(gridpar_t),1)

  $endif

end module gyre_gridpar
