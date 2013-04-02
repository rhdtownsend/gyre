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

module gyre_gridpar

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: gridpar_t
     private
     real(WP), public          :: alpha_osc
     real(WP), public          :: alpha_exp
     integer, public           :: n_center
     integer, public           :: n_floor
     real(WP), public          :: s
     integer, public           :: n_grid 
     character(LEN=64), public :: grid_type
   contains
     private
     procedure, public :: init
  end type gridpar_t

  ! Interfaces

  $if($MPI)

  interface bcast
     module procedure bcast_gp
  end interface bcast

  $endif

  ! Access specifiers

  private

  public :: gridpar_t
  $if($MPI)
  public :: bcast
  $endif

  ! Procedures

contains

  subroutine init (this, alpha_osc, alpha_exp, n_center, n_floor, s, n_grid, grid_type)

    class(gridpar_t), intent(out) :: this
    real(WP), intent(in)          :: alpha_osc
    real(WP), intent(in)          :: alpha_exp
    integer, intent(in)           :: n_center
    integer, intent(in)           :: n_floor
    real(WP), intent(in)          :: s
    integer, intent(in)           :: n_grid
    character(LEN=*), intent(in)  :: grid_type

    ! Initialize the gridpar

    this%alpha_osc = alpha_osc
    this%alpha_exp = alpha_exp

    this%n_center = n_center
    this%n_floor = n_floor

    this%s = s
    this%n_grid = n_grid

    this%grid_type = grid_type

    ! Finish

    return

  end subroutine init

!****

  $if($MPI)

  subroutine bcast_gp (gp, root_rank)

    type(gridpar_t), intent(inout) :: gp
    integer, intent(in)            :: root_rank

    ! Broadcast the gridpar

    call bcast(gp%alpha_osc, root_rank)
    call bcast(gp%alpha_exp, root_rank)

    call bcast(gp%n_center, root_rank)
    call bcast(gp%n_floor, root_rank)

    call bcast(gp%s, root_rank)
    call bcast(gp%n_grid, root_rank)

    call bcast(gp%grid_type, root_rank)

    ! Finish

    return

  end subroutine bcast_gp

  $endif

end module gyre_gridpar
