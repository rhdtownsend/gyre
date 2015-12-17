! Module   : gyre_model_coords
! Purpose  : coordinates within a stellar model 
!
! Copyright 2015 Rich Townsend
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

module gyre_model_coords

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: model_coords_t
     real(WP) :: x
     integer  :: s
  end type model_coords_t

  ! Access specifiers

  private

  public :: model_coord_t
  public :: partition

  ! Procedures

contains

  subroutine partition (x, mc, n_s)

    real(WP), intent(in)             :: x(:)
    type(model_coord_t), intent(out) :: mc(:)
    integer, intent(out)             :: n_s

    integer :: s
    integer :: i

    $CHECK_BOUNDS(SIZE(mp),SIZE(x))

    ! Partition the array x into strictly-monotonic-increasing
    ! segments, by splitting at double points; return model
    ! coordinates in mc, and the number of segments in n_s

    s = 1

    mc(1) = model_coords_t(x(1), s)

    x_loop : do i = 2, SIZE(x)
       
       if (x(i) == x(i-1)) s = s + 1

       mc(i) = model_coords_t(x(i), s)

    end do x_loop

    n_s = s

    ! Finish

    return

  end subroutine partition
    
end module gyre_model_coords
