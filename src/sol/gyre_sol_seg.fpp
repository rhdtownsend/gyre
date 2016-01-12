! Module   : gyre_sol_seg
! Purpose  : solution data segment
!
! Copyright 2016 Rich Townsend
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

module gyre_sol_seg

  ! Uses

  use core_kinds

  use gyre_spline

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: sol_seg_t
     private
     type(c_spline_t) :: sp_y(6)
     logical          :: df_y(6) = .FALSE.
   contains
     private
     procedure, public :: set_y
     procedure, public :: y
  end type sol_seg_t

  ! Interfaces

  interface sol_seg_t
     module procedure sol_seg_t_
  end interface sol_seg_t

  ! Access specifiers

  private

  public :: sol_seg_t

  ! Procedures

contains

  function sol_seg_t_ () result (ss)

    type(sol_seg_t) :: ss

    ! Construct the sol_seg_t

    ! Finish

    return

  end function sol_seg_t_

  !****

  subroutine set_y (this, i, x, y, dy_dx)

    class(sol_seg_t), intent(inout) :: this
    integer, intent(in)             :: i
    real(WP), intent(in)            :: x(:)
    complex(WP), intent(in)         :: y(:)
    complex(WP), intent(in)         :: dy_dx(:)

    $ASSERT_DEBUG(i >= 1,Invalid index)
    $ASSERT_DEBUG(i <= 6,Invalid index)

    $CHECK_BOUNDS(SIZE(y),SIZE(x))
    $CHECK_BOUNDS(SIZE(dy_dx),SIZE(x))

    ! Set the data for y(i)

    this%sp_y(i) = c_spline_t(x, y, dy_dx)
    
    this%df_y(i) = .TRUE.

    ! Finish

    return

  end subroutine set_y

  !****

  function y (this, i, x)

    class(sol_seg_t), intent(in) :: this
    integer, intent(in)          :: i
    real(WP), intent(in)         :: x
    complex(WP)                  :: y

    ! Interpolate y(i)

    if (this%df_y(i)) then

       y = this%sp_y(i)%f(x)
       
    else

       write(ERROR_UNIT, *) i

       $ABORT(No solution data provided for y)

    endif
       
    ! Finish

    return

  end function y

end module gyre_sol_seg
