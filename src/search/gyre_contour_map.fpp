! Module   : gyre_contour_map
! Purpose  : contour mapping
!
! Copyright 2015-2016 Rich Townsend
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

module gyre_contour_map

  ! Uses

  use core_kinds

  use gyre_contour_seg
  use gyre_ext

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: contour_map_t
     private
     type(r_ext_t), allocatable :: f(:,:)
     type(r_ext_t), allocatable :: x(:)
     type(r_ext_t), allocatable :: y(:)
     integer, public            :: n_x
     integer, public            :: n_y
   contains
     private
     procedure, public :: get_segs => get_segs_
  end type contour_map_t

  ! Interfaces

  interface contour_map_t
     module procedure contour_map_t_
  end interface contour_map_t
  
  ! Access specifiers

  private

  public :: contour_map_t

  ! Procedures

contains

  function contour_map_t_ (x, y, f) result (cm)

    type(r_ext_t), intent(in) :: x(:)
    type(r_ext_t), intent(in) :: y(:)
    type(r_ext_t), intent(in) :: f(:,:)
    type(contour_map_t)       :: cm

    $CHECK_BOUNDS(SIZE(f, 1),SIZE(x))
    $CHECK_BOUNDS(SIZE(f, 2),SIZE(y))

    ! Construct the contour_map_t

    cm%f = f

    cm%x = x
    cm%y = y

    cm%n_x = SIZE(x)
    cm%n_y = SIZE(y)
    
    ! Finish

    return

  end function contour_map_t_

  !****

  subroutine get_segs_ (this, i_x, i_y, cs)

    class(contour_map_t), intent(in)              :: this
    integer, intent(in)                           :: i_x
    integer, intent(in)                           :: i_y
    type(contour_seg_t), allocatable, intent(out) :: cs(:)

    integer :: code
    integer :: n_cs

    $ASSERT_DEBUG(i_x >= 1,Invalid index)
    $ASSERT_DEBUG(i_y >= 1,Invalid index)
    
    $ASSERT_DEBUG(i_x < this%n_x,Invalid index)
    $ASSERT_DEBUG(i_y < this%n_y,Invalid index)
    
    ! Evaluate the cell code

    code = 0

    if (this%f(i_x  ,i_y  ) > 0._WP) code = code + 1
    if (this%f(i_x+1,i_y  ) > 0._WP) code = code + 2
    if (this%f(i_x+1,i_y+1) > 0._WP) code = code + 4
    if (this%f(i_x  ,i_y+1) > 0._WP) code = code + 8

    ! Add segments based on the code

    n_cs = 0

    allocate(cs(2))

    select case (code)
    case (0)
    case (1)
       call add_segment_(1, 4)
    case (2)
       call add_segment_(2, 1)
    case (3)
       call add_segment_(2, 4)
    case (4)
       call add_segment_(3, 2)
    case (5)
       call add_segment_(1, 2)
       call add_segment_(3, 4)
    case (6)
       call add_segment_(3, 1)
    case (7)
       call add_segment_(3, 4)
    case (8)
       call add_segment_(4, 3)
    case (9)
       call add_segment_(1, 3)
    case (10)
       call add_segment_(4, 1)
       call add_segment_(2, 3)
    case (11)
       call add_segment_(2, 3)
    case (12)
       call add_segment_(4, 2)
    case (13)
       call add_segment_(1, 2)
    case (14)
       call add_segment_(4, 1)
    case (15)
    case default
       $ABORT(Invalid code)
    end select

    cs = cs(:n_cs)

    ! Finish

  contains

    subroutine add_segment_ (e_a, e_b)

      integer, intent(in) :: e_a
      integer, intent(in) :: e_b

      type(r_ext_t) :: x_a
      type(r_ext_t) :: x_b
      type(r_ext_t) :: y_a
      type(r_ext_t) :: y_b

      ! Add a segment between edges e_a and e_b

      call interp_isect_(e_a, x_a, y_a)
      call interp_isect_(e_b, x_b, y_b)

      n_cs = n_cs + 1

      cs(n_cs) = contour_seg_t([x_a,x_b], [y_a,y_b])

      ! Finish

      return

    end subroutine add_segment_

    subroutine interp_isect_(e, x, y)

      integer, intent(in)        :: e
      type(r_ext_t), intent(out) :: x
      type(r_ext_t), intent(out) :: y

      type(r_ext_t) :: x_a
      type(r_ext_t) :: x_b
      type(r_ext_t) :: y_a
      type(r_ext_t) :: y_b
      type(r_ext_t) :: f_a
      type(r_ext_t) :: f_b
      type(r_ext_t) :: w

      ! Interpolate the contour intersection with edge e

      select case (e)

      case (1)

         x_a = this%x(i_x  )
         x_b = this%x(i_x+1)

         y_a = this%y(i_y  )
         y_b = this%y(i_y  )

         f_a = this%f(i_x  ,i_y  )
         f_b = this%f(i_x+1,i_y  )

      case (2)

         x_a = this%x(i_x+1)
         x_b = this%x(i_x+1)

         y_a = this%y(i_y  )
         y_b = this%y(i_y+1)

         f_a = this%f(i_x+1,i_y  )
         f_b = this%f(i_x+1,i_y+1)

      case (3)

         x_a = this%x(i_x+1)
         x_b = this%x(i_x  )

         y_a = this%y(i_y+1)
         y_b = this%y(i_y+1)

         f_a = this%f(i_x+1,i_y+1)
         f_b = this%f(i_x  ,i_y+1)

      case (4)

         x_a = this%x(i_x  )
         x_b = this%x(i_x  )

         y_a = this%y(i_y+1)
         y_b = this%y(i_y  )

         f_a = this%f(i_x  ,i_y+1)
         f_b = this%f(i_x  ,i_y  )

      case default

         $ABORT(Invalid edge)

      end select

      w = -f_a/(f_b - f_a)

      x = x_a*(1._WP-w) + w*x_b
      y = y_a*(1._WP-w) + w*y_b

      ! Finish

      return

    end subroutine interp_isect_

  end subroutine get_segs_

end module gyre_contour_map
