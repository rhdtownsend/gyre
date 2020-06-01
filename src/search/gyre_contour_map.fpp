! Module   : gyre_contour_map
! Purpose  : contour mapping
!
! Copyright 2015-2020 Rich Townsend & The GYRE team
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
  use core_memory
  $if ($HDF5)
  use core_hgroup
  $endif

  use gyre_contour_path
  use gyre_ext

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: contour_map_t
     private
     type(c_ext_t), allocatable :: f(:,:)
     real(WP), allocatable      :: z_re(:)
     real(WP), allocatable      :: z_im(:)
     integer, public            :: n_re
     integer, public            :: n_im
   contains
     private
     procedure, public :: trace_paths
  end type contour_map_t

  ! Interfaces

  interface contour_map_t
     module procedure contour_map_t_
  end interface contour_map_t
  
  $if ($HDF5)
  interface write
     module procedure write_
  end interface write
  $endif

  ! Access specifiers

  private

  public :: contour_map_t
  $if ($HDF5)
  public :: write
  $endif

  ! Procedures

contains

  function contour_map_t_ (z_re, z_im, f) result (cm)

    real(WP), intent(in) :: z_re(:)
    real(WP), intent(in) :: z_im(:)
    type(c_ext_t)        :: f(:,:)
    type(contour_map_t)  :: cm

    $CHECK_BOUNDS(SIZE(f, 1),SIZE(z_re))
    $CHECK_BOUNDS(SIZE(f, 2),SIZE(z_im))

    ! Construct the contour_map_t

    cm%f = f

    cm%z_re = z_re
    cm%z_im = z_im

    cm%n_re = SIZE(z_re)
    cm%n_im = SIZE(z_im)
    
    ! Finish

    return

  end function contour_map_t_

  !****

  subroutine trace_paths (this, cp_re, cp_im, process_isect)

    class(contour_map_t), intent(in)               :: this
    type(contour_path_t), allocatable, intent(out) :: cp_re(:)
    type(contour_path_t), allocatable, intent(out) :: cp_im(:)
    interface
       subroutine process_isect (z_a_re, z_b_re, z_a_im, z_b_im)
         use core_kinds
         complex(WP), intent(in) :: z_a_re
         complex(WP), intent(in) :: z_b_re
         complex(WP), intent(in) :: z_a_im
         complex(WP), intent(in) :: z_b_im
       end subroutine process_isect
    end interface

    integer, parameter :: D = 128

    integer                  :: n_cp_re
    integer                  :: n_cp_im
    integer                  :: i_re
    integer                  :: i_im
    integer                  :: code_re
    integer                  :: code_im
    complex(WP), allocatable :: z_a_re(:)
    complex(WP), allocatable :: z_b_re(:)
    complex(WP), allocatable :: z_a_im(:)
    complex(WP), allocatable :: z_b_im(:)
    
    ! Loop through the cells, adding segments to the path lists

    allocate(cp_re(D))
    allocate(cp_im(D))

    n_cp_re = 0
    n_cp_im = 0

    cell_x_loop : do i_re = 1, this%n_re-1
       cell_y_loop : do i_im = 1, this%n_im-1

          ! Evaluate cell codes

          code_re = cell_code_(real_part(this%f(i_re:i_re+1,i_im:i_im+1)))
          code_im = cell_code_(imag_part(this%f(i_re:i_re+1,i_im:i_im+1)))

          ! Add segments based on the codes

          call add_segments_(code_re, cp_re, n_cp_re, z_a_re, z_b_re, 're')
          call add_segments_(code_im, cp_im, n_cp_im, z_a_im, z_b_im, 'im')

          ! Check whether an intersection is at all possible

          if (code_re /= 0 .AND. code_re /= 15 .AND. &
              code_im /= 0 .AND. code_im /= 15) then

             ! Perform detailed intersection check

             call check_isect(z_a_re, z_b_re, z_a_im, z_b_im, process_isect)

          end if

       end do cell_y_loop
    end do cell_x_loop
          
    ! Shrink the path lists

    call reallocate(cp_re, [n_cp_re])
    call reallocate(cp_im, [n_cp_im])

    ! Finish

    return

  contains

    subroutine add_segments_ (code, cp, n_cp, z_a, z_b, part)

      integer, intent(in)                              :: code
      type(contour_path_t), allocatable, intent(inout) :: cp(:)
      integer, intent(inout)                           :: n_cp
      complex(WP), allocatable, intent(out)            :: z_a(:)
      complex(WP), allocatable, intent(out)            :: z_b(:)
      character(*), intent(in)                         :: part

      complex(WP)   :: z_a_1
      complex(WP)   :: z_a_2
      complex(WP)   :: z_b_1
      complex(WP)   :: z_b_2

      ! Add segments based on the code

      select case (code)

      case (0)

         allocate(z_a(0))
         allocate(z_b(0))

      case (1)

         call add_segment_(1, 4, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (2)

         call add_segment_(2, 1, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (3)

         call add_segment_(2, 4, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (4)

         call add_segment_(3, 2, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (5)

         call add_segment_(1, 2, cp, n_cp, z_a_1, z_b_1, part)
         call add_segment_(3, 4, cp, n_cp, z_a_2, z_b_2, part)

         z_a = [z_a_1,z_a_2]
         z_b = [z_b_1,z_b_2]

      case (6)

         call add_segment_(3, 1, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (7)

         call add_segment_(3, 4, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (8)

         call add_segment_(4, 3, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (9)

         call add_segment_(1, 3, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (10)

         call add_segment_(4, 1, cp, n_cp, z_a_1, z_b_1, part)
         call add_segment_(2, 3, cp, n_cp, z_a_2, z_b_2, part)

         z_a = [z_a_1,z_a_2]
         z_b = [z_b_1,z_b_2]

      case (11)

         call add_segment_(2, 3, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (12)

         call add_segment_(4, 2, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (13)

         call add_segment_(1, 2, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

     case (14)

         call add_segment_(4, 1, cp, n_cp, z_a_1, z_b_1, part)

         z_a = [z_a_1]
         z_b = [z_b_1]

      case (15)

         allocate(z_a(0))
         allocate(z_b(0))

      case default

         $ABORT(Invalid code)

      end select

      ! Finish

      return

    end subroutine add_segments_

    !****

    subroutine add_segment_ (e_a, e_b, cp, n_cp, z_a, z_b, part)

      integer, intent(in)                              :: e_a
      integer, intent(in)                              :: e_b
      type(contour_path_t), allocatable, intent(inout) :: cp(:)
      integer, intent(inout)                           :: n_cp
      complex(WP), intent(out)                         :: z_a
      complex(WP), intent(out)                         :: z_b
      character(*), intent(in)                         :: part

      integer              :: tag_a
      integer              :: tag_b
      integer              :: i_cp_a
      integer              :: i_cp_b
      logical              :: tail_a
      logical              :: tail_b
      integer              :: i_cp
      type(contour_path_t) :: cp_new

      ! Add a segment between edges with positions e_a and e_b

      ! Map the edges indices into unique tags

      tag_a = edge_tag_(i_re, i_im, this%n_re, this%n_im, e_a)
      tag_b = edge_tag_(i_re, i_im, this%n_re, this%n_im, e_b)

      ! See if the tags connect to existing paths

      i_cp_a = 0
      i_cp_b = 0

      tail_a = .FALSE.
      tail_b = .FALSE.

      do i_cp = 1, n_cp

         if (cp(i_cp)%check_tag(tag_a, tail=.FALSE.)) then
            i_cp_a = i_cp
         elseif (cp(i_cp)%check_tag(tag_a, tail=.TRUE.)) then
            i_cp_a = i_cp
            tail_a = .TRUE.
         endif

         if (cp(i_cp)%check_tag(tag_b, tail=.FALSE.)) then
            i_cp_b = i_cp
         elseif (cp(i_cp)%check_tag(tag_b, tail=.TRUE.)) then
            i_cp_b = i_cp
            tail_b = .TRUE.
         endif

      end do

      ! Add the segment

      if (i_cp_a /= 0 .AND. i_cp_b /= 0) then

         ! Join paths i_cp_a and i_cp_b

         z_a = cp(i_cp_a)%get_point(tail_a)
         z_b = cp(i_cp_b)%get_point(tail_b)

         cp_new = contour_path_t(cp(i_cp_a), cp(i_cp_b), tail_a=tail_a, head_b=.NOT. tail_b)

         call delete_paths_(cp, n_cp, [i_cp_a,i_cp_b])
         
         call add_path_(cp, n_cp, cp_new)

      elseif (i_cp_a /= 0) then

         ! Add to path i_cp_a

         z_a = cp(i_cp_a)%get_point(tail_a)
         z_b = interp_point_(e_b, part)
 
         call cp(i_cp_a)%add_point(z_b, tag_b, tail_a)

      elseif (i_cp_b /= 0) then

         ! Add to path i_cp_b

         z_a = interp_point_(e_a, part)
         z_b = cp(i_cp_b)%get_point(tail_b)

         call cp(i_cp_b)%add_point(z_a, tag_a, tail_b)

      else

         ! Create a new path

         z_a = interp_point_(e_a, part)
         cp_new = contour_path_t(z_a, tag_a)

         z_b = interp_point_(e_b, part)
         call cp_new%add_point(z_b, tag_b)

         call add_path_(cp, n_cp, cp_new)

      endif
            
      ! Finish

      return

    end subroutine add_segment_

    !****

    function interp_point_ (e, part) result (z)

      integer, intent(in)      :: e
      character(*), intent(in) :: part
      complex(WP)              :: z

      complex(WP)   :: z_a
      complex(WP)   :: z_b
      type(c_ext_t) :: f_a
      type(c_ext_t) :: f_b
      real(WP)      :: w

      ! Interpolate the contour intersection with edge e

      select case (e)

      case (1)

         z_a = CMPLX(this%z_re(i_re  ), this%z_im(i_im  ), WP)
         z_b = CMPLX(this%z_re(i_re+1), this%z_im(i_im  ), WP)

         f_a = this%f(i_re  ,i_im  )
         f_b = this%f(i_re+1,i_im  )

      case (2)

         z_a = CMPLX(this%z_re(i_re+1), this%z_im(i_im  ), WP)
         z_b = CMPLX(this%z_re(i_re+1), this%z_im(i_im+1), WP)

         f_a = this%f(i_re+1,i_im  )
         f_b = this%f(i_re+1,i_im+1)

      case (3)

         z_a = CMPLX(this%z_re(i_re+1), this%z_im(i_im+1), WP)
         z_b = CMPLX(this%z_re(i_re  ), this%z_im(i_im+1), WP)

         f_a = this%f(i_re+1,i_im+1)
         f_b = this%f(i_re  ,i_im+1)

      case (4)

         z_a = CMPLX(this%z_re(i_re  ), this%z_im(i_im+1), WP)
         z_b = CMPLX(this%z_re(i_re  ), this%z_im(i_im  ), WP)

         f_a = this%f(i_re  ,i_im+1)
         f_b = this%f(i_re  ,i_im  )

      case default

         $ABORT(Invalid edge)

      end select

      select case (part)
      case ('re')
         w = real(-real_part(f_a)/(real_part(f_b) - real_part(f_a)))
      case ('im')
         w = real(-imag_part(f_a)/(imag_part(f_b) - imag_part(f_a)))
      case default
         $ABORT(Invalid part)
      end select

      $ASSERT_DEBUG(w >= 0._WP)
      $ASSERT_DEBUG(w <= 1._WP)

      z = (1._WP-w)*z_a + w*z_b

      ! Finish

      return

    end function interp_point_

  end subroutine trace_paths

  !****

  function cell_code_ (f) result (code)

    class(r_ext_t), intent(in) :: f(:,:)
    integer                    :: code

    $CHECK_BOUNDS(SIZE(f, 1),2)
    $CHECK_BOUNDS(SIZE(f, 2),2)

    ! Evaluate the cell code

    code = 0

    if (f(1,1) > 0._WP) code = code + 1
    if (f(2,1) > 0._WP) code = code + 2
    if (f(2,2) > 0._WP) code = code + 4
    if (f(1,2) > 0._WP) code = code + 8

    ! Finish

    return

  end function cell_code_

  !****

  function edge_tag_ (i_re, i_im, n_re, n_im, e) result (tag)

    integer, intent(in) :: i_re
    integer, intent(in) :: i_im
    integer, intent(in) :: n_re
    integer, intent(in) :: n_im
    integer, intent(in) :: e
    integer             :: tag

    ! Map an edge index into a tag

    select case (e)

    case (1)

       tag = i_re + (i_im-1)*n_re

    case (2)

       tag = n_re*n_im + i_re+1 + (i_im-1)*n_re

    case (3)

       tag = i_re + i_im*n_re

    case (4)

       tag = n_re*n_im + i_re + (i_im-1)*n_re

    case default

       $ABORT(Invalid edge)

    end select

    ! Finish

    return

  end function edge_tag_

  !****
  
  subroutine add_path_ (cp, n_cp, cp_new)

    type(contour_path_t), allocatable, intent(inout) :: cp(:)
    integer, intent(inout)                           :: n_cp
    type(contour_path_t), intent(in)                 :: cp_new

    integer :: d_cp

    ! Add a new path to the list

    d_cp = SIZE(cp)

    if (n_cp >= d_cp) then
       call reallocate(cp, [2*d_cp])
    endif
    
    n_cp = n_cp + 1

    cp(n_cp) = cp_new

    ! Finish

    return

  end subroutine add_path_

  !****

  subroutine delete_paths_ (cp, n_cp, i_cp)

    type(contour_path_t), intent(inout) :: cp(:)
    integer, intent(inout)              :: n_cp
    integer, intent(in)                 :: i_cp(:)

    integer :: j
    integer :: k

    ! Delete paths from the list

    j = 0

    do k = 1, n_cp
       if (.NOT. ANY(i_cp == k)) then
          j = j + 1
          cp(j) = cp(k)
       endif
    end do

    n_cp = j

    ! Finish

    return

  end subroutine delete_paths_

  !****

  subroutine check_isect (z_a_re, z_b_re, z_a_im, z_b_im, process_isect)

    complex(WP), intent(in)   :: z_a_re(:)
    complex(WP), intent(in)   :: z_b_re(:)
    complex(WP), intent(in)   :: z_a_im(:)
    complex(WP), intent(in)   :: z_b_im(:)
    interface
       subroutine process_isect (z_a_re, z_b_re, z_a_im, z_b_im)
         use core_kinds
         complex(WP), intent(in) :: z_a_re
         complex(WP), intent(in) :: z_b_re
         complex(WP), intent(in) :: z_a_im
         complex(WP), intent(in) :: z_b_im
       end subroutine process_isect
    end interface
    
    integer     :: i_re
    integer     :: i_im
    real(WP)    :: M(2,2)
    real(WP)    :: R(2)
    real(WP)    :: D
    real(WP)    :: w_re
    real(WP)    :: w_im

    $CHECK_BOUNDS(SIZE(z_a_re), SIZE(z_b_re))
    $CHECK_BOUNDS(SIZE(z_a_im), SIZE(z_b_im))

    ! Check for an intersection between the line segments
    ! [z_a_re,z_b_re] and [z_a_im,z_b_im]
      
    do i_re = 1, SIZE(z_a_re)
       do i_im = 1, SIZE(z_a_im)

          ! Solve for the intersection point

          M(1,1) = REAL(z_b_re(i_re)) - REAL(z_a_re(i_re))
          M(1,2) = REAL(z_a_im(i_im)) - REAL(z_b_im(i_im))
          
          M(2,1) = AIMAG(z_b_re(i_re)) - AIMAG(z_a_re(i_re))
          M(2,2) = AIMAG(z_a_im(i_im)) - AIMAG(z_b_im(i_im))

          R(1) = -REAL(z_a_re(i_re)) + REAL(z_a_im(i_im))
          R(2) = -AIMAG(z_a_re(i_re)) + AIMAG(z_a_im(i_im))

          D = M(1,1)*M(2,2) - M(1,2)*M(2,1)

          if (D /= 0._WP) then

             ! Lines intersect

             w_re = (M(2,2)*R(1) - M(1,2)*R(2))/D
             w_im = (M(1,1)*R(2) - M(2,1)*R(1))/D

             if (w_re >= 0._WP .AND. w_re <= 1._WP .AND. &
                 w_im >= 0._WP .AND. w_im <= 1._WP) then

                ! Lines intersect within their endpoints

                call process_isect(z_a_re(i_re), z_b_re(i_re), z_a_im(i_im), z_b_im(i_im))

             end if

          endif

       end do
    end do

    ! Finish

    return

  end subroutine check_isect

  !****

  $if ($HDF5)

  subroutine write_ (hg, cm)

    type(hgroup_t), intent(inout)   :: hg
    type(contour_map_t), intent(in) :: cm

    ! Write the contour_map_t

    call write_dset(hg, 'z_re', cm%z_re)
    call write_dset(hg, 'z_im', cm%z_im)

    call write_dset(hg, 'f_f', fraction(cm%f))
    call write_dset(hg, 'f_e', exponent(cm%f))

    ! Finish

    return

  end subroutine write_

  $endif

end module gyre_contour_map
