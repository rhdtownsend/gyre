! Module   : gyre_seg
! Purpose  : grid segment
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

module gyre_seg

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: seg_t
     real(WP), allocatable :: x(:)
     integer               :: i_a
     integer               :: i_b
     integer               :: n
  end type seg_t

  ! Interfaces

  interface seg_t
     module procedure seg_t_part_
  end interface seg_t

  ! Access specifiers

  private

  public :: seg_t

contains

  function seg_t_part_ (x) result (sg)

    DONT YET KNOW HOW TO DO THIS

    ! real(WP), intent(in)     :: x(:)
    ! type(seg_t), allocatable :: sg(:)

    ! logical :: mask(SIZE(x)-1)
    ! integer :: n

    ! ! Construct an array of segments, by partitioning x at double points

    ! ! Set up a mask for the double points

    ! n = SIZE(x)

    ! mask = x(2:) == x(:n-1)

    ! ! Allocate segments

    ! n_s = COUNT(mask) + 1
    
    ! allocate(sg(n_s))

    ! ! Set up index ranges giving the position of the segment in x

    ! sg%i_a = [1,PACK([(i,i=2,n)], MASK=mask)]
    ! sg%i_b = [PACK([(i,i=1,n-1)], MASK=mask),n]

    ! ! Set up segment x values

    ! seg_loop : do s = 1, n_s
    !    sg%x = x(sg%i_a(s):sg%i_b(s))
    ! end do seg_loop

    ! ! Finish

    ! return

  end function seg_t_part_

end module gyre_seg
