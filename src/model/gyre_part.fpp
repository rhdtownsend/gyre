! Module   : gyre_part
! Purpose  : array partitioning
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

module gyre_part

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: part_t
     private
     integer, allocatable :: i_a(:)
     integer, allocatable :: i_b(:)
     integer              :: n
   contains
     private
     procedure, public :: i_dbl => i_dbl_
     procedure, public :: i_unq => i_unq_
     procedure, public :: i_seg => i_seg_
     procedure, public :: n_seg => n_seg_
  end type part_t
 
  ! Interfaces

  interface part_t
     module procedure part_t_
  end interface part_t

  ! Access specifiers

  private

  public :: part_t

  ! Procedures

contains

  function part_t_ (x) result (pt)

    real(WP), intent(in) :: x(:)
    type(part_t)         :: pt

    integer :: n
    logical :: mask(SIZE(x)-1)
    integer :: s
    integer :: i

    ! Construct the part_t

    ! Create a mask for the array subintervals bounded by double
    ! points

    n = SIZE(x)

    mask = x(2:) == x(:n-1)

    ! Set up indices fr the segments delineated by double points

    pt%n = COUNT(mask) + 1

    allocate(pt%i_a(pt%n))
    allocate(nt%i_b(pt%n))

    s = 1

    pt%i_a(s) = 1

    mask_loop : do i = 1, n-1
       if (mask(i)) then
          pt%i_b(s) = i
          s = s + 1
          pt%i_a(s) = i + 1
       end if
    end do mask_loop
          
    pt%i_b(s) = n

    ! Finish

    return

  end function part_t_

  !****

  function i_dbl_ (this) result (i_dbl)

    class(part_t), intent(in) :: this
    integer, allocatable      :: i_dbl(:)

    ! Return indices for the left-hand elements of double-point pairs

    ! Not correct!

    i_dbl = this%i_b

    ! Finish

    return

  end function i_double_

  !****

  function i_unq_ (this) result (i_unq)

    class(part_t), intent(in) :: this
    integer, allocatable      :: i_unq(:)

    integer :: s

    ! Return indices for unique elements

    i_unq = [(i,i=this%i_a(1),this%i_b(1))]

    do s = 2, this%n
       i_unq = [i_unq,[(i,i=this%i_a(s)+1,this%i_b(s))]]
    end do

    ! Finish

    return

  end function i_unq_

  !****

  function i_seg_ (this, s) result (i_seg)

    class(part_t), intent(in) :: this
    integer, intent(in)       :: s
    integer, allocatable      :: i_seg(:)

    integer :: i

    ! Return indices for elements comprising the s'th
    ! strictly-monotonic segment

    i_seg = [(i,i=this%i_a(s),this%i_b(s))]

    ! Finish

    return

  end function i_seg_

  !****

  function n_seg_ (this) result (n_seg)

    class(part_t), intent(in) :: this
    integer                   :: n_seg

    ! Return the number of strictly-monotonic segments

    n_seg = this%n

    ! Finish

    return

  end function n_seg_

end module gyre_part
