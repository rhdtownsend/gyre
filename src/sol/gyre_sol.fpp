! Module   : gyre_sol
! Purpose  : solution data in canonical form
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

module gyre_sol

  ! Uses

  use core_kinds

  use gyre_ext
  use gyre_sol_seg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: sol_t
     private
     type(sol_seg_t), allocatable  :: ss(:)
     integer, allocatable, public  :: s(:)
     real(WP), allocatable, public :: x(:)
     complex(WP), public           :: omega
     type(c_ext_t), public         :: discrim
     integer, public               :: n_s
     integer, public               :: n_k
   contains
     private
     procedure, public :: set_y
     procedure         :: y_1_
     procedure         :: y_v_
     generic, public   :: y => y_1_, y_v_
  end type sol_t

  ! Interfaces

  interface sol_t
     module procedure sol_t_
  end interface sol_t

  ! Access specifiers

  private

  public :: sol_t

  ! Procedures

contains

  function sol_t_ (s, x, omega, discrim) result (sl)

    integer, intent(in)       :: s(:)
    real(WP), intent(in)      :: x(:)
    complex(WP), intent(in)   :: omega
    type(c_ext_t), intent(in) :: discrim
    type(sol_t)               :: sl

    $CHECK_BOUNDS(SIZE(x),SIZE(s))

    ! Construct the sol_t

    sl%s = s
    sl%x = x

    sl%n_k = SIZE(sl%x)
    sl%n_s = MAXVAL(sl%s)

    allocate(sl%ss(sl%n_s), SOURCE=sol_seg_t())

    sl%omega = omega
    sl%discrim = discrim

    ! Finish

    return

  end function sol_t_

  !****

  subroutine set_y (this, i, y, dy_dx)

    class(sol_t), intent(inout) :: this
    integer, intent(in)         :: i
    complex(WP), intent(in)     :: y(:)
    complex(WP), intent(in)     :: dy_dx(:)

    integer :: s
    logical :: mask(this%n_k)

    $CHECK_BOUNDS(SIZE(y),this%n_k)
    $CHECK_BOUNDS(SIZE(dy_dx),this%n_k)

    ! Set the data for y(i)

    seg_loop : do s = 1, this%n_s

       mask = this%s == s

       call this%ss(s)%set_y(i, PACK(this%x, mask), PACK(y, mask), PACK(dy_dx, mask))

     end do seg_loop

    ! Finish

    return

  end subroutine set_y

  !****

  function y_1_ (this, i, s, x) result (y)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: i
    integer, intent(in)      :: s
    real(WP), intent(in)     :: x
    complex(WP)              :: y

    $ASSERT_DEBUG(s >= 1,Invalid segment index)
    $ASSERT_DEBUG(s <= this%n_s,Invalid segment index)

    ! Evaluate y(i)

    y = this%ss(s)%y(i, x)

    ! Finish

    return

  end function y_1_
  
  !****

  function y_v_ (this, i, s, x) result (y)

    class(sol_t), intent(in) :: this
    integer, intent(in)      :: i
    integer, intent(in)      :: s(:)
    real(WP), intent(in)     :: x(:)
    complex(WP)              :: y(SIZE(s))

    integer :: k

    $CHECK_BOUNDS(SIZE(x),SIZE(s))

    ! Evaluate y(i)

    !$OMP PARALLEL DO
    do k = 1, SIZE(s)
       y(k) = this%y(i, s(k), x(k))
    end do

    ! Finish

    return

  end function y_v_

end module gyre_sol
