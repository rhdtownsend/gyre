! Module   : gyre_soln
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
$include 'core_parallel.inc'

module gyre_soln

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_ext
  use gyre_grid
  use gyre_point
  use gyre_soln_seg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: soln_t
     private
     type(grid_t)                  :: gr
     type(soln_seg_t), allocatable :: ss(:)
     integer                       :: s_i
     integer                       :: s_o
     integer                       :: n_k
     complex(WP), public           :: omega
     type(c_ext_t), public         :: discrim
   contains
     private
     procedure, public :: set_y
     procedure         :: y_1_
     procedure         :: y_v_
     generic, public   :: y => y_1_, y_v_
     procedure, public :: grid
  end type soln_t

  ! Interfaces

  interface soln_t
     module procedure soln_t_
  end interface soln_t

  $if ($MPI)
  interface bcast
     module procedure bcast_0_
  end interface bcast
  interface bcast_alloc
     module procedure bcast_alloc_0_
  end interface bcast_alloc
  $endif

  ! Access specifiers

  private

  public :: soln_t
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  function soln_t_ (gr, omega, discrim) result (sl)

    type(grid_t), intent(in)  :: gr
    complex(WP), intent(in)   :: omega
    type(c_ext_t), intent(in) :: discrim
    type(soln_t)              :: sl

    ! Construct the soln_t

    sl%gr = gr

    sl%s_i = gr%s_i()
    sl%s_o = gr%s_o()

    sl%n_k = gr%n_k

    allocate(sl%ss(sl%s_i:sl%s_o))

    sl%omega = omega
    sl%discrim = discrim

    ! Finish

    return

  end function soln_t_

  !****

  $if ($MPI)

  subroutine bcast_0_ (sl, root_rank)

    type(soln_t), intent(inout) :: sl
    integer, intent(in)         :: root_rank

    ! Broadcast the soln_t

    call bcast(sl%gr, root_rank)
    call bcast_alloc(sl%ss, root_rank)

    call bcast(sl%omega, root_rank)
    call bcast(sl%discrim, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST_ALLOC(type(soln_t),0)

  $endif

  !****

  subroutine set_y (this, i, y, dy_dx)

    class(soln_t), intent(inout) :: this
    integer, intent(in)          :: i
    complex(WP), intent(in)      :: y(:)
    complex(WP), intent(in)      :: dy_dx(:)

    integer :: s
    integer :: k_i
    integer :: k_o

    ! Set the data for y(i)

    $CHECK_BOUNDS(SIZE(y),this%n_k)
    $CHECK_BOUNDS(SIZE(dy_dx),this%n_k)

    seg_loop : do s = this%s_i, this%s_o

       k_i = this%gr%k_i(s)
       k_o = this%gr%k_o(s)

       associate (pt => this%gr%pt)
         call this%ss(s)%set_y(i, pt(k_i:k_o)%x, y(k_i:k_o), dy_dx(k_i:k_o))
       end associate

    end do seg_loop

    ! Finish

    return

  end subroutine set_y

  !****

  function y_1_ (this, i, pt) result (y)

    class(soln_t), intent(in) :: this
    integer, intent(in)       :: i
    type(point_t), intent(in) :: pt
    complex(WP)               :: y

    $ASSERT_DEBUG(pt%s >= this%s_i,Invalid segment)
    $ASSERT_DEBUG(pt%s <= this%s_o,Invalid segment)

    ! Evaluate y(i)

    y = this%ss(pt%s)%y(i, pt%x)

    ! Finish

    return

  end function y_1_
  
  !****

  function y_v_ (this, i, pt) result (y)

    class(soln_t), intent(in) :: this
    integer, intent(in)       :: i
    type(point_t), intent(in) :: pt(:)
    complex(WP)               :: y(SIZE(pt))

    integer :: j

    ! Evaluate y(i)

    !$OMP PARALLEL DO
    do j = 1, SIZE(pt)
       y(j) = this%y(i, pt(j))
    end do

    ! Finish

    return

  end function y_v_

  !****

  function grid (this) result (gr)

    class(soln_t), intent(in) :: this
    type(grid_t)              :: gr

    ! Return the grid

    gr = this%gr

    ! Finish

    return

  end function grid

end module gyre_soln
