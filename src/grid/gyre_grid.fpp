! Module   : gyre_grid
! Purpose  : segmented grids
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_grid

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_point

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: grid_t
     type(point_t), allocatable :: pt(:)
     integer                    :: n_k
   contains
     private
     procedure, public :: pt_i
     procedure, public :: pt_o
     procedure, public :: pt_x
     procedure, public :: s_i
     procedure, public :: s_o
     procedure, public :: s_x
     procedure, public :: x_i
     procedure, public :: x_o
     procedure, public :: k_s_i
     procedure, public :: k_s_o
  end type grid_t

  ! Interfaces

  interface grid_t
     module procedure grid_t_x_
     module procedure grid_t_nest_
     module procedure grid_t_resamp_
  end interface grid_t

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

  public :: grid_t
  $if ($MPI)
  public :: bcast
  public :: bcast_alloc
  $endif

  ! Procedures

contains

  function grid_t_x_ (x) result (gr)

    real(WP), intent(in) :: x(:)
    type(grid_t)         :: gr

    integer :: n_k
    integer :: k

    ! Construct a grid_t from the input abscissae x (with segment
    ! boundaries delineated by double points)

    n_k = SIZE(x)

    if (n_k > 0) then

       $ASSERT_DEBUG(ALL(x(2:) >= x(:n_k-1)),Non-monotonic data)
       
       allocate(gr%pt(n_k))

       gr%pt(1)%x = x(1)
       gr%pt(1)%s = 1

       do k = 2, n_k

          gr%pt(k)%x = x(k)

          if (x(k) == x(k-1)) then
             gr%pt(k)%s = gr%pt(k-1)%s + 1
          else
             gr%pt(k)%s = gr%pt(k-1)%s
          endif

       end do

    end if

    gr%n_k = n_k

    ! Finish

    return

  end function grid_t_x_

  !****

  function grid_t_nest_ (gr_base, x_i, x_o) result (gr)

    type(grid_t), intent(in) :: gr_base
    real(WP), intent(in)     :: x_i
    real(WP), intent(in)     :: x_o
    type(grid_t)             :: gr

    real(WP)                   :: x_i_
    real(WP)                   :: x_o_
    type(point_t)              :: pt_i
    type(point_t)              :: pt_o
    type(point_t), allocatable :: pt_int(:)

    ! Construct a grid_t as a nesting of gr_base, with limits
    ! (x_i,x_o)

    ! First, restrict the limits to the range of gr_base

    x_i_ = MAX(x_i, gr_base%x_i())
    x_o_ = MIN(x_o, gr_base%x_o())

    ! Set up inner and outer points

    pt_i = gr_base%pt_x(x_i_)
    pt_o = gr_base%pt_x(x_o_, back=.TRUE.)

    ! Select internal points from gr_base to include (excluding
    ! boundary points)

    pt_int = PACK(gr_base%pt, MASK=(gr_base%pt%x > x_i_ .AND. gr_base%pt%x < x_o_))

    ! Create the grid_t

    gr%pt = [pt_i,pt_int,pt_o]
    gr%n_k = SIZE(gr%pt)

    ! Finish

    return

  end function grid_t_nest_

  !****

  function grid_t_resamp_ (gr_base, dn) result (gr)

    type(grid_t), intent(in) :: gr_base
    integer, intent(in)      :: dn(:)
    type(grid_t)             :: gr

    integer  :: n_k_base
    integer  :: n_k
    integer  :: k
    integer  :: j
    integer  :: i
    real(WP) :: w

    $CHECK_BOUNDS(SIZE(dn),gr_base%n_k-1)

    ! Construct a grid_t by resampling gr_base, with dn additional
    ! points placed uniformly in each cell

    n_k_base = gr_base%n_k
    n_k = n_k_base + SUM(dn)

    allocate(gr%pt(n_k))

    k = 1

    cell_loop : do j = 1, n_k_base-1

       associate (pt_a => gr_base%pt(j), &
                  pt_b => gr_base%pt(j+1))

         if (pt_a%s == pt_b%s) then

            do i = 1, dn(j)+1

               w = REAL(i-1, WP)/REAL(dn(j)+1, WP)

               gr%pt(k)%s = pt_a%s
               gr%pt(k)%x = (1._WP-w)*pt_a%x + w*pt_b%x

               k = k + 1

            end do

         else

            $ASSERT(dn(j) == 0,Attempt to add points at cell boundary)

            gr%pt(k)%s = pt_a%s
            gr%pt(k)%x = pt_a%x

            k = k + 1

         endif

       end associate

    end do cell_loop

    gr%pt(k) = gr_base%pt(n_k_base)

    gr%n_k = n_k

    ! Finish

    return

  end function grid_t_resamp_

  !****

  $if ($MPI)

  subroutine bcast_0_ (gr, root_rank)

    type(grid_t), intent(inout) :: gr
    integer, intent(in)         :: root_rank

    ! Broadcast the grid_t

    call bcast_alloc(gr%pt, root_rank)

    call bcast(gr%n_k, root_rank)

    ! Finish

    return

  end subroutine bcast_0_

  $BCAST_ALLOC(type(grid_t),0)

  $endif

  !****

  function pt_i (this)

    class(grid_t), intent(in) :: this
    type(point_t)             :: pt_i

    ! Return the innermost point of the grid

    pt_i = this%pt(1)

    ! Finish

    return

  end function pt_i

  !****

  function pt_o (this)

    class(grid_t), intent(in) :: this
    type(point_t)             :: pt_o

    ! Return the outermost point of the grid

    pt_o = this%pt(this%n_k)

    ! Finish

    return

  end function pt_o

  !****

  function pt_x (this, x, back) result (pt)

    class(grid_t), intent(in)     :: this
    real(WP), intent(in)          :: x
    logical, intent(in), optional :: back
    type(point_t)                 :: pt

    ! Return a point containing the abcissa x. If back is present and
    ! .TRUE., the segment search is done outside-in; otherwise, it is
    ! inside-out (see s_x)

    pt%s = this%s_x(x, back)
    pt%x = x

    ! Finish

    return

  end function pt_x

  !****

  function s_i (this)

    class(grid_t), intent(in) :: this
    integer                   :: s_i

    type(point_t) :: pt_i

    ! Return the innermost segment of the grid

    pt_i = this%pt_i()

    s_i = pt_i%s

    ! Finish

    return

  end function s_i

  !****

  function s_o (this)

    class(grid_t), intent(in) :: this
    integer                   :: s_o

    type(point_t) :: pt_o

    ! Return the outermost segment of the grid

    pt_o = this%pt_o()

    s_o = pt_o%s

    ! Finish

    return

  end function s_o

  !****

  function s_x (this, x, back) result (s)

    class(grid_t), intent(in)     :: this
    real(WP), intent(in)          :: x
    logical, intent(in), optional :: back
    integer                       :: s

    logical :: back_
    integer :: s_a
    integer :: s_b
    integer :: ds
    integer :: k_i
    integer :: k_o

    if (PRESENT(back)) then
       back_ = back
    else
       back_ = .FALSE.
    endif

    ! Locate the segment which brackets the abcissa x. If back is
    ! present and .TRUE., the search is done outside-in; otherwise, it
    ! is inside-out. If x is not inside the grid, then return a
    ! segment index just above/below the maximum/minimum segment of
    ! the grid

    if (back_) then
       s_a = this%pt(this%n_k)%s
       s_b = this%pt(1)%s
       ds = -1
    else
       s_a = this%pt(1)%s
       s_b = this%pt(this%n_k)%s
       ds = 1
    endif

    seg_loop : do s = s_a, s_b, ds

       k_i = this%k_s_i(s)
       k_o = this%k_s_o(s)
         
       if (x >= this%pt(k_i)%x .AND. x <= this%pt(k_o)%x) exit seg_loop
       
    end do seg_loop

    ! Finish

    return

  end function s_x

  !****

  function x_i (this)

    class(grid_t), intent(in) :: this
    real(WP)                  :: x_i

    type(point_t) :: pt_i

    ! Return the innermost abscissa of the grid

    pt_i = this%pt_i()

    x_i = pt_i%x

    ! Finish

    return

  end function x_i

  !****

  function x_o (this)

    class(grid_t), intent(in) :: this
    real(WP)                  :: x_o

    type(point_t) :: pt_o

    ! Return the outermost abscissa of the grid

    pt_o = this%pt_o()

    x_o = pt_o%x

    ! Finish

    return

  end function x_o

  !****

  function k_s_i (this, s) result (k_i)

    class(grid_t), intent(in) :: this
    integer, intent(in)       :: s
    integer                   :: k_i

    $ASSERT_DEBUG(s >= this%s_i(),Invalid segment)
    $ASSERT_DEBUG(s <= this%s_o(),Invalid segment)

    ! Return the index of the innermost point in segment s

    do k_i = 1, this%n_k
       if (this%pt(k_i)%s == s) exit
    end do

    ! Finish

    return

  end function k_s_i

  !****

  function k_s_o (this, s) result (k_o)

    class(grid_t), intent(in) :: this
    integer, intent(in)       :: s
    integer                   :: k_o

    $ASSERT_DEBUG(s >= this%s_i(),Invalid segment)
    $ASSERT_DEBUG(s <= this%s_o(),Invalid segment)

    ! Return the index of the outermost point in segment s

    do k_o = this%n_k, 1, -1
       if (this%pt(k_o)%s == s) exit
    end do

    ! Finish

    return

  end function k_s_o

end module gyre_grid
