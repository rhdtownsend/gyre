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
     procedure, public :: k_i
     procedure, public :: k_o
     procedure, public :: s_i
     procedure, public :: s_o
     procedure, public :: s_x
     procedure, public :: x_i
     procedure, public :: x_o
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
    logical                    :: mask(gr_base%n_k)
    integer                    :: n_mask
    type(point_t), allocatable :: pt_int(:)
    type(point_t)              :: pt_i
    type(point_t)              :: pt_o

    ! Construct a grid_t as a nesting of gr_base, with limits
    ! (x_i,x_o)

    ! First, restrict the limits to the range of gr_base

    x_i_ = MAX(x_i, gr_base%pt(1)%x)
    x_o_ = MIN(x_o, gr_base%pt(gr_base%n_k)%x)

    ! Set up a mask for the points from gr_base to include (excluding
    ! boundary points)

    mask = gr_base%pt%x > x_i_ .AND. gr_base%pt%x < x_o_

    n_mask = COUNT(mask)

    ! Set up the internal and boundary points

    if (n_mask > 0) then

       pt_int = PACK(gr_base%pt, mask)
       
       pt_i = point_t(pt_int(1)%s, x_i_)
       pt_o = point_t(pt_int(n_mask)%s, x_o_)

    else

       allocate(pt_int(0))

       pt_i = point_t(gr_base%s_x(x_i_, back=.TRUE.), x_i_)
       pt_o = point_t(gr_base%s_x(x_o_, back=.FALSE.), x_o_)

    end if

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

  function k_i (this, s)

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

  end function k_i

  !****

  function k_o (this, s)

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

  end function k_o

  !****

  function s_i (this)

    class(grid_t), intent(in) :: this
    integer                   :: s_i

    ! Return the innermost segment of the grid

    s_i = this%pt(1)%s

    ! Finish

    return

  end function s_i

  !****

  function s_o (this)

    class(grid_t), intent(in) :: this
    integer                   :: s_o

    ! Return the outermost segment of the grid

    s_o = this%pt(this%n_k)%s

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
    ! is inside-out

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

       k_i = this%k_i(s)
       k_o = this%k_o(s)
         
       if (x >= this%pt(k_i)%x .AND. x <= this%pt(k_o)%x) exit seg_loop
       
    end do seg_loop

    ! Finish

    return

  end function s_x

  !****

  function x_i (this)

    class(grid_t), intent(in) :: this
    real(WP)                  :: x_i

    ! Return the innermost abscissa

    x_i = this%pt(1)%x

    ! Finish

    return

  end function x_i

  !****

  function x_o (this)

    class(grid_t), intent(in) :: this
    real(WP)                  :: x_o

    ! Return the outermost abscissa

    x_o = this%pt(this%n_k)%x

    ! Finish

    return

  end function x_o

end module gyre_grid
