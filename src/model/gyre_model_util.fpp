! Module   : gyre_model_util
! Purpose  : stellar model utilities
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

module gyre_model_util

  ! Uses

  use core_kinds

  use gyre_model_par
  use gyre_util

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: seg_indices
  public :: snap_points

  ! Procedures

contains

  subroutine seg_indices (x, k_i, k_o)

    real(WP), intent(in)              :: x(:)
    integer, allocatable, intent(out) :: k_i(:)
    integer, allocatable, intent(out) :: k_o(:)

    integer :: n_k
    logical :: mask(SIZE(x)-1)
    integer :: n_s
    integer :: k

    ! Partition the monotonic array x into strictly-monotonic
    ! segments, by splitting at duplicate points; return the index
    ! range of the segments in k_i/k_o

    n_k = SIZE(x)

    $ASSERT_DEBUG(ALL(x(2:) >= x(:n_k-1)),Non-monotonic data)
    
    mask = x(:n_k-1) == x(2:)

    n_s = COUNT(mask)

    k_i = [1,PACK([(k+1,k=1,n_k)], mask)]
    k_o = [PACK([(k,k=1,n_k)], mask),n_k]

    ! Finish

    return

  end subroutine seg_indices

  !****

  subroutine snap_points (dx_snap, x, m)

    real(WP), intent(in)              :: dx_snap
    real(WP), intent(inout)           :: x(:)
    real(WP), optional, intent(inout) :: m(:)

    integer  :: i
    real(WP) :: x_snap
    real(WP) :: m_snap

    if (PRESENT(m)) then
       $CHECK_BOUNDS(SIZE(m),SIZE(x))
    endif

    ! Snap model points to fix possible numerical issues

    ! Central point

    if (x(1) < dx_snap) then

       x(1) = 0._WP

       if (PRESENT(m)) then
          m(1) = 0._WP
       endif

       if (check_log_level('INFO')) then
          write(OUTPUT_UNIT, 100) 'Snapping central point to x=0'
100       format(3X,A)
       endif

    endif

    ! Other points

    snap_loop : do i = 2, SIZE(x)-1

       if (x(i+1) - x(i) > 0._WP .AND. x(i+1) - x(i) < dx_snap) then
          
          x_snap = 0.5_WP*(x(i+1) + x(i))
          x(i:i+1) = x_snap

          if (PRESENT(m)) then
             m_snap = 0.5_WP*(m(i+1) + m(i))
             m(i:i+1) = m_snap
          endif

          if (check_log_level('INFO')) then
             write(OUTPUT_UNIT, 110) 'Snapping points', i, 'and', i+1, 'to x=', x_snap
110          format(3X,A,1X,I0,1X,A,1X,I0,1X,A,F6.4)
          endif
             
       end if

    end do snap_loop

    ! Finish

    return

  end subroutine snap_points

end module gyre_model_util
