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

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: seg_indices

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

    mask = x(:n_k-1) == x(2:)

    n_s = COUNT(mask)

    k_i = [1,PACK([(k+1,k=1,n_k)], mask)]
    k_o = [PACK([(k,k=1,n_k)], mask),n_k]

    ! Finish

    return

  end subroutine seg_indices

end module gyre_model_util
