! Module   : gyre_out_util
! Purpose  : output utility routines
!
! Copyright 2020 Rich Townsend & The GYRE Team
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
  
module gyre_out_util

  ! Uses

  use core_kinds

  use gyre_util
  use gyre_wave

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: filter_wave

contains

  function filter_wave (wv, filter_list)

    class(wave_t), intent(in) :: wv
    character(*), intent(in)  :: filter_list
    logical                   :: filter_wave

    character(LEN(filter_list)), allocatable :: filters(:)
    integer                                  :: i

    ! Decide whether to filter the wave

    filters = split_list(filter_list, ',')

    filter_wave = .FALSE.

    filter_loop : do i = 1, SIZE(filters)

       select case (filters(i))
       case ('stable')
          filter_wave = filter_wave .OR. AIMAG(wv%omega) <= 0._WP
       case ('unstable')
          filter_wave = filter_wave .OR. AIMAG(wv%omega) > 0._WP
       case default
          $ABORT(Unrecognized filter in filter list)
       end select

    end do filter_loop

    ! Finish

    return

  end function filter_wave

end module gyre_out_util
