! Module   : gyre_status
! Purpose  : status codes for numerical algorithms
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

module gyre_status

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: STATUS_OK = 0
  integer, parameter :: STATUS_OMEGA_DOMAIN = 1
  integer, parameter :: STATUS_ITER_MAX = 2

  ! Interfaces

  private

  public :: STATUS_OK
  public :: STATUS_OMEGA_DOMAIN
  public :: STATUS_ITER_MAX
  public :: status_str

  ! Procedures

contains

  function status_str (status)

    integer, intent(in)       :: status
    character(:), allocatable :: status_str

    ! Return a string representation of the status

    select case (status)
    case (STATUS_OK)
       status_str = 'OK'
    case (STATUS_OMEGA_DOMAIN)
       status_str = 'out-of-domain frequency'
    case (STATUS_ITER_MAX)
       status_str = 'maximum iterations exceeded'
    case default
       status_str = 'unrecognized status'
    end select

    ! Finish

    return

  end function status_str
  
end module gyre_status
