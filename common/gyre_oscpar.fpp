! Module   : gyre_ad_oscpar
! Purpose  : adiabatic oscillation parameters
!
! Copyright 2013 Rich Townsend
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

module gyre_ad_oscpar

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: ad_oscpar_t
     private
     real(WP), public              :: lambda_0
     integer, public               :: l
!    character(LEN=:), allocatable :: outer_bound_type
     character(LEN=256), public    :: outer_bound_type
   contains
     private
     procedure, public :: init
  end type ad_oscpar_t

 ! Access specifiers

  private

  public :: ad_oscpar_t

  ! Procedures

contains

  subroutine init (this, l, outer_bound_type)

    class(ad_oscpar_t), intent(out) :: this
    integer, intent(in)             :: l
    character(LEN=*), intent(in)    :: outer_bound_type

    ! Initialize the ad_oscpar

!    if(l == 0) then
!       this%lambda_0 = 0._WP
!    else
       this%lambda_0 = l - 2._WP
!    endif

    this%l = l

    this%outer_bound_type = outer_bound_type

    ! Finish

    return

  end subroutine init

end module gyre_ad_oscpar
