! Module   : gyre_r_discrim_func
! Purpose  : discriminant function (real)
!
! Copyright 2013-2017 Rich Townsend
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

module gyre_r_discrim_func

  ! Uses

  use core_kinds

  use gyre_bvp
  use gyre_ext
  use gyre_ext_func
  use gyre_status
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_ext_func_t) :: r_discrim_func_t
     private
     class(r_bvp_t), pointer       :: bp
     class(r_state_t), allocatable :: st
     real(WP)                      :: omega_min
     real(WP)                      :: omega_max
   contains 
     private
     procedure, public :: eval => eval_
  end type r_discrim_func_t

  ! Interfaces

  interface r_discrim_func_t
     module procedure r_discrim_func_t_
  end interface r_discrim_func_t

  ! Access specifiers

  private

  public :: r_discrim_func_t

  ! Procedures

contains

  function r_discrim_func_t_ (bp, st, omega_min, omega_max) result (df)

    class(r_bvp_t), pointer, intent(in) :: bp
    class(r_state_t), intent(in)        :: st
    real(WP), intent(in)                :: omega_min
    real(WP), intent(in)                :: omega_max
    type(r_discrim_func_t)              :: df

    ! Construct the r_discrim_func_t

    df%bp => bp

    allocate(df%st, SOURCE=st)

    df%omega_min = omega_min
    df%omega_max = omega_max

    ! Finish

    return

  end function r_discrim_func_t_

  !****

  subroutine eval_ (this, rx, f_rx, status)

    class(r_discrim_func_t), intent(inout) :: this
    type(r_ext_t), intent(in)              :: rx
    type(r_ext_t), intent(out)             :: f_rx
    integer, intent(out)                   :: status

    real(WP) :: omega

    ! Evaluate the discriminant function

    omega = real(rx)

    if (omega >= this%omega_min .AND. omega <= this%omega_max) then

       this%st%omega = omega

       call this%bp%build(this%st)
       call this%bp%factor()

       f_rx = this%bp%det()

       status = STATUS_OK

    else

       status = STATUS_OMEGA_DOMAIN

    endif
    
    ! Finish

    return

  end subroutine eval_

end module gyre_r_discrim_func
