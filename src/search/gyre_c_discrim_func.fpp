! Module   : gyre_c_discrim_func
! Purpose  : discriminant function (complex)
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

module gyre_c_discrim_func

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

  type, extends (c_ext_func_t) :: c_discrim_func_t
     private
     class(c_bvp_t), pointer       :: bp
     class(c_state_t), allocatable :: st
     real(WP)                      :: omega_min
     real(WP)                      :: omega_max
     complex(WP), allocatable      :: omega_def(:)
   contains 
     private
     procedure, public :: eval => eval_
  end type c_discrim_func_t

  ! Interfaces

  interface c_discrim_func_t
     module procedure c_discrim_func_t_
  end interface c_discrim_func_t

  ! Access specifiers

  private

  public :: c_discrim_func_t

  ! Procedures

contains

  function c_discrim_func_t_ (bp, st, omega_min, omega_max, omega_def) result (df)

    class(c_bvp_t), pointer, intent(in) :: bp
    class(c_state_t), intent(in)        :: st
    real(WP), intent(in)                :: omega_min
    real(WP), intent(in)                :: omega_max
    complex(WP), intent(in), optional   :: omega_def(:)
    type(c_discrim_func_t)              :: df

    ! Construct the c_discrim_func_t

    df%bp => bp

    allocate(df%st, SOURCE=st)

    df%omega_min = omega_min
    df%omega_max = omega_max

    if (PRESENT(omega_def)) then
       df%omega_def = omega_def
    else
       allocate(df%omega_def(0))
    endif

    ! Finish

    return

  end function c_discrim_func_t_

  !****

  subroutine eval_ (this, cx, f_cx, status)

    class(c_discrim_func_t), intent(inout) :: this
    type(c_ext_t), intent(in)              :: cx
    type(c_ext_t), intent(out)             :: f_cx
    integer, intent(out)                   :: status

    complex(WP)     :: omega
    type(c_state_t) :: st

    ! Evaluate the discriminant function

    omega = cmplx(cx)

    if (REAL(omega) >= this%omega_min .AND. REAL(omega) <= this%omega_max) then

       this%st%omega = omega

       call this%bp%build(this%st)
       call this%bp%factor()

       f_cx = this%bp%det()

       status = STATUS_OK

       ! Handle root deflation

       f_cx = f_cx*PRODUCT(omega/(omega - this%omega_def))

    else

       status = STATUS_OMEGA_DOMAIN

    endif
    
    ! Finish

    return

  end subroutine eval_

end module gyre_c_discrim_func
