! Module   : gyre_nad_findiff_ivp
! Purpose  : initial-value solvers (finite difference, nonadiabatic)
!
! Copyright 2013-2014 Rich Townsend
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

module gyre_nad_findiff_ivp

  ! Uses

  use core_kinds

  use gyre_ext
  use gyre_findiff_ivp
  use gyre_jacob
  use gyre_model

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_findiff_ivp_t) ::  nad_findiff_ivp_t
     private
     class(model_t), pointer :: ml => null()
   contains
     procedure, public :: shoot => shoot_
  end type nad_findiff_ivp_t

  ! Interfaces

  interface nad_findiff_ivp_t
     module procedure nad_findiff_ivp_t_
  end interface nad_findiff_ivp_t

  ! Access specifiers

  private

  public :: nad_findiff_ivp_t

  ! Procedures

contains

  function nad_findiff_ivp_t_ (ml, jc) result (iv)

    class(model_t), pointer, intent(in) :: ml
    class(c_jacob_t), intent(in)        :: jc
    type(nad_findiff_ivp_t)             :: iv

    ! Construct the nad_findiff_ivp_t

    iv%c_findiff_ivp_t = c_findiff_ivp_t(jc)

    iv%ml => ml

    ! Finish

    return
    
  end function nad_findiff_ivp_t_

!****

  subroutine shoot_ (this, omega, x_a, x_b, E_l, E_r, S)

    class(nad_findiff_ivp_t), intent(in) :: this
    complex(WP), intent(in)              :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    complex(WP), intent(out)             :: E_l(:,:)
    complex(WP), intent(out)             :: E_r(:,:)
    type(c_ext_t), intent(out)           :: S

    real(WP), allocatable :: x(:)
    real(WP)              :: w(6)

    ! Set up the shooting matrices and scales

    x = this%abscissa(x_a, x_b)

    if (ANY(this%ml%c_thm(x) > 1.E4*this%ml%c_rad(x))) then
       w = [0.5_WP,0.5_WP,0.5_WP,0.5_WP,1._WP,0._WP]
    else
       w = [0.5_WP,0.5_WP,0.5_WP,0.5_WP,0.5_WP,0.5_WP]
    endif

    call this%c_findiff_ivp_t%shoot_w(w, omega, x_a, x_b, E_l, E_r, S)

    ! Finish

    return

  end subroutine shoot_

end module gyre_nad_findiff_ivp
