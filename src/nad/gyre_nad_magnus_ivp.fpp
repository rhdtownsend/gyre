! Module   : gyre_nad_magnus_ivp
! Purpose  : initial-value solvers (Magnus method, nonadiabatic)
!
! Copyright 2013-2015 Rich Townsend
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

module gyre_nad_magnus_ivp

  ! Uses

  use core_kinds

  use gyre_eqns
  use gyre_ext
  use gyre_magnus_ivp
  use gyre_model

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_magnus_ivp_t) ::  nad_magnus_ivp_t
     private
     class(model_t), pointer      :: ml => null()
     class(c_eqns_t), allocatable :: eq_
   contains
     private
     procedure, public :: shoot => shoot_
  end type nad_magnus_ivp_t

  ! Interfaces

  interface nad_magnus_ivp_t
     module procedure nad_magnus_ivp_t_
  end interface nad_magnus_ivp_t

  ! Access specifiers

  private

  public :: nad_magnus_ivp_t

  ! Procedures

contains

  function nad_magnus_ivp_t_ (ml, eq, scheme) result (iv)

    class(model_t), pointer, intent(in) :: ml
    class(c_eqns_t), intent(in)         :: eq
    character(*), intent(in)            :: scheme
    type(nad_magnus_ivp_t)              :: iv

    ! Construct the nad_magnus_ivp_t

    iv%c_magnus_ivp_t = c_magnus_ivp_t(eq, scheme)

    iv%ml => ml
    allocate(iv%eq_, SOURCE=eq)

    ! Finish

    return

  end function nad_magnus_ivp_t_

!****

  subroutine shoot_ (this, omega, x_a, x_b, E_l, E_r, S)

    class(nad_magnus_ivp_t), intent(in) :: this
    complex(WP), intent(in)             :: omega
    real(WP), intent(in)                :: x_a
    real(WP), intent(in)                :: x_b
    complex(WP), intent(out)            :: E_l(:,:)
    complex(WP), intent(out)            :: E_r(:,:)
    type(c_ext_t), intent(out)          :: S

    logical, parameter :: RESCALE_EIGEN = .TRUE.
!    logical, parameter :: RESCALE_EIGEN = .FALSE.

    real(WP)    :: x
    complex(WP) :: A(this%n_e,this%n_e)
    complex(WP) :: lambda_1
    complex(WP) :: lambda_2
    complex(WP) :: lambda_3

    ! Set up the shooting matrices and scales

    call this%c_magnus_ivp_t%shoot(omega, x_a, x_b, E_l, E_r, S)

    ! Rescale by the uncoupled eigenvalues, in order to help the root
    ! finder

    if (RESCALE_EIGEN) then

       x = 0.5_WP*(x_a + x_b)

       A = this%eq_%A(x, omega)

       lambda_1 = SQRT(A(2,1)*A(1,2))
       lambda_2 = SQRT(A(4,3)*A(3,4))
       lambda_3 = SQRT(A(6,5)*A(5,6))

!       S = S*exp(c_ext_t(-(lambda_1+lambda_2+lambda_3)*(x_b - x_a)))
       S = S*exp(c_ext_t(-(lambda_3)*(x_b - x_a)))

    else

       S = c_ext_t(1._WP)

    endif

    ! Finish

    return

  end subroutine shoot_

end module gyre_nad_magnus_ivp
