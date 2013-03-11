! Module   : gyre_ad_jacobian
! Purpose  : adiabatic Jacobian evaluation
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

module gyre_ad_jacobian

  ! Uses

  use core_kinds

  use gyre_jacobian
  use gyre_mech_coeffs
  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(jacobian_t) :: ad_jacobian_t
     private
     class(mech_coeffs_t), pointer :: mc => null()
     class(oscpar_t), pointer      :: op => null()
   contains
     private
     procedure, public :: init
     procedure, public :: eval
     procedure, public :: eval_logx
  end type ad_jacobian_t

  ! Access specifiers

  private

  public :: ad_jacobian_t

  ! Procedures

contains

  subroutine init (this, mc, op)

    class(ad_jacobian_t), intent(out)        :: this
    class(mech_coeffs_t), intent(in), target :: mc
    class(oscpar_t), intent(in), target      :: op

    ! Initialize the ad_jacobian

    this%mc => mc
    this%op => op

    this%n_e = 4

    ! Finish

    return

  end subroutine init

!****

  subroutine eval (this, omega, x, A)

    class(ad_jacobian_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x
    complex(WP), intent(out)         :: A(:,:)
    
    ! Evaluate the Jacobian matrix

    call this%eval_logx(omega, x, A)

    A = A/x

    ! Finish

    return

  end subroutine eval

!****

  subroutine eval_logx (this, omega, x, A)

    class(ad_jacobian_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x
    complex(WP), intent(out)         :: A(:,:)
    
    $CHECK_BOUNDS(SIZE(A, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(A, 2),this%n_e)

    ! Evaluate the log(x)-space Jacobian matrix

    associate(V_g => this%mc%V(x)/this%mc%Gamma_1(x), U => this%mc%U(x), &
              As => this%mc%As(x), c_1 => this%mc%c_1(x), &
              lambda_0 => this%op%lambda_0, l => this%op%l)

      A(1,1) = V_g - 3._WP - lambda_0
      A(1,2) = l*(l+1)/(c_1*omega**2) - V_g
      A(1,3) = V_g
      A(1,4) = 0._WP
      
      A(2,1) = c_1*omega**2 - As
      A(2,2) = As - U + 1._WP - lambda_0
      A(2,3) = -As
      A(2,4) = 0._WP
      
      A(3,1) = 0._WP
      A(3,2) = 0._WP
      A(3,3) = 1._WP - U - lambda_0
      A(3,4) = 1._WP
      
      A(4,1) = U*As
      A(4,2) = U*V_g
      A(4,3) = l*(l+1) - U*V_g
      A(4,4) = -U - lambda_0

    end associate

    ! Finish

    return

  end subroutine eval_logx

end module gyre_ad_jacobian
