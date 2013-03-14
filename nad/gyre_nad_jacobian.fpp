! Module   : gyre_nad_jacobian
! Purpose  : nonadiabatic Jacobian evaluation
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

module gyre_nad_jacobian

  ! Uses

  use core_kinds

  use gyre_jacobian
  use gyre_mech_coeffs
  use gyre_therm_coeffs
  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(jacobian_t) :: nad_jacobian_t
     private
     class(mech_coeffs_t), pointer  :: mc => null()
     class(therm_coeffs_t), pointer :: tc => null()
     class(oscpar_t), pointer       :: op => null()
   contains
     private
     procedure, public :: init
     procedure, public :: eval
     procedure, public :: eval_logx
  end type nad_jacobian_t

  ! Access specifiers

  private

  public :: nad_jacobian_t

  ! Procedures

contains

  subroutine init (this, mc, tc, op)

    class(nad_jacobian_t), intent(out)        :: this
    class(mech_coeffs_t), intent(in), target  :: mc
    class(therm_coeffs_t), intent(in), target :: tc
    class(oscpar_t), intent(in), target       :: op

    ! Initialize the nad_jacobian

    this%mc => mc
    this%tc => tc
    this%op => op

    this%n_e = 6

    ! Finish

    return

  end subroutine init

!****

  subroutine eval (this, omega, x, A)

    class(nad_jacobian_t), intent(in) :: this
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

    class(nad_jacobian_t), intent(in) :: this
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x
    complex(WP), intent(out)         :: A(:,:)

    real(WP) :: V_x2
    real(WP) :: c_2
    
    $CHECK_BOUNDS(SIZE(A, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(A, 2),this%n_e)

    ! Evaluate the log(x)-space Jacobian matrix

    associate(V=> this%mc%V(x), V_g => this%mc%V(x)/this%mc%Gamma_1(x), U => this%mc%U(x), &
              As => this%mc%As(x), c_1 => this%mc%c_1(x), &
              c_rad => this%tc%c_rad(x), dc_rad => this%tc%dc_rad(x), &
              c_gen => this%tc%c_gen(x), c_thm => this%tc%c_thm(x), &
              nabla => this%tc%nabla(x), nabla_ad => this%tc%nabla_ad(x), dnabla_ad => this%tc%dnabla_ad(x), &
              alpha_T => this%tc%alpha_T(x), &
              kappa_ad => this%tc%kappa_ad(x), kappa_S => this%tc%kappa_S(x), &
              epsilon_ad => this%tc%epsilon_ad(x), epsilon_S => this%tc%epsilon_S(x), &
              lambda_0 => this%op%lambda_0, l => this%op%l)

      if(x /= 0._WP) then
         V_x2 = V/x**2
      else
         V_x2 = this%tc%V_x2_0
      endif

      c_2 = (kappa_ad-4._WP*nabla_ad)*V*nabla + nabla_ad*(dnabla_ad+V)

      A(1,1) = V_g - 3._WP - lambda_0
      A(1,2) = l*(l+1)/(c_1*omega**2) - V_g
      A(1,3) = V_g
      A(1,4) = 0._WP
      A(1,5) = alpha_T*x**2
      A(1,6) = 0._WP

      A(2,1) = c_1*omega**2 - As
      A(2,2) = As - U + 1 - lambda_0
      A(2,3) = -As
      A(2,4) = 0._WP
      A(2,5) = alpha_T*x**2
      A(2,6) = 0._WP

      A(3,1) = 0._WP
      A(3,2) = 0._WP
      A(3,3) = 1._WP - U - lambda_0
      A(3,4) = 1._WP
      A(3,5) = 0._WP
      A(3,6) = 0._WP

      A(4,1) = U*As
      A(4,2) = U*V_g
      A(4,3) = l*(l+1) - U*V_g
      A(4,4) = -U - lambda_0
      A(4,5) = -U*alpha_T*x**2
      A(4,6) = 0._WP

      A(5,1) = V_x2*(nabla_ad*(U - c_1*omega**2) - 4._WP*(nabla_ad - nabla) + c_2)
      A(5,2) = V_x2*(l*(l+1)/(c_1*omega**2)*(nabla_ad - nabla) - c_2)
      A(5,3) = V_x2*c_2
      A(5,4) = V_x2*nabla_ad
      A(5,5) = V*nabla*(4._WP - kappa_S) - 2._WP - lambda_0
      A(5,6) = -V_x2*nabla/c_rad

      A(6,1) = l*(l+1)*(nabla_ad/nabla - 1._WP)*c_rad - epsilon_ad*V*c_gen
      A(6,2) = epsilon_ad*V*c_gen - l*(l+1)*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega**2))
      A(6,3) = l*(l+1)*nabla_ad/nabla*c_rad - epsilon_ad*V*c_gen
      A(6,4) = 0._WP
      A(6,5) = epsilon_S*c_gen*x**2 - l*(l+1)*c_rad/(nabla*V_x2) - (0._WP,1._WP)*omega*c_thm*x**2
      A(6,6) = -3._WP - lambda_0

    end associate

    ! Finish

    return

  end subroutine eval_logx

end module gyre_nad_jacobian
