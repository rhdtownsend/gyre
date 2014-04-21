! Module   : gyre_nad_dziem_jacobian
! Purpose  : nonadiabatic Jacobian evaluation (DZIEM variables)
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

module gyre_nad_dziem_jacobian

  ! Uses

  use core_kinds

  use gyre_jacobian
  use gyre_model
  use gyre_modepar
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(jacobian_t) :: nad_dziem_jacobian_t
     private
     class(model_t), pointer :: ml => null()
     type(modepar_t)         :: mp
   contains
     private
     procedure, public :: eval => eval_
     procedure, public :: eval_logx => eval_logx_
     procedure, public :: trans_matrix => trans_matrix_
  end type nad_dziem_jacobian_t

  ! Interfaces

  interface nad_dziem_jacobian_t
     module procedure nad_dziem_jacobian_t_
  end interface nad_dziem_jacobian_t

  ! Access specifiers

  private

  public :: nad_dziem_jacobian_t

  ! Procedures

contains

  function nad_dziem_jacobian_t_ (ml, mp) result (jc)

    class(model_t), pointer, intent(in) :: ml
    type(modepar_t), intent(in)         :: mp
    type(nad_dziem_jacobian_t)          :: jc

    ! Construct the nad_dziem_jacobian_t

    jc%ml => ml
    jc%mp = mp

    jc%n_e = 6

    ! Finish

    return

  end function nad_dziem_jacobian_t_

!****

  subroutine eval_ (this, x, omega, A)

    class(nad_dziem_jacobian_t), intent(in) :: this
    real(WP), intent(in)                    :: x
    complex(WP), intent(in)                 :: omega
    complex(WP), intent(out)                :: A(:,:)
    
    ! Evaluate the Jacobian matrix

    call this%eval_logx(x, omega, A)

    A = A/x

    ! Finish

    return

  end subroutine eval_

!****

  subroutine eval_logx_ (this, x, omega, A)

    class(nad_dziem_jacobian_t), intent(in) :: this
    real(WP), intent(in)                    :: x
    complex(WP), intent(in)                 :: omega
    complex(WP), intent(out)                :: A(:,:)

    $CHECK_BOUNDS(SIZE(A, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(A, 2),this%n_e)

    ! Evaluate the log(x)-space Jacobian matrix

    associate(V => this%ml%V(x), V_g => this%ml%V(x)/this%ml%Gamma_1(x), &
              U => this%ml%U(x), As => this%ml%As(x), c_1 => this%ml%c_1(x), &
              nabla_ad => this%ml%nabla_ad(x), delta => this%ml%delta(x), &
              nabla => this%ml%nabla(x), &
              c_rad => this%ml%c_rad(x), dc_rad => this%ml%dc_rad(x), &
              c_thm => this%ml%c_thm(x), c_dif => this%ml%c_dif(x), &
              c_eps_ad => this%ml%c_eps_ad(x), c_eps_S => this%ml%c_eps_S(x), &
              kappa_ad => this%ml%kappa_ad(x), kappa_S => this%ml%kappa_S(x), &
              l => this%mp%l, omega_c => this%ml%omega_c(x, this%mp%m, omega))

      A(1,1) = V_g - 1._WP - l
      A(1,2) = l*(l+1)/(c_1*omega_c**2) - V_g
      A(1,3) = V_g
      A(1,4) = 0._WP
      A(1,5) = delta
      A(1,6) = 0._WP

      A(2,1) = c_1*omega_c**2 - As
      A(2,2) = As - U + 3._WP - l
      A(2,3) = -As
      A(2,4) = 0._WP
      A(2,5) = delta
      A(2,6) = 0._WP

      A(3,1) = 0._WP
      A(3,2) = 0._WP
      A(3,3) = 3._WP - U - l
      A(3,4) = 1._WP
      A(3,5) = 0._WP
      A(3,6) = 0._WP

      A(4,1) = U*As
      A(4,2) = U*V_g
      A(4,3) = l*(l+1) - U*V_g
      A(4,4) = -U - l + 2._WP
      A(4,5) = -U*delta
      A(4,6) = 0._WP

      A(5,1) = V*(nabla_ad*(U - c_1*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif)
      A(5,2) = V*(l*(l+1)/(c_1*omega_c**2)*(nabla_ad - nabla) - c_dif)
      A(5,3) = V*c_dif
      A(5,4) = V*nabla_ad
      A(5,5) = V*nabla*(4._WP - kappa_S) - (l - 2._WP)
      A(5,6) = -V*nabla/c_rad

      A(6,1) = l*(l+1)*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
      A(6,2) = V*c_eps_ad - l*(l+1)*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega_c**2))
      A(6,3) = l*(l+1)*nabla_ad/nabla*c_rad - V*c_eps_ad
      A(6,4) = 0._WP
      A(6,5) = c_eps_S - l*(l+1)*c_rad/(nabla*V) - (0._WP,1._WP)*omega_c*c_thm
      A(6,6) = -1._WP - l

    end associate

    ! Finish

    return

  end subroutine eval_logx_

!****

  function trans_matrix_ (this, x, omega, to_canon) result (M)

    class(nad_dziem_jacobian_t), intent(in) :: this
    real(WP), intent(in)                    :: x
    complex(WP), intent(in)                 :: omega
    logical, intent(in)                     :: to_canon
    $if ($GFORTRAN_PR_58007)
    complex(WP), allocatable                :: M(:,:)
    $else
    complex(WP)                             :: M(this%n_e,this%n_e)
    $endif

    $if ($GFORTRAN_PR_58007)
    allocate(M(this%n_e,this%n_e))
    $endif

    ! Calculate the transformation matrix to convert variables between
    ! the canonical formulation and the Dziembowski formulation

    if (to_canon) then
       M = identity_matrix(this%n_e)
    else
       M = identity_matrix(this%n_e)
    endif

    ! Finish

    return

  end function trans_matrix_

end module gyre_nad_dziem_jacobian
