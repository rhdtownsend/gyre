! Module   : gyre_jacobian_rad_dziem
! Purpose  : radial adiabatic Jacobian evaluation (Dziembowski variables)
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

module gyre_jacobian_rad_dziem

  ! Uses

  use core_kinds

  use gyre_jacobian
  use gyre_model
  use gyre_oscpar
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (jacobian_t) :: jacobian_rad_dziem_t
     private
     class(model_t), pointer :: cf => null()
     type(oscpar_t)          :: op
   contains
     private
     procedure, public :: eval => eval_
     procedure, public :: eval_logx => eval_logx_
     procedure, public :: trans_matrix => trans_matrix_
  end type jacobian_rad_dziem_t

  ! Interfaces

  interface jacobian_rad_dziem_t
     module procedure jacobian_rad_dziem_t_
  end interface jacobian_rad_dziem_t

  ! Access specifiers

  private

  public :: jacobian_rad_dziem_t

  ! Procedures

contains

  function jacobian_rad_dziem_t_ (cf, op) result (jc)

    class(model_t), pointer, intent(in) :: cf
    type(oscpar_t), intent(in)          :: op
    type(jacobian_rad_dziem_t)          :: jc

    ! Construct the jacobian_rad_dziem_t

    jc%cf => cf
    jc%op = op

    jc%n_e = 2

    ! Finish

    return

  end function jacobian_rad_dziem_t_

!****

  subroutine eval_ (this, x, omega, A)

    class(jacobian_rad_dziem_t), intent(in) :: this
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

    class(jacobian_rad_dziem_t), intent(in) :: this
    real(WP), intent(in)                    :: x
    complex(WP), intent(in)                 :: omega
    complex(WP), intent(out)                :: A(:,:)
    
    $CHECK_BOUNDS(SIZE(A, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(A, 2),this%n_e)

    ! Evaluate the log(x)-space Jacobian matrix

    associate(V_g => this%cf%V(x)/this%cf%Gamma_1(x), U => this%cf%U(x), &
              As => this%cf%As(x), c_1 => this%cf%c_1(x), &
              omega_c => this%cf%omega_c(x, this%op%m, omega))

      A(1,1) = V_g - 1._WP
      A(1,2) = -V_g
      
      A(2,1) = c_1*omega_c**2 + U - As
      A(2,2) = As - U + 3._WP

    end associate

    ! Finish

    return

  end subroutine eval_logx_

!****

  function trans_matrix_ (this, x, omega, to_canon) result (M)

    class(jacobian_rad_dziem_t), intent(in) :: this
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

    ! Calculate the transformation matrix to convert variables between the
    ! canonical formulation and the Dziembowski formulation

    if (to_canon) then
       M = identity_matrix(this%n_e)
    else
       M = identity_matrix(this%n_e)
    endif

    ! Finish

    return

  end function trans_matrix_

end module gyre_jacobian_rad_dziem
