! Module   : gyre_jacobian_rad_jcd
! Purpose  : radial adiabatic Jacobian evaluation (JCD variables)
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

module gyre_jacobian_rad_jcd

  ! Uses

  use core_kinds

  use gyre_jacobian
  use gyre_coeffs
  use gyre_oscpar

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (jacobian_t) :: jacobian_rad_jcd_t
     private
     class(coeffs_t), pointer :: cf => null()
     type(oscpar_t), pointer  :: op => null()
   contains
     private
     procedure, public :: init
     procedure, public :: eval
     procedure, public :: eval_logx
     procedure, public :: trans_matrix
  end type jacobian_rad_jcd_t

  ! Access specifiers

  private

  public :: jacobian_rad_jcd_t

  ! Procedures

contains

  subroutine init (this, cf, op)

    class(jacobian_rad_jcd_t), intent(out) :: this
    class(coeffs_t), intent(in), target    :: cf
    type(oscpar_t), intent(in), target     :: op

    ! Initialize the jacobian

    this%cf => cf
    this%op => op

    this%n_e = 2

    ! Finish

    return

  end subroutine init

!****

  subroutine eval (this, x, omega, A)

    class(jacobian_rad_jcd_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    complex(WP), intent(in)               :: omega
    complex(WP), intent(out)              :: A(:,:)
    
    ! Evaluate the Jacobian matrix

    call this%eval_logx(x, omega, A)

    A = A/x

    ! Finish

    return

  end subroutine eval

!****

  subroutine eval_logx (this, x, omega, A)

    class(jacobian_rad_jcd_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    complex(WP), intent(in)               :: omega
    complex(WP), intent(out)              :: A(:,:)
    
    $CHECK_BOUNDS(SIZE(A, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(A, 2),this%n_e)

    ! Evaluate the log(x)-space Jacobian matrix

    associate(V_g => this%cf%V(x)/this%cf%Gamma_1(x), U => this%cf%U(x), &
              As => this%cf%As(x), c_1 => this%cf%c_1(x), &
              omega_c => omega)

      A(1,1) = V_g - 1._WP
      A(1,2) = -V_g*c_1*omega_c**2
      
      A(2,1) = 1._WP - (As - U)/(c_1*omega_c**2)
      A(2,2) = As + 1._WP

    end associate

    ! Finish

    return

  end subroutine eval_logx

!****

  function trans_matrix (this, x, omega, to_canon)

    class(jacobian_rad_jcd_t), intent(in) :: this
    real(WP), intent(in)                  :: x
    complex(WP), intent(in)               :: omega
    logical, intent(in)                   :: to_canon
    complex(WP)                           :: trans_matrix(this%n_e,this%n_e)

    ! Calculate the transformation matrix to convert variables between the
    ! canonical formulation and the Dziembowski formulation

    if (to_canon) then

       associate(c_1 => this%cf%c_1(x), &
                 l => this%op%l)

         trans_matrix(1,1) = 1._WP
         trans_matrix(1,2) = 0._WP

         trans_matrix(2,1) = 0._WP
         trans_matrix(2,2) = c_1*omega**2

       end associate

    else

       associate(c_1 => this%cf%c_1(x), &
                 l => this%op%l)

         trans_matrix(1,1) = 1._WP
         trans_matrix(1,2) = 0._WP

         trans_matrix(2,1) = 0._WP
         trans_matrix(2,2) = 1._WP/(c_1*omega**2)

       end associate

    endif

    ! Finish

    return

  end function trans_matrix

end module gyre_jacobian_rad_jcd
