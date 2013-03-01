! Module   : gyre_ad_jacobian
! Purpose  : Adiabatic Jacobian evaluation

$include 'core.inc'

module gyre_ad_jacobian

  ! Uses

  use core_kinds

  use gyre_jacobian
  use gyre_mech_coeffs

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends(jacobian_t) :: ad_jacobian_t
     private
     class(mech_coeffs_t), pointer :: mc => null()
     real(WP)                      :: lambda_0
     integer                       :: l
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

  subroutine init (this, mc, lambda_0, l)

    class(ad_jacobian_t), intent(out)        :: this
    class(mech_coeffs_t), intent(in), target :: mc
    real(WP), intent(in)                     :: lambda_0
    integer, intent(in)                      :: l

    ! Initialize the ad_jacobian

    this%mc => mc

    this%lambda_0 = lambda_0
    this%l = l

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
              lambda_0 => this%lambda_0, l => this%l)

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
