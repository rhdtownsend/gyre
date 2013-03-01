! Module   : gyre_ivp_findiff
! Purpose  : solve IVPs across fixed intervals using finite differences

$include 'core.inc'

module gyre_ivp_findiff

  ! Uses

  use core_kinds

  use gyre_jacobian
  use gyre_ext_arith
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: solve_findiff
  public :: recon_findiff

  ! Procedures

contains

  subroutine solve_findiff (jc, omega, x_a, x_b, E_l, E_r, scale)

    class(jacobian_t), intent(in)    :: jc
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    complex(WP), intent(out)         :: E_l(:,:)
    complex(WP), intent(out)         :: E_r(:,:)
    type(ext_complex_t), intent(out) :: scale

    real(WP)    :: x
    real(WP)    :: dx
    complex(WP) :: A(jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),jc%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),jc%n_e)

    ! Solve the IVP across the interval x_a -> x_b using centered
    ! (2nd-order) finite diffferences

    ! Evaluate the Jacobian

    x = 0.5_WP*(x_a + x_b)
    dx = x_b - x_a

    call jc%eval(omega, x, A)

    ! Set up the solution matrices and scales

    E_l = 0.5_WP*dx*A + identity_matrix(jc%n_e)
    E_r = 0.5_WP*dx*A - identity_matrix(jc%n_e)

    scale = ext_complex(1._WP)

    ! Finish

  end subroutine solve_findiff

!****

  subroutine recon_findiff (jc, omega, x_a, x_b, y_a, y_b, x, y)

    class(jacobian_t), intent(in) :: jc
    complex(WP), intent(in)       :: omega
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(in)       :: y_a(:)
    complex(WP), intent(in)       :: y_b(:)
    real(WP), intent(in)          :: x(:)
    complex(WP), intent(out)      :: y(:,:)

    integer  :: i
    real(WP) :: w

    $CHECK_BOUNDS(SIZE(y_a),jc%n_e)
    $CHECK_BOUNDS(SIZE(y_b),jc%n_e)
    
    $CHECK_BOUNDS(SIZE(y, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))
    
    ! Reconstruct the solution within the interval x_a -> x_b using
    ! linear interpolation

    recon_loop : do i = 1,SIZE(x)

       w = (x(i) - x_a)/(x_b - x_a)

       y(:,i) = y_a*(1._WP-w) + y_b*w

    end do recon_loop

    ! Finish

    return

  end subroutine recon_findiff

end module gyre_ivp_findiff
