! Module   : gyre_ivp_findiff
! Purpose  : solve IVPs across fixed intervals using finite differences
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
  public :: abscissa_findiff

  ! Procedures

contains

  subroutine solve_findiff (jc, omega, x_a, x_b, E_l, E_r, S)

    class(jacobian_t), intent(in)    :: jc
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    complex(WP), intent(out)         :: E_l(:,:)
    complex(WP), intent(out)         :: E_r(:,:)
    type(ext_complex_t), intent(out) :: S

    real(WP)    :: dx
    real(WP)    :: x(1)
    complex(WP) :: A(jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),jc%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),jc%n_e)

    ! Solve the IVP across the interval x_a -> x_b using centered
    ! (2nd-order) finite diffferences

    ! Evaluate the Jacobian

    x = abscissa_findiff(x_a, x_b)
    dx = x_b - x_a

    call jc%eval(omega, x(1), A)

    ! Set up the solution matrices and scales

    E_l = 0.5_WP*dx*A + identity_matrix(jc%n_e)
    E_r = 0.5_WP*dx*A - identity_matrix(jc%n_e)

    S = ext_complex(1._WP)

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

!****

  function abscissa_findiff (x_a, x_b) result (x)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP)             :: x(1)

    real(WP) :: dx

    ! Set up the abscissa for centered finite-differences

    dx = x_b - x_a

    x = x_a + [0.5_WP]*dx

    ! Finish

    return

  end function abscissa_findiff

end module gyre_ivp_findiff
