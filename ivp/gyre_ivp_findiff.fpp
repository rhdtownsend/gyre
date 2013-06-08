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

  public :: solve_findiff_GL2
  public :: solve_findiff_GL4
  public :: recon_findiff_GL2
  public :: recon_findiff_GL4
  public :: abscissa_findiff_GL2
  public :: abscissa_findiff_GL4

  ! Procedures

contains

  subroutine solve_findiff_GL2 (jc, omega, x_a, x_b, E_l, E_r, S)

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

    x = abscissa_findiff_GL2(x_a, x_b)
    dx = x_b - x_a

    call jc%eval(omega, x(1), A)

    ! Set up the solution matrices and scales

    E_l = 0.5_WP*dx*A + identity_matrix(jc%n_e)
    E_r = 0.5_WP*dx*A - identity_matrix(jc%n_e)

    S = ext_complex(1._WP)

    ! Finish

  end subroutine solve_findiff_GL2

!****

  subroutine solve_findiff_GL4 (jc, omega, x_a, x_b, E_l, E_r, S)

    class(jacobian_t), intent(in)    :: jc
    complex(WP), intent(in)          :: omega
    real(WP), intent(in)             :: x_a
    real(WP), intent(in)             :: x_b
    complex(WP), intent(out)         :: E_l(:,:)
    complex(WP), intent(out)         :: E_r(:,:)
    type(ext_complex_t), intent(out) :: S

    real(WP), parameter :: ALPHA_11 = 0.25_WP
    real(WP), parameter :: ALPHA_21 = 0.25_WP + SQRT(3._WP)/6._WP
    real(WP), parameter :: ALPHA_12 = 0.25_WP - SQRT(3._WP)/6._WP
    real(WP), parameter :: ALPHA_22 = 0.25_WP
    real(WP), parameter :: BETA_1 = 0.5_WP
    real(WP), parameter :: BETA_2 = 0.5_WP

    real(WP)    :: dx
    real(WP)    :: x(2)
    complex(WP) :: A(jc%n_e,jc%n_e,2)
    integer     :: n_e
    complex(WP) :: W(2*jc%n_e,2*jc%n_e)
    integer     :: i
    complex(WP) :: V(2*jc%n_e,jc%n_e)
    complex(WP) :: W_V(2*jc%n_e,jc%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),jc%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),jc%n_e)

    ! Solve the IVP across the interval x_a -> x_b using 4th-order
    ! Gauss-Legendre finite diffferences (see Ascher et al. 1995,
    ! Chap. 5)

    ! Evaluate the Jacobian

    x = abscissa_findiff_GL4(x_a, x_b)
    dx = x_b - x_a

    call jc%eval(omega, x(1), A(:,:,1))
    call jc%eval(omega, x(2), A(:,:,2))

    ! Set up the solution matrices and scales

    n_e = jc%n_e

    W(:n_e,:n_e) = -dx*ALPHA_11*A(:,:,1)
    W(n_e+1:,:n_e) = -dx*ALPHA_21*A(:,:,2)

    W(:n_e,n_e+1:) = -dx*ALPHA_12*A(:,:,1)
    W(n_e+1:,n_e+1:) = -dx*ALPHA_22*A(:,:,2)

    do i = 1,2*n_e
       W(i,i) = W(i,i) + 1._WP
    end do

    V(:n_e,:) = A(:,:,1)
    V(n_e+1:,:) = A(:,:,2)

    W_V = linear_solve(W, V)

    ! Set up the solution matrices and scales

    E_l = -dx*(BETA_1*W_V(:n_e,:) + BETA_2*W_V(n_e+1:,:))
    E_r = 0._WP

    do i = 1,n_e
       E_l(i,i) = E_l(i,i) - 1._WP
       E_r(i,i) = E_r(i,i) + 1._WP
    end do

    S = ext_complex(1._WP)

    ! Finish

  end subroutine solve_findiff_GL4

!****

  subroutine recon_findiff_GL2 (jc, omega, x_a, x_b, y_a, y_b, x, y)

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

  end subroutine recon_findiff_GL2

!****

  subroutine recon_findiff_GL4 (jc, omega, x_a, x_b, y_a, y_b, x, y)

    class(jacobian_t), intent(in) :: jc
    complex(WP), intent(in)       :: omega
    real(WP), intent(in)          :: x_a
    real(WP), intent(in)          :: x_b
    complex(WP), intent(in)       :: y_a(:)
    complex(WP), intent(in)       :: y_b(:)
    real(WP), intent(in)          :: x(:)
    complex(WP), intent(out)      :: y(:,:)

    $CHECK_BOUNDS(SIZE(y_a),jc%n_e)
    $CHECK_BOUNDS(SIZE(y_b),jc%n_e)
    
    $CHECK_BOUNDS(SIZE(y, 1),jc%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))

    call recon_findiff_GL2(jc, omega, x_a, x_b, y_a, y_b, x, y)

    ! Finish

    return

  end subroutine recon_findiff_GL4

!****

  function abscissa_findiff_GL2 (x_a, x_b) result (x)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP)             :: x(1)

    real(WP) :: dx

    ! Set up the abscissa for centered finite-differences

    dx = x_b - x_a

    x = x_a + [0.5_WP]*dx

    ! Finish

    return

  end function abscissa_findiff_GL2

!****

  function abscissa_findiff_GL4 (x_a, x_b) result (x)

    real(WP), intent(in) :: x_a
    real(WP), intent(in) :: x_b
    real(WP)             :: x(2)

    real(WP) :: dx

    ! Set up the abscissa for centered finite-differences

    dx = x_b - x_a

    x = x_a + (0.5_WP+[-1._WP,1._WP]*SQRT(3._WP)/6._WP)*dx

    ! Finish

    return

  end function abscissa_findiff_GL4

end module gyre_ivp_findiff
