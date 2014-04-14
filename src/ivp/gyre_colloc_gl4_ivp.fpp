! Module   : gyre_colloc_gl4_ivp
! Purpose  : solve initial-value problems (4th-order Gauss-Legendre collocation)
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

module gyre_colloc_gl4_ivp

  ! Uses

  use core_kinds
  use core_linalg

  use gyre_jacobian
  use gyre_ivp_colloc
  use gyre_ext_arith

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (ivp_colloc_t) :: colloc_gl4_ivp_t
     private
     class(jacobian_t), allocatable :: jc
   contains
     private
     procedure, public :: solve => solve_
     procedure, public :: recon => recon_
     procedure, public :: abscissa => abscissa_
  end type colloc_gl4_ivp_t

  ! Interfaces

  interface colloc_gl4_ivp_t
     module procedure colloc_gl4_ivp_t_
  end interface colloc_gl4_ivp_t

  ! Access specifiers

  private

  public :: colloc_gl4_ivp_t

contains

  function colloc_gl4_ivp_t_ (jc) result (iv)

    class(jacobian_t), intent(in) :: jc
    type(colloc_gl4_ivp_t)        :: iv
    
    ! Construct the colloc_gl4_ivp_t

    allocate(iv%jc, SOURCE=jc)

    iv%n_e = jc%n_e

    ! Finish

    return
    
  end function colloc_gl4_ivp_t_

!****

  subroutine solve_ (this, omega, x_a, x_b, E_l, E_r, S, use_real)

    class(colloc_gl4_ivp_t), intent(in) :: this
    complex(WP), intent(in)             :: omega
    real(WP), intent(in)                :: x_a
    real(WP), intent(in)                :: x_b
    complex(WP), intent(out)            :: E_l(:,:)
    complex(WP), intent(out)            :: E_r(:,:)
    type(ext_complex_t), intent(out)    :: S
    logical, optional, intent(in)       :: use_real

    real(WP), parameter :: ALPHA_11 = 0.25_WP
    real(WP), parameter :: ALPHA_21 = 0.25_WP + SQRT(3._WP)/6._WP
    real(WP), parameter :: ALPHA_12 = 0.25_WP - SQRT(3._WP)/6._WP
    real(WP), parameter :: ALPHA_22 = 0.25_WP
    real(WP), parameter :: BETA_1 = 0.5_WP
    real(WP), parameter :: BETA_2 = 0.5_WP

    real(WP)    :: dx
    real(WP)    :: x(2)
    complex(WP) :: A(this%n_e,this%n_e,2)
    integer     :: n_e
    complex(WP) :: W(2*this%n_e,2*this%n_e)
    integer     :: i
    complex(WP) :: V(2*this%n_e,this%n_e)
    integer     :: ipiv(2*this%n_e)
    integer     :: info

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Solve the IVP across the interval x_a -> x_b

    ! Evaluate the Jacobians

    x = this%abscissa(x_a, x_b)
    dx = x_b - x_a

    call this%jc%eval(x(1), omega, A(:,:,1))
    call this%jc%eval(x(2), omega, A(:,:,2))

    ! Set up the solution matrices and scales

    n_e = this%n_e

    W(:n_e,:n_e) = -dx*ALPHA_11*A(:,:,1)
    W(n_e+1:,:n_e) = -dx*ALPHA_21*A(:,:,2)

    W(:n_e,n_e+1:) = -dx*ALPHA_12*A(:,:,1)
    W(n_e+1:,n_e+1:) = -dx*ALPHA_22*A(:,:,2)

    do i = 1,2*n_e
       W(i,i) = W(i,i) + 1._WP
    end do

    V(:n_e,:) = A(:,:,1)
    V(n_e+1:,:) = A(:,:,2)

    call XGESV(2*n_e, n_e, W, 2*n_e, ipiv, V, 2*n_e, info)
    $ASSERT(info == 0,Non-zero return from XGESV)
    
    ! Set up the solution matrices and scales

    E_l = -dx*(BETA_1*V(:n_e,:) + BETA_2*V(n_e+1:,:))
    E_r = 0._WP

    do i = 1,n_e
       E_l(i,i) = E_l(i,i) - 1._WP
       E_r(i,i) = E_r(i,i) + 1._WP
    end do

    S = ext_complex_t(1._WP)

    ! Finish

  end subroutine solve_

!****

  subroutine recon_ (this, omega, x_a, x_b, y_a, y_b, x, y, use_real)

    class(colloc_gl4_ivp_t), intent(in) :: this
    complex(WP), intent(in)             :: omega
    real(WP), intent(in)                :: x_a
    real(WP), intent(in)                :: x_b
    complex(WP), intent(in)             :: y_a(:)
    complex(WP), intent(in)             :: y_b(:)
    real(WP), intent(in)                :: x(:)
    complex(WP), intent(out)            :: y(:,:)
    logical, optional, intent(in)       :: use_real

    integer  :: i
    real(WP) :: w

    $CHECK_BOUNDS(SIZE(y_a),this%n_e)
    $CHECK_BOUNDS(SIZE(y_b),this%n_e)
    
    $CHECK_BOUNDS(SIZE(y, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))
    
    ! Reconstruct the solution within the interval x_a -> x_b

    ! (This just uses straight linear interpolation, and needs to be updated)

    recon_loop : do i = 1,SIZE(x)

       w = (x(i) - x_a)/(x_b - x_a)

       y(:,i) = y_a*(1._WP-w) + y_b*w

    end do recon_loop

    ! Finish

    return

  end subroutine recon_

!****

  function abscissa_ (this, x_a, x_b) result (x)

    class(colloc_gl4_ivp_t), intent(in) :: this
    real(WP), intent(in)                :: x_a
    real(WP), intent(in)                :: x_b
    real(WP), allocatable               :: x(:)

    real(WP) :: dx

    ! Set up the abscissa

    dx = x_b - x_a

    x = x_a + (0.5_WP+[-1._WP,1._WP]*SQRT(3._WP)/6._WP)*dx

    ! Finish

    return

  end function abscissa_

end module gyre_colloc_gl4_ivp
