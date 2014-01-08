! Module   : gyre_ivp_magnus_GL4
! Purpose  : solve initial-value problems (4th-order Gauss-Legendre Magnus)
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

module gyre_ivp_magnus_GL4

  ! Uses

  use core_kinds
  use core_constants

  use gyre_jacobian
  use gyre_linalg
  use gyre_ivp_magnus

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (ivp_magnus_t) :: ivp_magnus_GL4_t
     private
     class(jacobian_t), allocatable :: jc
   contains
     private
     procedure, public :: eval_dOmega
     procedure, public :: abscissa
  end type ivp_magnus_GL4_t

  ! Interfaces

  interface ivp_magnus_GL4_t
     module procedure init_iv
  end interface ivp_magnus_GL4_t

  ! Access specifiers

  private

  public :: ivp_magnus_GL4_t

  ! Procedures

contains

  function init_iv (jc) result (iv)

    class(jacobian_t), intent(in) :: jc
    type(ivp_magnus_GL4_t)        :: iv

    ! Construct the ivp_magnus_GL4

    allocate(iv%jc, SOURCE=jc)

    iv%n_e = jc%n_e

    ! Finish

    return
    
  end function init_iv

!****

  subroutine eval_dOmega (this, omega, x_a, x_b, dOmega)

    class(ivp_magnus_GL4_t), intent(in) :: this
    complex(WP), intent(in)             :: omega
    real(WP), intent(in)                :: x_a
    real(WP), intent(in)                :: x_b
    complex(WP), intent(out)            :: dOmega(:,:)

    real(WP)    :: dx
    real(WP)    :: x(2)
    complex(WP) :: A(this%n_e,this%n_e,2)
    complex(WP) :: dalpha(this%n_e,this%n_e,2)

    $CHECK_BOUNDS(SIZE(dOmega, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(dOmega, 2),this%n_e)

    ! Evaluate the Magnus slope matrix

    ! Calculate the Jacobians

    dx = x_b - x_a
    x = this%abscissa(x_a, x_b)

    call this%jc%eval(x(1), omega, A(:,:,1))
    call this%jc%eval(x(2), omega, A(:,:,2))

    ! Set up the Magnus slope matrix (Blanes et al. 2009, eqns. 243
    ! and 253; note that the 12 in the denominator of their expression
    ! for alpha_2 is erroneous)

    dalpha(:,:,1) = 0.5_WP*(A(:,:,1) + A(:,:,2))
    dalpha(:,:,2) = SQRT(3._WP)*(A(:,:,2) - A(:,:,1))

    dOmega = dalpha(:,:,1) - dx*commutator(dalpha(:,:,1), dalpha(:,:,2))/12._WP
    
    ! Finish

    return

  end subroutine eval_dOmega

!****

  function abscissa (this, x_a, x_b) result (x)

    class(ivp_magnus_GL4_t), intent(in) :: this
    real(WP), intent(in)                :: x_a
    real(WP), intent(in)                :: x_b
    real(WP), allocatable               :: x(:)

    real(WP) :: dx

    ! Set up the abscissa

    dx = x_b - x_a

    x = x_a + (0.5_WP+[-1._WP,1._WP]*SQRT(3._WP)/6._WP)*dx

    ! Finish

    return

  end function abscissa

end module gyre_ivp_magnus_GL4
