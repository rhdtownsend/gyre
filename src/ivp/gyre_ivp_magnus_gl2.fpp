! Module   : gyre_ivp_magnus_GL2
! Purpose  : solve initial-value problems (2nd-order Gauss-Legendre Magnus)
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

module gyre_ivp_magnus_GL2

  ! Uses

  use core_kinds
  use core_constants

  use gyre_jacobian
  use gyre_ivp_magnus

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (ivp_magnus_t) :: ivp_magnus_GL2_t
     private
     class(jacobian_t), pointer :: jc => null()
   contains
     private
     procedure, public :: init
     procedure, public :: eval_dOmega
     procedure, public :: abscissa
  end type ivp_magnus_GL2_t

  ! Access specifiers

  private

  public :: ivp_magnus_GL2_t

  ! Procedures

contains

  subroutine init (this, jc)

    class(ivp_magnus_GL2_t), intent(out)  :: this
    class(jacobian_t), intent(in), target :: jc

    ! Initialize the ivp_t

    this%jc => jc

    this%n_e = jc%n_e

    ! Finish

    return
    
  end subroutine init

!****

  subroutine eval_dOmega (this, omega, x_a, x_b, dOmega)

    class(ivp_magnus_GL2_t), intent(in) :: this
    complex(WP), intent(in)             :: omega
    real(WP), intent(in)                :: x_a
    real(WP), intent(in)                :: x_b
    complex(WP), intent(out)            :: dOmega(:,:)

    real(WP)    :: dx
    real(WP)    :: x(1)
    complex(WP) :: A(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(dOmega, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(dOmega, 2),this%n_e)

    ! Evaluate the Magnus slope matrix

    ! Calculate the Jacobian

    dx = x_b - x_a
    x = this%abscissa(x_a, x_b)

    call this%jc%eval(x(1), omega, A)

    ! Set up the slope matrix

    dOmega = A

    ! Finish

    return

  end subroutine eval_dOmega

!****

  function abscissa (this, x_a, x_b) result (x)

    class(ivp_magnus_GL2_t), intent(in) :: this
    real(WP), intent(in)                :: x_a
    real(WP), intent(in)                :: x_b
    real(WP), allocatable               :: x(:)

    real(WP) :: dx

    ! Set up the abscissa

    dx = x_b - x_a

    x = x_a + [0.5_WP]*dx

    ! Finish

    return
    
  end function abscissa

end module gyre_ivp_magnus_GL2
