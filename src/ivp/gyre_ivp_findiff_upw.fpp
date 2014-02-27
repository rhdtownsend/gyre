! Module   : gyre_ivp_findiff_upw
! Purpose  : solve initial-value problems (finite differences w/ upwinding)
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

module gyre_ivp_findiff_upw

  ! Uses

  use core_kinds
  use core_linalg

  use gyre_jacobian
  use gyre_ivp_colloc
  use gyre_ext_arith
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (ivp_colloc_t) :: ivp_findiff_upw_t
     private
     integer, allocatable, public   :: stencil(:)
     class(jacobian_t), allocatable :: jc
   contains
     private
     procedure, public :: solve => solve_
     procedure, public :: recon => recon_
     procedure, public :: abscissa => abscissa_
  end type ivp_findiff_upw_t

  ! Interfaces

  interface ivp_findiff_upw_t
     module procedure ivp_findiff_upw_t_
  end interface ivp_findiff_upw_t

  ! Access specifiers

  private

  public :: ivp_findiff_upw_t

contains

  function ivp_findiff_upw_t_ (jc) result (iv)

    class(jacobian_t), intent(in) :: jc
    type(ivp_findiff_upw_t)       :: iv

    ! Construct the ivp_findiff_upw_t

    allocate(iv%jc, SOURCE=jc)

    allocate(iv%stencil(jc%n_e))

    iv%n_e = jc%n_e

    ! Finish

    return
    
  end function ivp_findiff_upw_t_

!****

  subroutine solve_ (this, omega, x_a, x_b, E_l, E_r, S, use_real)

    class(ivp_findiff_upw_t), intent(in) :: this
    complex(WP), intent(in)             :: omega
    real(WP), intent(in)                :: x_a
    real(WP), intent(in)                :: x_b
    complex(WP), intent(out)            :: E_l(:,:)
    complex(WP), intent(out)            :: E_r(:,:)
    type(ext_complex_t), intent(out)    :: S
    logical, optional, intent(in)       :: use_real

    real(WP)    :: dx
    real(WP)    :: x(3)
    complex(WP) :: A(this%n_e,this%n_e,3)
    integer     :: i

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
    call this%jc%eval(x(3), omega, A(:,:,3))

    ! Set up the solution matrices and scales

    E_l =  identity_matrix(this%n_e)
    E_r = -identity_matrix(this%n_e)

    do i = 1,this%n_e

       select case (this%stencil(i))
       case (-1)
          E_l(i,:) = E_l(i,:) + dx*A(i,:,1)
       case (0)
          E_l(i,:) = E_l(i,:) + 0.5_WP*dx*A(i,:,2)
          E_r(i,:) = E_r(i,:) + 0.5_WP*dx*A(i,:,2)
       case (1)
          E_r(i,:) = E_r(i,:) + dx*A(i,:,3)
       case default 
          $ABORT(Invalid stencil)
       end select

    end do

    S = ext_complex_t(1._WP)

    ! Finish

  end subroutine solve_

!****

  subroutine recon_ (this, omega, x_a, x_b, y_a, y_b, x, y, use_real)

    class(ivp_findiff_upw_t), intent(in) :: this
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

    recon_loop : do i = 1,SIZE(x)

       w = (x(i) - x_a)/(x_b - x_a)

       y(:,i) = y_a*(1._WP-w) + y_b*w

    end do recon_loop

    ! Finish

    return

  end subroutine recon_

!****

  function abscissa_ (this, x_a, x_b) result (x)

    class(ivp_findiff_upw_t), intent(in) :: this
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    real(WP), allocatable                :: x(:)

    ! Set up the abscissa

    x = [x_a,0.5*(x_a+x_b),x_b]

    ! Finish

    return

  end function abscissa_

end module gyre_ivp_findiff_upw
