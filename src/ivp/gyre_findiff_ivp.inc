! Incfile  : gyre_findiff_ivp
! Purpose  : initial-value solvers (finite difference, template)
!
! Copyright 2013-2014 Rich Townsend
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_findiff_ivp

  ! Uses

  use core_kinds
  use core_linalg

  use gyre_eqns
  use gyre_ivp
  use gyre_ext
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (${T}_ivp_t) :: ${T}_findiff_ivp_t
     private
     class(${T}_eqns_t), allocatable :: eq
   contains
     private
     procedure, public :: shoot => shoot_
     procedure, public :: shoot_w => shoot_w_
     procedure, public :: recon => recon_
  end type ${T}_findiff_ivp_t

  ! Interfaces

  interface ${T}_findiff_ivp_t
     module procedure ${T}_findiff_ivp_t_
  end interface ${T}_findiff_ivp_t

  ! Access specifiers

  private

  public :: ${T}_findiff_ivp_t

contains

  function ${T}_findiff_ivp_t_ (eq) result (iv)

    class(${T}_eqns_t), intent(in) :: eq
    type(${T}_findiff_ivp_t)       :: iv

    ! Construct the findiff_ivp_t

    allocate(iv%eq, SOURCE=eq)

    iv%n_e = eq%n_e

    ! Finish

    return
    
  end function ${T}_findiff_ivp_t_

!****

  subroutine shoot_ (this, omega, x_a, x_b, E_l, E_r, S)

    class(${T}_findiff_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                 :: omega
    real(WP), intent(in)                  :: x_a
    real(WP), intent(in)                  :: x_b
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: S

    real(WP) :: w(this%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Set up the standard-weight shooting matrices and scales

    w = 0.5_WP

    call this%shoot_w(w, omega, x_a, x_b, E_l, E_r, S)

    ! Finish

  end subroutine shoot_

!****

  subroutine shoot_w_ (this, w, omega, x_a, x_b, E_l, E_r, S)

    class(${T}_findiff_ivp_t), intent(in) :: this
    real(WP), intent(in)                  :: w(:)
    $TYPE(WP), intent(in)                 :: omega
    real(WP), intent(in)                  :: x_a
    real(WP), intent(in)                  :: x_b
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: S

    real(WP)  :: dx
    real(WP)  :: x(2)
    $TYPE(WP) :: xA(this%n_e,this%n_e,2)
    real(WP)  :: x_m

    $CHECK_BOUNDS(SIZE(w),this%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Evaluate the RHS matrix

    dx = x_b - x_a

    x = [x_a,x_b]

    xA(:,:,1) = this%eq%xA(x(1), omega)
    xA(:,:,2) = this%eq%xA(x(2), omega)

    ! Set up the shooting matrices and scales

    x_m = 0.5_WP*(x_a + x_b)

    E_l = -x_m*identity_matrix(this%n_e) - dx*MATMUL(diagonal_matrix(1._WP-w), xA(:,:,1))
    E_r =  x_m*identity_matrix(this%n_e) - dx*MATMUL(diagonal_matrix(      w), xA(:,:,2))

    S = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine shoot_w_
  
!****

  subroutine recon_ (this, omega, x_a, x_b, y_a, y_b, x, y)

    class(${T}_findiff_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                 :: omega
    real(WP), intent(in)                  :: x_a
    real(WP), intent(in)                  :: x_b
    $TYPE(WP), intent(in)                 :: y_a(:)
    $TYPE(WP), intent(in)                 :: y_b(:)
    real(WP), intent(in)                  :: x(:)
    $TYPE(WP), intent(out)                :: y(:,:)

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

end module gyre_${T}_findiff_ivp

$endsub

