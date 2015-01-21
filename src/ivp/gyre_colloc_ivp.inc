! Incfile  : gyre_colloc_ivp
! Purpose  : initial-value solvers (collocation, template)
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

module gyre_${T}_colloc_ivp

  ! Uses

  use core_kinds
  use core_linalg

  use gyre_eqns
  use gyre_ext
  use gyre_ivp
  use gyre_linalg

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: GL2_SCHEME = 1
  integer, parameter :: GL4_SCHEME = 2

  ! Derived-type definitions

  type, extends (${T}_ivp_t) :: ${T}_colloc_ivp_t
     private
     class(${T}_eqns_t), allocatable :: eq
     integer                         :: scheme
   contains
     private
     procedure, public :: shoot => shoot_
     procedure, public :: recon => recon_
     procedure, public :: abscissa => abscissa_
  end type ${T}_colloc_ivp_t

  ! Interfaces

  interface ${T}_colloc_ivp_t
     module procedure ${T}_colloc_ivp_t_
  end interface ${T}_colloc_ivp_t

  ! Access specifiers

  private

  public :: ${T}_colloc_ivp_t

contains

  function ${T}_colloc_ivp_t_ (eq, scheme) result (iv)

    class(${T}_eqns_t), intent(in) :: eq
    character(*), intent(in)       :: scheme
    type(${T}_colloc_ivp_t)        :: iv

    ! Construct the colloc_ivp_t

    allocate(iv%eq, SOURCE=eq)

    select case (scheme)
    case ('GL2')
       iv%scheme = GL2_SCHEME
    case ('GL4')
       iv%scheme = GL4_SCHEME
    case default
       $ABORT(Invalid scheme)
    end select
       
    iv%n_e = eq%n_e

    ! Finish

    return
    
  end function ${T}_colloc_ivp_t_

!****

  subroutine shoot_ (this, omega, x_a, x_b, E_l, E_r, S)

    class(${T}_colloc_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP), intent(out)               :: E_l(:,:)
    $TYPE(WP), intent(out)               :: E_r(:,:)
    type(${T}_ext_t), intent(out)        :: S

    ! Set up the shooting matrices and scales

    select case (this%scheme)
    case (GL2_SCHEME)
       call shoot_GL2_(this, omega, x_a, x_b, E_l, E_r, S)
    case (GL4_SCHEME)
       call shoot_GL4_(this, omega, x_a, x_b, E_l, E_r, S)
    case default
       $ABORT(Invalid scheme)
    end select

    ! Finish

    return

  end subroutine shoot_

!****

  subroutine shoot_GL2_ (this, omega, x_a, x_b, E_l, E_r, S)

    class(${T}_colloc_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP), intent(out)               :: E_l(:,:)
    $TYPE(WP), intent(out)               :: E_r(:,:)
    type(${T}_ext_t), intent(out)        :: S

    real(WP)  :: dx
    real(WP)  :: x(1)
    $TYPE(WP) :: A(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Evaluate the RHS matrix

    x = this%abscissa(x_a, x_b)
    dx = x_b - x_a

    A = this%eq%A(x(1), omega)

    ! Set up the GL2 shooting matrices and scales

    E_l = 0.5_WP*dx*A + identity_matrix(this%n_e)
    E_r = 0.5_WP*dx*A - identity_matrix(this%n_e)

    S = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine shoot_GL2_

!****

  subroutine shoot_GL4_ (this, omega, x_a, x_b, E_l, E_r, S)

    class(${T}_colloc_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP), intent(out)               :: E_l(:,:)
    $TYPE(WP), intent(out)               :: E_r(:,:)
    type(${T}_ext_t), intent(out)        :: S

    real(WP), parameter :: ALPHA_11 = 0.25_WP
    real(WP), parameter :: ALPHA_21 = 0.25_WP + SQRT(3._WP)/6._WP
    real(WP), parameter :: ALPHA_12 = 0.25_WP - SQRT(3._WP)/6._WP
    real(WP), parameter :: ALPHA_22 = 0.25_WP
    real(WP), parameter :: BETA_1 = 0.5_WP
    real(WP), parameter :: BETA_2 = 0.5_WP

    real(WP)  :: dx
    real(WP)  :: x(2)
    $TYPE(WP) :: A(this%n_e,this%n_e,2)
    integer   :: n_e
    $TYPE(WP) :: W(2*this%n_e,2*this%n_e)
    integer   :: i
    $TYPE(WP) :: V(2*this%n_e,this%n_e)
    integer   :: ipiv(2*this%n_e)
    integer   :: info

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Evaluate the RHS matrices

    x = this%abscissa(x_a, x_b)
    dx = x_b - x_a

    A(:,:,1) = this%eq%A(x(1), omega)
    A(:,:,2) = this%eq%A(x(2), omega)

    ! Set up the GL4 shooting matrices and scales

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

    S = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine shoot_GL4_

!****

  subroutine recon_ (this, omega, x_a, x_b, y_a, y_b, x, y)

    class(${T}_colloc_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP), intent(in)                :: y_a(:)
    $TYPE(WP), intent(in)                :: y_b(:)
    real(WP), intent(in)                 :: x(:)
    $TYPE(WP), intent(out)               :: y(:,:)

    ! Reconstruct the solution 

    select case (this%scheme)
    case (GL2_SCHEME)
       call recon_GL2_ (this, omega, x_a, x_b, y_a, y_b, x, y)
    case (GL4_SCHEME)
       call recon_GL4_ (this, omega, x_a, x_b, y_a, y_b, x, y)
    case default
       $ABORT(Invalid scheme)
    end select

    ! Finish

    return

  end subroutine recon_

!****

  subroutine recon_GL2_ (this, omega, x_a, x_b, y_a, y_b, x, y)

    class(${T}_colloc_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP), intent(in)                :: y_a(:)
    $TYPE(WP), intent(in)                :: y_b(:)
    real(WP), intent(in)                 :: x(:)
    $TYPE(WP), intent(out)               :: y(:,:)

    integer  :: i
    real(WP) :: w

    $CHECK_BOUNDS(SIZE(y_a),this%n_e)
    $CHECK_BOUNDS(SIZE(y_b),this%n_e)
    
    $CHECK_BOUNDS(SIZE(y, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))
    
    ! Reconstruct the GL2 solution

    recon_loop : do i = 1,SIZE(x)

       w = (x(i) - x_a)/(x_b - x_a)

       y(:,i) = y_a*(1._WP-w) + y_b*w

    end do recon_loop

    ! Finish

    return

  end subroutine recon_GL2_

!****

  subroutine recon_GL4_ (this, omega, x_a, x_b, y_a, y_b, x, y)

    class(${T}_colloc_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP), intent(in)                :: y_a(:)
    $TYPE(WP), intent(in)                :: y_b(:)
    real(WP), intent(in)                 :: x(:)
    $TYPE(WP), intent(out)               :: y(:,:)

    integer  :: i
    real(WP) :: w

    $CHECK_BOUNDS(SIZE(y_a),this%n_e)
    $CHECK_BOUNDS(SIZE(y_b),this%n_e)
    
    $CHECK_BOUNDS(SIZE(y, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(y, 2),SIZE(x))
    
    ! Reconstruct the GL4 solution

    ! (This just uses straight linear interpolation, and needs to be updated)

    recon_loop : do i = 1,SIZE(x)

       w = (x(i) - x_a)/(x_b - x_a)

       y(:,i) = y_a*(1._WP-w) + y_b*w

    end do recon_loop

    ! Finish

    return

  end subroutine recon_GL4_

!****

  function abscissa_ (this, x_a, x_b) result (x)

    class(${T}_colloc_ivp_t), intent(in) :: this
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    real(WP), allocatable                :: x(:)

    ! Set up the abscissa

    select case (this%scheme)
    case (GL2_SCHEME)
       x = abscissa_GL2_(this, x_a, x_b)
    case (GL4_SCHEME)
       x = abscissa_GL4_(this, x_a, x_b)
    case default
       $ABORT(Invalid scheme)
    end select

    return

  end function abscissa_

!****

  function abscissa_GL2_ (this, x_a, x_b) result (x)

    class(${T}_colloc_ivp_t), intent(in) :: this
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    real(WP), allocatable                :: x(:)

    real(WP) :: dx

    ! Set up the GL2 abscissa

    dx = x_b - x_a

    x = x_a + [0.5_WP]*dx

    ! Finish

    return

  end function abscissa_GL2_

!****

  function abscissa_GL4_ (this, x_a, x_b) result (x)

    class(${T}_colloc_ivp_t), intent(in) :: this
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    real(WP), allocatable                :: x(:)

    real(WP) :: dx

    ! Set up the GL4 abscissa

    dx = x_b - x_a

    x = x_a + (0.5_WP+[-1._WP,1._WP]*SQRT(3._WP)/6._WP)*dx

    ! Finish

    return

  end function abscissa_GL4_

end module gyre_${T}_colloc_ivp

$endsub

