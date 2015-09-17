! Incfile  : gyre_modcolloc_ivp
! Purpose  : initial-value solvers (collocation, template)
!
! Copyright 2013-2015 Rich Townsend
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
  integer, parameter :: GL6_SCHEME = 3

  ! Derived-type definitions

  type, extends (${T}_ivp_t) :: ${T}_colloc_ivp_t
     private
     class(${T}_eqns_t), allocatable :: eq
     integer                         :: scheme
   contains
     private
     procedure, public :: shoot => shoot_
     procedure, public :: recon => recon_
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
    case ('GL6')
       iv%scheme = GL6_SCHEME
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
    case (GL6_SCHEME)
       call shoot_GL6_(this, omega, x_a, x_b, E_l, E_r, S)
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
    $TYPE(WP) :: A(this%n_e,this%n_e,1)

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Evaluate the RHS matrix

    dx = x_b - x_a

    x = x_a + [0.5_WP]*dx

    A(:,:,1) = this%eq%A(x(1), omega)

    ! Set up the GL2 shooting matrices and scales

    E_l = 0.5_WP*dx*A(:,:,1) + identity_matrix(this%n_e)
    E_r = 0.5_WP*dx*A(:,:,1) - identity_matrix(this%n_e)

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

    real(WP), parameter :: A_IJ(2,2) = RESHAPE([0.25_WP,&
                                                0.25_WP+SQRT(3._WP)/6._WP, &
                                                0.25_WP-SQRT(3._WP)/6._WP, &
                                                0.25_WP], [2,2])
    real(WP), parameter :: B_I(2) = [0.5_WP, &
                                     0.5_WP]
    real(WP), parameter :: C_I(2) = [0.5_WP-SQRT(3._WP)/6._WP, &
                                     0.5_WP+SQRT(3._WP)/6._WP]

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)
 
    ! Set up the GL4 shooting matrices and scales

    if (x_a < 0.5_WP) then
       call shoot_irk_(this, omega, 2, x_a, x_b, A_IJ, B_I, C_I, E_l, E_r)
    else
       call shoot_irk_(this, omega, 2, x_b, x_a, A_IJ, B_I, C_I, E_r, E_l)
    endif

    S = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine shoot_GL4_

!****

  subroutine shoot_GL6_ (this, omega, x_a, x_b, E_l, E_r, S)

    class(${T}_colloc_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    $TYPE(WP), intent(out)               :: E_l(:,:)
    $TYPE(WP), intent(out)               :: E_r(:,:)
    type(${T}_ext_t), intent(out)        :: S

    real(WP), parameter :: A_IJ(3,3) = RESHAPE([5._WP/36._WP, &
                                                5._WP/36._WP+SQRT(15._WP)/24._WP, &
                                                5._WP/36._WP+SQRT(15._WP)/30._WP, &
                                                2._WP/9._WP-SQRT(15._WP)/15._WP, &
                                                2._WP/9._WP, &
                                                2._WP/9._WP+SQRT(15._WP)/15._WP, &
                                                5._WP/36._WP-SQRT(15._WP)/30._WP, &
                                                5._WP/36._WP-SQRT(15._WP)/24._WP, &
                                                5._WP/36._WP], [3,3])
    real(WP), parameter :: B_I(3) = [5._WP/18._WP, &
                                     4._WP/9._WP, &
                                     5._WP/18._WP]
    real(WP), parameter :: C_I(3) = [0.5_WP-SQRT(15._WP)/10._WP, &
                                     0.5_WP, &
                                     0.5_WP+SQRT(15._WP)/10._WP]

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Set up the GL6 shooting matrices and scales

    if (x_a < 0.5_WP) then
       call shoot_irk_(this, omega, 3, x_a, x_b, A_IJ, B_I, C_I, E_l, E_r)
    else
       call shoot_irk_(this, omega, 3, x_b, x_a, A_IJ, B_I, C_I, E_r, E_l)
    endif

    S = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine shoot_GL6_

!****

  subroutine shoot_irk_ (this, omega, n_s, x_a, x_b, a_ij, b_i, c_i, E_a, E_b)

    class(${T}_colloc_ivp_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    integer, intent(in)                  :: n_s
    real(WP), intent(in)                 :: x_a
    real(WP), intent(in)                 :: x_b
    real(WP), intent(in)                 :: a_ij(:,:)
    real(WP), intent(in)                 :: b_i(:)
    real(WP), intent(in)                 :: c_i(:)
    $TYPE(WP), intent(out)               :: E_a(:,:)
    $TYPE(WP), intent(out)               :: E_b(:,:)

    $TYPE(WP) :: K(this%n_e,this%n_e,n_s)
    real(WP)  :: dx
    integer   :: i
    
    $CHECK_BOUNDS(SIZE(a_ij, 1),n_s)
    $CHECK_BOUNDS(SIZE(a_ij, 2),n_s)

    $CHECK_BOUNDS(SIZE(b_i),n_s)

    $CHECK_BOUNDS(SIZE(c_i),n_s)

    $CHECK_BOUNDS(SIZE(E_a, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_a, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_b, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_b, 2),this%n_e)

    ! Shoot across the interval [x_a,x_b] using an n_s-step implicit
    ! Runge Kutta (IRK) scheme

    ! Calculate the intermediate-state matrices

    K = irk_K_(this%eq, omega, n_s, x_a, x_b, a_ij, c_i)

    ! Set up the shooting matrices

    dx = x_b - x_a

    E_a = identity_matrix(this%n_e)

    do i = 1, n_s
       E_a = E_a + dx*b_i(i)*K(:,:,i)
    enddo

    E_b = -identity_matrix(this%n_e)

    ! Finish

    return

  end subroutine shoot_irk_

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
       call recon_GL2_(this, omega, x_a, x_b, y_a, y_b, x, y)
    case (GL4_SCHEME)
       call recon_GL2_(this, omega, x_a, x_b, y_a, y_b, x, y)
    case (GL6_SCHEME)
       call recon_GL2_(this, omega, x_a, x_b, y_a, y_b, x, y)
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

  function irk_K_ (eq, omega, n_s, x_a, x_b, a_ij, c_i) result (K)

    class(${T}_eqns_t), intent(in) :: eq
    $TYPE(WP), intent(in)          :: omega
    integer, intent(in)            :: n_s
    real(WP), intent(in)           :: x_a
    real(WP), intent(in)           :: x_b
    real(WP), intent(in)           :: a_ij(:,:)
    real(WP), intent(in)           :: c_i(:)
    $TYPE(WP)                      :: K(eq%n_e,eq%n_e,n_s)

    real(WP)  :: dx
    $TYPE(WP) :: A(eq%n_e,eq%n_e,n_s)
    $TYPE(WP) :: M(n_s*eq%n_e,n_s*eq%n_e)
    integer   :: i
    integer   :: i_a
    integer   :: i_b
    integer   :: j
    integer   :: j_a
    integer   :: j_b
    $TYPE(WP) :: B(n_s*eq%n_e,eq%n_e)
    $TYPE(WP) :: X(n_s*eq%n_e,eq%n_e)

    $CHECK_BOUNDS(SIZE(a_ij, 1),n_s)
    $CHECK_BOUNDS(SIZE(a_ij, 2),n_s)

    $CHECK_BOUNDS(SIZE(c_i),n_s)

    ! Calculate the intermediate-state matrices for the n_s-step IRK

    ! First, evaluate the RHS matrices

    dx = x_b - x_a
    
    do i = 1, n_s
       A(:,:,i) = eq%A(x_a + dx*c_i(i), omega)
    end do

    ! Set up the linear problem for K

    do i = 1, n_s

       i_a = (i-1)*eq%n_e + 1
       i_b = (i-1)*eq%n_e + eq%n_e

       do j = 1, n_s

          j_a = (j-1)*eq%n_e + 1
          j_b = (j-1)*eq%n_e + eq%n_e

          M(i_a:i_b,j_a:j_b) = -dx*a_ij(i,j)*A(:,:,i)

          if (i == j) then
             M(i_a:i_b,j_a:j_b) = M(i_a:i_b,j_a:j_b) + identity_matrix(eq%n_e)
          endif

       end do

       B(i_a:i_b,:) = A(:,:,i)

    end do

    ! Solve the linear system

    X = linear_solve(M, B)

    ! Extract the K matrices

    do i = 1, n_s

       i_a = (i-1)*eq%n_e + 1
       i_b = (i-1)*eq%n_e + eq%n_e

       K(:,:,i) = X(i_a:i_b,:)

    end do

    ! Finish

    return

  end function irk_K_

end module gyre_${T}_colloc_ivp

$endsub

