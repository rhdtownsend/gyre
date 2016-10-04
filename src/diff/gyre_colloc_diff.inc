! Incfile  : gyre_colloc_diff
! Purpose  : difference equations (local collocation, template)
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

module gyre_${T}_colloc_diff

  ! Uses

  use core_kinds
  use core_linalg

  use gyre_diff
  use gyre_eqns
  use gyre_ext
  use gyre_linalg
  use gyre_point

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: GL2_SCHEME = 1
  integer, parameter :: GL4_SCHEME = 2
  integer, parameter :: GL6_SCHEME = 3

  ! Derived-type definitions

  type, extends (${T}_diff_t) :: ${T}_colloc_diff_t
     private
     class(${T}_eqns_t), allocatable :: eq
     type(point_t)                   :: pt_a
     type(point_t)                   :: pt_b
     integer                         :: scheme
   contains
     private
     procedure, public :: build
     procedure         :: build_GL2_
     procedure         :: build_GL4_
     procedure         :: build_GL6_
     procedure         :: build_irk_
  end type ${T}_colloc_diff_t

  ! Interfaces

  interface ${T}_colloc_diff_t
     module procedure ${T}_colloc_diff_t_
  end interface ${T}_colloc_diff_t

  ! Access specifiers

  private

  public :: ${T}_colloc_diff_t

contains

  function ${T}_colloc_diff_t_ (eq, pt_a, pt_b, scheme) result (df)

    class(${T}_eqns_t), intent(in) :: eq
    type(point_t), intent(in)      :: pt_a
    type(point_t), intent(in)      :: pt_b
    character(*), intent(in)       :: scheme
    type(${T}_colloc_diff_t)       :: df

    $ASSERT_DEBUG(pt_a%s == pt_b%s,Segment mismatch)

    ! Construct the colloc_diff_t

    allocate(df%eq, SOURCE=eq)

    df%pt_a = pt_a
    df%pt_b = pt_b
       
    select case (scheme)
    case ('GL2')
       df%scheme = GL2_SCHEME
    case ('GL4')
       df%scheme = GL4_SCHEME
    case ('GL6')
       df%scheme = GL6_SCHEME
    case default
       $ABORT(Invalid scheme)
    end select

    df%n_e = eq%n_e

    ! Finish

    return
    
  end function ${T}_colloc_diff_t_

  !****

  subroutine build (this, omega, E_l, E_r, scl)

    class(${T}_colloc_diff_t), intent(in) :: this
    $TYPE(WP), intent(in)                 :: omega
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: scl

    ! Build the difference equations

    select case (this%scheme)
    case (GL2_SCHEME)
       call this%build_GL2_(omega, E_l, E_r, scl)
    case (GL4_SCHEME)
       call this%build_GL4_(omega, E_l, E_r, scl)
    case (GL6_SCHEME)
       call this%build_GL6_(omega, E_l, E_r, scl)
    case default
       $ABORT(Invalid scheme)
    end select

    ! Finish

    return

  end subroutine build
  
  !****

  subroutine build_GL2_ (this, omega, E_l, E_r, scl)

    class(${T}_colloc_diff_t), intent(in) :: this
    $TYPE(WP), intent(in)                 :: omega
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: scl

    real(WP)      :: dx
    type(point_t) :: pt
    $TYPE(WP)     :: A(this%n_e,this%n_e,1)

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Evaluate the RHS matrix

    dx = this%pt_b%x - this%pt_a%x

    pt%s = this%pt_a%s
    pt%x = this%pt_a%x + dx

    A(:,:,1) = this%eq%A(pt, omega)

    ! Build the GL2 difference equations

    E_l = 0.5_WP*dx*A(:,:,1) + identity_matrix(this%n_e)
    E_r = 0.5_WP*dx*A(:,:,1) - identity_matrix(this%n_e)

    scl = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine build_GL2_

  !****

  subroutine build_GL4_ (this, omega, E_l, E_r, scl)

    class(${T}_colloc_diff_t), intent(in) :: this
    $TYPE(WP), intent(in)                 :: omega
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: scl

    logical :: BALANCE = .FALSE.

    real(WP), parameter :: A_IJ(2,2) = RESHAPE([0.25_WP,&
                                                0.25_WP+SQRT(3._WP)/6._WP, &
                                                0.25_WP-SQRT(3._WP)/6._WP, &
                                                0.25_WP], [2,2])
    real(WP), parameter :: B_I(2) = [0.5_WP, &
                                     0.5_WP]
    real(WP), parameter :: C_I(2) = [0.5_WP-SQRT(3._WP)/6._WP, &
         0.5_WP+SQRT(3._WP)/6._WP]

    $TYPE(WP) :: E_l_a(this%n_e,this%n_e)
    $TYPE(WP) :: E_r_a(this%n_e,this%n_e)
    $TYPE(WP) :: E_l_b(this%n_e,this%n_e)
    $TYPE(WP) :: E_r_b(this%n_e,this%n_e)

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)
 
    ! Build the GL4 difference equations

    if (BALANCE) then
       
       call this%build_irk_(omega, 2, this%pt_a, this%pt_b, A_IJ, B_I, C_I, E_l_a, E_r_a)
       call this%build_irk_(omega, 2, this%pt_b, this%pt_a, A_IJ, B_I, C_I, E_l_b, E_r_b)

       E_l =  E_l_a - E_r_b
       E_r = -E_l_b + E_r_a

    else

       if (this%pt_a%x < 0.5_WP) then
          call this%build_irk_(omega, 2, this%pt_a, this%pt_b, A_IJ, B_I, C_I, E_l, E_r)
       else
          call this%build_irk_(omega, 2, this%pt_b, this%pt_a, A_IJ, B_I, C_I, E_r, E_l)
       endif

    endif

    scl = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine build_GL4_

  !****

  subroutine build_GL6_ (this, omega, E_l, E_r, scl)

    class(${T}_colloc_diff_t), intent(in) :: this
    $TYPE(WP), intent(in)                 :: omega
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: scl

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

    ! Build the GL6 difference equations

    if (this%pt_a%x < 0.5_WP) then
       call this%build_irk_(omega, 3, this%pt_a, this%pt_b, A_IJ, B_I, C_I, E_l, E_r)
    else
       call this%build_irk_(omega, 3, this%pt_b, this%pt_a, A_IJ, B_I, C_I, E_r, E_l)
    endif

    scl = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine build_GL6_

  !****

  subroutine build_irk_ (this, omega, n_s, pt_a, pt_b, a_ij, b_i, c_i, E_a, E_b)

    class(${T}_colloc_diff_t), intent(in) :: this
    $TYPE(WP), intent(in)                 :: omega
    integer, intent(in)                   :: n_s
    type(point_t), intent(in)             :: pt_a
    type(point_t), intent(in)             :: pt_b
    real(WP), intent(in)                  :: a_ij(:,:)
    real(WP), intent(in)                  :: b_i(:)
    real(WP), intent(in)                  :: c_i(:)
    $TYPE(WP), intent(out)                :: E_a(:,:)
    $TYPE(WP), intent(out)                :: E_b(:,:)

    real(WP)      :: dx
    type(point_t) :: pt
    $TYPE(WP)     :: A(this%n_e,this%n_e,n_s)
    $TYPE(WP)     :: M(n_s*this%n_e,n_s*this%n_e)
    integer       :: i
    integer       :: i_a
    integer       :: i_b
    integer       :: j
    integer       :: j_a
    integer       :: j_b
    $TYPE(WP)     :: B(n_s*this%n_e,this%n_e)
    $TYPE(WP)     :: X(n_s*this%n_e,this%n_e)
    $TYPE(WP)     :: K(this%n_e,this%n_e,n_s)
    
    $CHECK_BOUNDS(SIZE(a_ij, 1),n_s)
    $CHECK_BOUNDS(SIZE(a_ij, 2),n_s)

    $CHECK_BOUNDS(SIZE(b_i),n_s)

    $CHECK_BOUNDS(SIZE(c_i),n_s)

    $CHECK_BOUNDS(SIZE(E_a, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_a, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_b, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_b, 2),this%n_e)

    ! Build the difference equations using an n_s-step implicit Runge
    ! Kutta (IRK) scheme with Butcher tableaux coefficients a_ij, b_i,
    ! and c_i

    ! Evaluate the RHS matrices

    dx = pt_b%x - pt_a%x

    do i = 1, n_s

       pt%s = pt_a%s
       pt%x = pt_a%x + dx*c_i(i)

       A(:,:,i) = this%eq%A(pt, omega)

    end do

    ! Calculate the intermediate-state matrices K

    ! First, set up the linear system

    do i = 1, n_s

       i_a = (i-1)*this%n_e + 1
       i_b = (i-1)*this%n_e + this%n_e

       do j = 1, n_s

          j_a = (j-1)*this%n_e + 1
          j_b = (j-1)*this%n_e + this%n_e

          M(i_a:i_b,j_a:j_b) = -dx*a_ij(i,j)*A(:,:,i)

          if (i == j) then
             M(i_a:i_b,j_a:j_b) = M(i_a:i_b,j_a:j_b) + identity_matrix(this%n_e)
          endif

       end do

       B(i_a:i_b,:) = A(:,:,i)

    end do

    ! Solve the linear system

    X = linear_solve(M, B)

    ! Extract the K matrices

    do i = 1, n_s

       i_a = (i-1)*this%n_e + 1
       i_b = (i-1)*this%n_e + this%n_e

       K(:,:,i) = X(i_a:i_b,:)

    end do

    ! Build the difference equations

    E_a = identity_matrix(this%n_e)

    do i = 1, n_s
       E_a = E_a + dx*b_i(i)*K(:,:,i)
    enddo

    E_b = -identity_matrix(this%n_e)

    ! Finish

    return

  end subroutine build_irk_

end module gyre_${T}_colloc_diff

$endsub

