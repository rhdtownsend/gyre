! Incfile  : gyre_colloc_diff
! Purpose  : difference equations (local collocation, template)
!
! Copyright 2013-2017 Rich Townsend
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
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: GL2_SCHEME = 1
  integer, parameter :: GL4_SCHEME = 2
  integer, parameter :: GL6_SCHEME = 3

  real(WP), parameter :: GL2_C_I(1) = [0.5_WP]

  real(WP), parameter :: GL4_A_IJ(2,2) = RESHAPE([ &
       0.25_WP,&
       0.25_WP+SQRT(3._WP)/6._WP, &
       0.25_WP-SQRT(3._WP)/6._WP, &
       0.25_WP], &
       [2,2])
  real(WP), parameter :: GL4_B_I(2) = [ &
       0.5_WP, &
       0.5_WP]
  real(WP), parameter :: GL4_C_I(2) = [ &
       0.5_WP-SQRT(3._WP)/6._WP, &
       0.5_WP+SQRT(3._WP)/6._WP]
  
  real(WP), parameter :: GL6_A_IJ(3,3) = RESHAPE([ &
       5._WP/36._WP, &
       5._WP/36._WP+SQRT(15._WP)/24._WP, &
       5._WP/36._WP+SQRT(15._WP)/30._WP, &
       2._WP/9._WP-SQRT(15._WP)/15._WP, &
       2._WP/9._WP, &
       2._WP/9._WP+SQRT(15._WP)/15._WP, &
       5._WP/36._WP-SQRT(15._WP)/30._WP, &
       5._WP/36._WP-SQRT(15._WP)/24._WP, &
       5._WP/36._WP], &
       [3,3])
  real(WP), parameter :: GL6_B_I(3) = [ &
       5._WP/18._WP, &
       4._WP/9._WP, &
       5._WP/18._WP]
  real(WP), parameter :: GL6_C_I(3) = [ &
       0.5_WP-SQRT(15._WP)/10._WP, &
       0.5_WP, &
       0.5_WP+SQRT(15._WP)/10._WP]
       
  ! Derived-type definitions

  type, extends (${T}_diff_t) :: ${T}_colloc_diff_t
     private
     class(${T}_eqns_t), allocatable :: eq
     real(WP)                        :: dx
     integer                         :: n_s
     integer                         :: scheme
     logical                         :: reverse
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

    real(WP), allocatable      :: c_i(:)
    type(point_t), allocatable :: pt(:)

    $ASSERT_DEBUG(pt_a%s == pt_b%s,Segment mismatch)

    ! Construct the colloc_diff_t

    select case (scheme)
    case ('GL2')
       df%scheme = GL2_SCHEME
       df%n_s = 1
       c_i = GL2_C_I
    case ('GL4')
       df%scheme = GL4_SCHEME
       df%n_s = 2
       c_i = GL4_C_I
    case ('GL6')
       df%scheme = GL6_SCHEME
       df%n_s = 3
       c_i = GL6_C_I
    case default
       $ABORT(Invalid scheme)
    end select

    df%reverse = pt_a%x >= 0.5_WP

    allocate(pt(df%n_s))

    pt%s = pt_a%s

    if (df%reverse) then
       df%dx = pt_a%x - pt_b%x
       pt%x = pt_b%x + c_i*df%dx
    else
       df%dx = pt_b%x - pt_a%x
       pt%x = pt_a%x + c_i*df%dx
    endif

    allocate(df%eq, SOURCE=eq)

    call df%eq%stencil(pt)

    df%n_e = eq%n_e

    ! Finish

    return
    
  end function ${T}_colloc_diff_t_

  !****

  subroutine build (this, st, E_l, E_r, scl)

    class(${T}_colloc_diff_t), intent(in) :: this
    class(${T}_state_t), intent(in)       :: st
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: scl

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Build the difference equations

    select case (this%scheme)
    case (GL2_SCHEME)
       call this%build_GL2_(st, E_l, E_r, scl)
    case (GL4_SCHEME)
       call this%build_GL4_(st, E_l, E_r, scl)
    case (GL6_SCHEME)
       call this%build_GL6_(st, E_l, E_r, scl)
    case default
       $ABORT(Invalid scheme)
    end select

    ! Finish

    return

  end subroutine build
  
  !****

  subroutine build_GL2_ (this, st, E_l, E_r, scl)

    class(${T}_colloc_diff_t), intent(in) :: this
    class(${T}_state_t), intent(in)       :: st
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: scl

    $TYPE(WP) :: A(this%n_e,this%n_e)

    ! Evaluate the RHS matrix

    A = this%eq%A(1, st)

    ! Build the GL2 difference equations

    if (this%reverse) then
       E_l = 0.5_WP*this%dx*A - identity_matrix(this%n_e)
       E_r = 0.5_WP*this%dx*A + identity_matrix(this%n_e)
    else
       E_l = 0.5_WP*this%dx*A + identity_matrix(this%n_e)
       E_r = 0.5_WP*this%dx*A - identity_matrix(this%n_e)
    endif

    scl = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine build_GL2_

  !****

  subroutine build_GL4_ (this, st, E_l, E_r, scl)

    class(${T}_colloc_diff_t), intent(in) :: this
    class(${T}_state_t), intent(in)       :: st
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: scl

    ! Build the GL4 difference equations

    if (this%reverse) then
       call this%build_irk_(st, GL4_A_IJ, GL4_B_I, E_r, E_l)
    else
       call this%build_irk_(st, GL4_A_IJ, GL4_B_I, E_l, E_r)
    endif

    scl = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine build_GL4_

  !****

  subroutine build_GL6_ (this, st, E_l, E_r, scl)

    class(${T}_colloc_diff_t), intent(in) :: this
    class(${T}_state_t), intent(in)       :: st
    $TYPE(WP), intent(out)                :: E_l(:,:)
    $TYPE(WP), intent(out)                :: E_r(:,:)
    type(${T}_ext_t), intent(out)         :: scl

    ! Build the GL6 difference equations

    if (this%reverse) then
       call this%build_irk_(st, GL6_A_IJ, GL6_B_I, E_r, E_l)
    else
       call this%build_irk_(st, GL6_A_IJ, GL6_B_I, E_l, E_r)
    endif

    scl = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine build_GL6_

  !****

  subroutine build_irk_ (this, st, a_ij, b_i, E_a, E_b)

    class(${T}_colloc_diff_t), intent(in) :: this
    class(${T}_state_t), intent(in)       :: st
    real(WP), intent(in)                  :: a_ij(:,:)
    real(WP), intent(in)                  :: b_i(:)
    $TYPE(WP), intent(out)                :: E_a(:,:)
    $TYPE(WP), intent(out)                :: E_b(:,:)

    $TYPE(WP) :: A(this%n_e,this%n_e,this%n_s)
    $TYPE(WP) :: M(this%n_s*this%n_e,this%n_s*this%n_e)
    integer   :: i
    integer   :: i_a
    integer   :: i_b
    integer   :: j
    integer   :: j_a
    integer   :: j_b
    $TYPE(WP) :: B(this%n_s*this%n_e,this%n_e)
    $TYPE(WP) :: X(this%n_s*this%n_e,this%n_e)
    $TYPE(WP) :: K(this%n_e,this%n_e,this%n_s)
    
    $CHECK_BOUNDS(SIZE(a_ij, 1),this%n_s)
    $CHECK_BOUNDS(SIZE(a_ij, 2),this%n_s)

    $CHECK_BOUNDS(SIZE(b_i),this%n_s)

    ! Build the difference equations using an n_s-step implicit Runge
    ! Kutta (IRK) scheme with Butcher tableaux coefficients a_ij and
    ! b_i (c_i has already been used in the constructor)

    ! Evaluate the RHS matrices

    do i = 1, this%n_s
       A(:,:,i) = this%eq%A(i, st)
    end do

    ! Calculate the intermediate-state matrices K

    ! First, set up the linear system

    do i = 1, this%n_s

       i_a = (i-1)*this%n_e + 1
       i_b = (i-1)*this%n_e + this%n_e

       do j = 1, this%n_s

          j_a = (j-1)*this%n_e + 1
          j_b = (j-1)*this%n_e + this%n_e

          M(i_a:i_b,j_a:j_b) = -this%dx*a_ij(i,j)*A(:,:,i)

          if (i == j) then
             M(i_a:i_b,j_a:j_b) = M(i_a:i_b,j_a:j_b) + identity_matrix(this%n_e)
          endif

       end do

       B(i_a:i_b,:) = A(:,:,i)

    end do

    ! Solve the linear system

    X = linear_solve(M, B)

    ! Extract the K matrices

    do i = 1, this%n_s

       i_a = (i-1)*this%n_e + 1
       i_b = (i-1)*this%n_e + this%n_e

       K(:,:,i) = X(i_a:i_b,:)

    end do

    ! Build the difference equations

    E_a = identity_matrix(this%n_e)

    do i = 1, this%n_s
       E_a = E_a + this%dx*b_i(i)*K(:,:,i)
    enddo

    E_b = -identity_matrix(this%n_e)

    ! Finish

    return

  end subroutine build_irk_

end module gyre_${T}_colloc_diff

$endsub

