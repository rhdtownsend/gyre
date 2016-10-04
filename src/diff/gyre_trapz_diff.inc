! Incfile  : gyre_trapz_diff
! Purpose  : difference equations (quasi-trapezoidal, template)
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

module gyre_${T}_trapz_diff

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

  ! Derived-type definitions

  type, extends (${T}_diff_t) :: ${T}_trapz_diff_t
     private
     class(${T}_eqns_t), allocatable :: eq
     type(point_t)                   :: pt_a
     type(point_t)                   :: pt_b
     real(WP), allocatable           :: w(:)
   contains
     private
     procedure, public :: build
  end type ${T}_trapz_diff_t

  ! Interfaces

  interface ${T}_trapz_diff_t
     module procedure ${T}_trapz_diff_t_
  end interface ${T}_trapz_diff_t

  ! Access specifiers

  private

  public :: ${T}_trapz_diff_t

contains

  function ${T}_trapz_diff_t_ (eq, pt_a, pt_b, w) result (df)

    class(${T}_eqns_t), intent(in) :: eq
    type(point_t), intent(in)      :: pt_a
    type(point_t), intent(in)      :: pt_b
    real(WP), intent(in)           :: w(:)
    type(${T}_trapz_diff_t)        :: df

    $ASSERT_DEBUG(pt_a%s == pt_b%s,Segment mismatch)

    $CHECK_BOUNDS(SIZE(w),eq%n_e)

    ! Construct the trapz_diff_t

    allocate(df%eq, SOURCE=eq)

    df%pt_a = pt_a
    df%pt_b = pt_b

    df%w = w

    df%n_e = eq%n_e

    ! Finish

    return
    
  end function ${T}_trapz_diff_t_

  !****

  subroutine build (this, omega, E_l, E_r, scl)

    class(${T}_trapz_diff_t), intent(in) :: this
    $TYPE(WP), intent(in)                :: omega
    $TYPE(WP), intent(out)               :: E_l(:,:)
    $TYPE(WP), intent(out)               :: E_r(:,:)
    type(${T}_ext_t), intent(out)        :: scl

    real(WP)  :: dx
    $TYPE(WP) :: xA(this%n_e,this%n_e,2)
    real(WP)  :: x_m

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    ! Evaluate the RHS matrix

    dx = this%pt_b%x - this%pt_a%x

    xA(:,:,1) = this%eq%xA(this%pt_a, omega)
    xA(:,:,2) = this%eq%xA(this%pt_b, omega)

    ! Build the difference equations

    x_m = 0.5_WP*(this%pt_a%x + this%pt_b%x)

    E_l = -x_m*identity_matrix(this%n_e) - dx*MATMUL(diagonal_matrix(1._WP-this%w), xA(:,:,1))
    E_r =  x_m*identity_matrix(this%n_e) - dx*MATMUL(diagonal_matrix(      this%w), xA(:,:,2))

    scl = ${T}_ext_t(1._WP)

    ! Finish

  end subroutine build
  
end module gyre_${T}_trapz_diff

$endsub

