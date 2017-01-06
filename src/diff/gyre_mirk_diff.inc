! Incfile  : gyre_mirk_diff
! Purpose  : difference equations (mono-implicit Runge-Kutta, template)
!
! Copyright 2017 Rich Townsend
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

module gyre_${T}_mirk_diff

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

  type, extends (${T}_diff_t) :: ${T}_mirk_diff_t
     private
     class(${T}_eqns_t), allocatable :: eq
     type(point_t)                   :: pt_a
     type(point_t)                   :: pt_b
   contains
     private
     procedure, public :: build
  end type ${T}_mirk_diff_t

  ! Interfaces

  interface ${T}_mirk_diff_t
     module procedure ${T}_mirk_diff_t_
  end interface ${T}_mirk_diff_t

  ! Access specifiers

  private

  public :: ${T}_mirk_diff_t

contains

  function ${T}_mirk_diff_t_ (eq, pt_a, pt_b) result (df)

    class(${T}_eqns_t), intent(in) :: eq
    type(point_t), intent(in)      :: pt_a
    type(point_t), intent(in)      :: pt_b
    type(${T}_mirk_diff_t)         :: df

    $ASSERT_DEBUG(pt_a%s == pt_b%s,Segment mismatch)

    ! Construct the mirk_diff_t

    allocate(df%eq, SOURCE=eq)

    df%pt_a = pt_a
    df%pt_b = pt_b

    df%n_e = eq%n_e

    ! Finish

    return
    
  end function ${T}_mirk_diff_t_

  !****

  subroutine build (this, omega, E_l, E_r, scl)

    class(${T}_mirk_diff_t), intent(in) :: this
    $TYPE(WP), intent(in)               :: omega
    $TYPE(WP), intent(out)              :: E_l(:,:)
    $TYPE(WP), intent(out)              :: E_r(:,:)
    type(${T}_ext_t), intent(out)       :: scl

    real(WP)      :: dx
    type(point_t) :: pt
    $TYPE(WP)     :: A(this%n_e,this%n_e,3)

    $CHECK_BOUNDS(SIZE(E_l, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_l, 2),this%n_e)

    $CHECK_BOUNDS(SIZE(E_r, 1),this%n_e)
    $CHECK_BOUNDS(SIZE(E_r, 2),this%n_e)

    if (this%pt_a%x == 0) then

       ! If we're at the origin, fudge things to avoid evaluating A there

       E_l = -identity_matrix(this%n_e)
       E_r =  identity_matrix(this%n_e)

       scl = ${T}_ext_t(1._WP)

    else
       
       ! Evaluate the RHS matrices

       dx = this%pt_b%x - this%pt_a%x

       A(:,:,1) = this%eq%A(this%pt_a, omega)

       pt%s = this%pt_a%s
       pt%x = this%pt_a%x + 0.5_WP*dx
       
       A(:,:,2) = this%eq%A(pt, omega)
    
       A(:,:,3) = this%eq%A(this%pt_b, omega)

       ! Build the difference equations

       E_l = -identity_matrix(this%n_e) - dx*A(:,:,1)/6._WP - dx*A(:,:,2)/3._WP - dx**2*MATMUL(A(:,:,2), A(:,:,1))/12._WP
       E_r =  identity_matrix(this%n_e) - dx*A(:,:,3)/6._WP - dx*A(:,:,2)/3._WP + dx**2*MATMUL(A(:,:,2), A(:,:,3))/12._WP
       
       scl = ${T}_ext_t(1._WP)

    endif

    ! Finish

  end subroutine build
  
end module gyre_${T}_mirk_diff

$endsub

