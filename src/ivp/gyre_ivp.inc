! Incfile  : gyre_ivp
! Purpose  : initial-value solvers (template)
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

module gyre_${T}_ivp

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ${T}_ivp_t
     private
     integer, public :: n_e
   contains
     private
     procedure(shoot_), deferred, public    :: shoot
     procedure(recon_), deferred, public    :: recon
     procedure(abscissa_), deferred, public :: abscissa
  end type ${T}_ivp_t

  ! Interfaces

  abstract interface

     subroutine shoot_ (this, omega, x_a, x_b, E_l, E_r, S)
       use core_kinds
       use gyre_ext
       import ${T}_ivp_t
       class(${T}_ivp_t), intent(in) :: this
       $TYPE(WP), intent(in)         :: omega
       real(WP), intent(in)          :: x_a
       real(WP), intent(in)          :: x_b
       $TYPE(WP), intent(out)        :: E_l(:,:)
       $TYPE(WP), intent(out)        :: E_r(:,:)
       type(${T}_ext_t), intent(out) :: S
     end subroutine shoot_

     subroutine recon_ (this, omega, x_a, x_b, y_a, y_b, x, y)
       use core_kinds
       import ${T}_ivp_t
       class(${T}_ivp_t), intent(in) :: this
       $TYPE(WP), intent(in)         :: omega
       real(WP), intent(in)          :: x_a
       real(WP), intent(in)          :: x_b
       $TYPE(WP), intent(in)         :: y_a(:)
       $TYPE(WP), intent(in)         :: y_b(:)
       real(WP), intent(in)          :: x(:)
       $TYPE(WP), intent(out)        :: y(:,:)
     end subroutine recon_

     function abscissa_ (this, x_a, x_b) result (x)
       use core_kinds
       import ${T}_ivp_t
       class(${T}_ivp_t), intent(in) :: this
       real(WP), intent(in)          :: x_a
       real(WP), intent(in)          :: x_b
       real(WP), allocatable         :: x(:)
     end function abscissa_

  end interface

  ! Access specifiers

  private

  public :: ${T}_ivp_t

end module gyre_${T}_ivp

$endsub
