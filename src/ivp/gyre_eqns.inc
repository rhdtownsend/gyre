! Incfile  : gyre_eqns
! Purpose  : differential equations evaluation (template)
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_eqns

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ${T}_eqns_t
     private
     integer, public :: n_e
   contains
     private
     procedure(A_), deferred, public :: A
     procedure(A_), deferred, public :: xA
  end type ${T}_eqns_t

  ! Interfaces

  abstract interface

     function A_ (this, x, omega) result (A)
       use core_kinds
       import ${T}_eqns_t
       class(${T}_eqns_t), intent(in) :: this
       real(WP), intent(in)           :: x
       $TYPE(WP), intent(in)          :: omega
       $TYPE(WP)                      :: A(this%n_e,this%n_e)
     end function A_

  end interface

  ! Access specifiers

  private

  public :: ${T}_eqns_t

end module gyre_${T}_eqns

$endsub
