! Incfile  : gyre_bound
! Purpose  : boundary conditions (template)
!
! Copyright 2015 Rich Townsend
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

module gyre_${T}_bound

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ${T}_bound_t
     private
     integer, public :: n
     integer, public :: n_e
   contains
     private
     procedure(B_), deferred, public :: B
   end type ${T}_bound_t

  ! Interfaces

  abstract interface

     function B_ (this, omega) result (B)
       use core_kinds
       import ${T}_bound_t
       class(${T}_bound_t), intent(in) :: this
       $TYPE(WP), intent(in)           :: omega
       $TYPE(WP)                       :: B(this%n,this%n_e)
     end function B_

  end interface

  ! Access specifiers

  private

  public :: ${T}_bound_t

end module gyre_${T}_bound

$endsub
