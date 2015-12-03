! Incfile  : gyre_match
! Purpose  : match conditions (template)
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

module gyre_${T}_match

  ! Uses

  use core_kinds

  use gyre_ext

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ${T}_match_t
     private
     integer, public :: n_e
   contains
     private
     procedure(build_), deferred, public :: build
   end type ${T}_match_t

  ! Interfaces

  abstract interface

     subroutine build_ (this, omega, E_l, E_r, S)
       use core_kinds
       import ${T}_match_t
       class(${T}_match_t), target, intent(in) :: this
       $TYPE(WP), intent(in)                   :: omega
       $TYPE(WP), intent(out)                  :: E_l(:,:)
       $TYPE(WP), intent(out)                  :: E_r(:,:)
       type(${T}_ext_t), intent(out)           :: S
     end subroutine build_

  end interface

  ! Access specifiers

  private

  public :: ${T}_match_t

end module gyre_${T}_match

$endsub
