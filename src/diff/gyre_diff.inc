! Incfile  : gyre_diff
! Purpose  : difference equations (template)
!
! Copyright 2015-2017 Rich Townsend
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

module gyre_${T}_diff

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ${T}_diff_t
     integer :: n_e
   contains
     procedure(build), deferred :: build
  end type ${T}_diff_t

  ! Interfaces

  abstract interface

     subroutine build (this, st, E_l, E_r, scl)
       use core_kinds
       use gyre_ext
       use gyre_state
       import ${T}_diff_t
       class(${T}_diff_t), intent(in)  :: this
       class(${T}_state_t), intent(in) :: st
       $TYPE(WP), intent(out)          :: E_l(:,:)
       $TYPE(WP), intent(out)          :: E_r(:,:)
       type(${T}_ext_t), intent(out)   :: scl
     end subroutine build

  end interface

  ! Access specifiers

  private

  public :: ${T}_diff_t

end module gyre_${T}_diff

$endsub
