! Incfile  : gyre_eqns
! Purpose  : differential equations evaluation (template)
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
     integer :: n_e
   contains
     procedure(stencil), deferred :: stencil
     procedure(A), deferred       :: A
     procedure(xA), deferred      :: xA
  end type ${T}_eqns_t

  ! Interfaces

  abstract interface

     subroutine stencil (this, pt)
       use gyre_point
       import ${T}_eqns_t
       class(${T}_eqns_t), intent(inout) :: this
       type(point_t), intent(in)         :: pt(:)
     end subroutine stencil
     
     function A (this, i, st)
       use core_kinds
       use gyre_state
       import ${T}_eqns_t
       class(${T}_eqns_t), intent(in)  :: this
       integer, intent(in)             :: i
       class(${T}_state_t), intent(in) :: st
       $TYPE(WP)                       :: A(this%n_e,this%n_e)
     end function A

     function xA (this, i, st)
       use core_kinds
       use gyre_state
       import ${T}_eqns_t
       class(${T}_eqns_t), intent(in)  :: this
       integer, intent(in)             :: i
       class(${T}_state_t), intent(in) :: st
       $TYPE(WP)                       :: xA(this%n_e,this%n_e)
     end function xA

  end interface

  ! Access specifiers

  private

  public :: ${T}_eqns_t

end module gyre_${T}_eqns

$endsub
