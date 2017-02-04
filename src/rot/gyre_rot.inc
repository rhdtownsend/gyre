! Incfile  : gyre_rot
! Purpose  : rotational effects (template)
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

module gyre_${T}_rot

  ! Uses

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ${T}_rot_t
   contains
     procedure(stencil), deferred :: stencil
     procedure(omega), deferred   :: omega
     procedure(omega_c), deferred :: omega_c
     procedure(l_i), deferred     :: l_i
     procedure(l_e), deferred     :: l_e
     procedure(lambda), deferred  :: lambda
  end type ${T}_rot_t

  ! Interfaces

  abstract interface

     subroutine stencil (this, pt)
       use gyre_point
       import ${T}_rot_t
       class(${T}_rot_t), intent(inout) :: this
       type(point_t), intent(in)        :: pt(:)
     end subroutine stencil

     function omega (this, i, omega_c)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: i
       $TYPE(WP), intent(in)         :: omega_c
       $TYPE(WP)                     :: omega
     end function omega

     function omega_c (this, i, omega)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: i
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: omega_c
     end function omega_c

     function l_i (this, omega)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: l_i
     end function l_i
     
     function l_e (this, i, omega)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: i
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: l_e
     end function l_e
     
     function lambda (this, i, omega)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: i
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: lambda
     end function lambda
     
  end interface

  ! Access specifiers

  private

  public :: ${T}_rot_t

end module gyre_${T}_rot

$endsub
