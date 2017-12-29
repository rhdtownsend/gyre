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
     procedure(l_e), deferred    :: l_e
     procedure(lambda), deferred :: lambda
  end type ${T}_rot_t

  ! Interfaces

  abstract interface

     function l_e (this, Omega_rot, omega)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: Omega_rot
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: l_e
     end function l_e
     
     function lambda (this, Omega_rot, omega)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: Omega_rot
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: lambda
     end function lambda
     
  end interface

  ! Access specifiers

  private

  public :: ${T}_rot_t

end module gyre_${T}_rot

$endsub
