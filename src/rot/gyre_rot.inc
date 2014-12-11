! Incfile  : gyre_rot
! Purpose  : rotational effects (template)
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

module gyre_${T}_rot

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ${T}_rot_t
   contains
     private
     procedure(l_e_), deferred, public     :: l_e
     procedure(l_0_), deferred, public     :: l_0
     procedure(omega_), deferred, public   :: omega
     procedure(omega_c_), deferred, public :: omega_c
  end type ${T}_rot_t

  ! Interfaces

  abstract interface

     function l_e_ (this, x, omega) result (l_e)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: l_e
     end function l_e_
     
     function l_0_ (this, omega) result (l_0)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: l_0
     end function l_0_
     
     function omega_ (this, x, omega_c) result (omega)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega_c
       $TYPE(WP)                     :: omega
     end function omega_

     function omega_c_ (this, x, omega) result (omega_c)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: omega_c
     end function omega_c_

  end interface

  ! Access specifiers

  private

  public :: ${T}_rot_t

end module gyre_${T}_rot

$endsub
