! Incfile  : gyre_rot
! Purpose  : rotational effects (template)
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
     procedure(omega_), deferred, public :: omega
     procedure(omega_c_1_), deferred     :: omega_c_1_
     procedure(omega_c_v_), deferred     :: omega_c_v_
     generic, public                     :: omega_c => omega_c_1_, omega_c_v_
     procedure(l_i_), deferred, public   :: l_i
     procedure(l_e_1_), deferred         :: l_e_1_
     procedure(l_e_v_), deferred         :: l_e_v_
     generic, public                     :: l_e => l_e_1_, l_e_v_
     procedure(lambda_1_), deferred      :: lambda_1_
     procedure(lambda_v_), deferred      :: lambda_v_
     generic, public                     :: lambda => lambda_1_, lambda_v_
  end type ${T}_rot_t

  ! Interfaces

  abstract interface

     function omega_ (this, s, x, omega_c) result (omega)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: s
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega_c
       $TYPE(WP)                     :: omega
     end function omega_

     function omega_c_1_ (this, s, x, omega) result (omega_c)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: s
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: omega_c
     end function omega_c_1_

     function omega_c_v_ (this, s, x, omega) result (omega_c)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: s(:)
       real(WP), intent(in)          :: x(:)
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: omega_c(SIZE(s))
     end function omega_c_v_

     function l_i_ (this, omega) result (l_i)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: l_i
     end function l_i_
     
     function l_e_1_ (this, s, x, omega) result (l_e)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: s
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: l_e
     end function l_e_1_
     
     function l_e_v_ (this, s, x, omega) result (l_e)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: s(:)
       real(WP), intent(in)          :: x(:)
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: l_e(SIZE(s))
     end function l_e_v_
     
     function lambda_1_ (this, s, x, omega) result (lambda)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: s
       real(WP), intent(in)          :: x
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: lambda
     end function lambda_1_
     
     function lambda_v_ (this, s, x, omega) result (lambda)
       use core_kinds
       import ${T}_rot_t
       class(${T}_rot_t), intent(in) :: this
       integer, intent(in)           :: s(:)
       real(WP), intent(in)          :: x(:)
       $TYPE(WP), intent(in)         :: omega
       $TYPE(WP)                     :: lambda(SIZE(s))
     end function lambda_v_
     
  end interface

  ! Access specifiers

  private

  public :: ${T}_rot_t

end module gyre_${T}_rot

$endsub
