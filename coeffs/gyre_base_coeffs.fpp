! Module   : gyre_base_coeffs
! Purpose  : base structure coefficients (interface)
!
! Copyright 2013 Rich Townsend
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

module gyre_base_coeffs

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure(c_1_i), deferred :: ${NAME}_1
    procedure(c_v_i), deferred :: ${NAME}_v
    generic, public            :: ${NAME} => ${NAME}_1, ${NAME}_v
  $endsub

  type, abstract :: base_coeffs_t
   contains
     private
     $PROC_DECL(V)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(delta)
     procedure(pi_c_i), deferred, public         :: pi_c
     procedure(enable_cache_i), deferred, public :: enable_cache
     procedure(enable_cache_i), deferred, public :: disable_cache
     procedure(fill_cache_i), deferred, public   :: fill_cache
  end type base_coeffs_t

  ! Interfaces

  abstract interface

     function c_1_i (this, x) result (c)
       use core_kinds
       import base_coeffs_t
       class(base_coeffs_t), intent(in) :: this
       real(WP), intent(in)             :: x
       real(WP)                         :: c
     end function c_1_i

     function c_v_i (this, x) result (c)
       use core_kinds
       import base_coeffs_t
       class(base_coeffs_t), intent(in) :: this
       real(WP), intent(in)             :: x(:)
       real(WP)                         :: c(SIZE(x))
     end function c_v_i

     function pi_c_i (this) result (pi_c)
       use core_kinds
       import base_coeffs_t
       class(base_coeffs_t), intent(in) :: this
       real(WP)                         :: pi_c
     end function pi_c_i

     subroutine enable_cache_i (this)
       import base_coeffs_t
       class(base_coeffs_t), intent(inout) :: this
     end subroutine enable_cache_i

     subroutine fill_cache_i (this, x)
       use core_kinds
       import base_coeffs_t
       class(base_coeffs_t), intent(inout) :: this
       real(WP), intent(in)                :: x(:)
     end subroutine fill_cache_i

  end interface

 ! Access specifiers

  private

  public :: base_coeffs_t

end module gyre_base_coeffs
