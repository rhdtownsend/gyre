! Module   : gyre_therm_coeffs
! Purpose  : thermal structure coefficients (interface)
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

module gyre_therm_coeffs

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure(get_1_i), deferred :: get_${NAME}_1
    procedure(get_v_i), deferred :: get_${NAME}_v
    generic, public              :: ${NAME} => get_${NAME}_1, get_${NAME}_v
  $endsub

  type, abstract :: therm_coeffs_t
     private
   contains
     private
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_gen)
     $PROC_DECL(c_thm)
     $PROC_DECL(c_dif)
     $PROC_DECL(nabla)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(epsilon_ad)
     $PROC_DECL(epsilon_S)
  end type therm_coeffs_t

  ! Interfaces

  abstract interface

     function get_1_i (this, x) result (y)
       use core_kinds
       import therm_coeffs_t
       class(therm_coeffs_t), intent(in) :: this
       real(WP), intent(in)              :: x
       real(WP)                          :: y
     end function get_1_i

     function get_v_i (this, x) result (y)
       use core_kinds
       import therm_coeffs_t
       class(therm_coeffs_t), intent(in) :: this
       real(WP), intent(in)              :: x(:)
       real(WP)                          :: y(SIZE(x))
     end function get_v_i

  end interface

 ! Access specifiers

  private

  public :: therm_coeffs_t

end module gyre_therm_coeffs
