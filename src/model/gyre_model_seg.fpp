! Module   : gyre_model_seg
! Purpose  : stellar model segment (interface)
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
$include 'core_parallel.inc'

module gyre_model_seg

  ! Uses

  use core_kinds

  use gyre_seg

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure(y_1_), deferred :: ${NAME}_1_
    procedure(y_v_), deferred :: ${NAME}_v_
    generic, public           :: ${NAME} => ${NAME}_1_, ${NAME}_v_
  $endsub

  type, abstract, extends(seg_t) :: model_seg_t
   contains
     private
     $PROC_DECL(V_2)
     $PROC_DECL(As)
     $PROC_DECL(U)
     $PROC_DECL(dU)
     $PROC_DECL(c_1)
     $PROC_DECL(Gamma_1)
     $PROC_DECL(delta)
     $PROC_DECL(nabla_ad)
     $PROC_DECL(dnabla_ad)
     $PROC_DECL(nabla)
     $PROC_DECL(c_rad)
     $PROC_DECL(dc_rad)
     $PROC_DECL(c_thm)
     $PROC_DECL(c_dif)
     $PROC_DECL(c_eps_ad)
     $PROC_DECL(c_eps_S)
     $PROC_DECL(kappa_ad)
     $PROC_DECL(kappa_S)
     $PROC_DECL(Omega_rot)
     $PROC_DECL(dOmega_rot)
     procedure(scaffold_), deferred, public :: scaffold
  end type model_seg_t

  ! Interfaces

  abstract interface

     function y_1_(this, x) result(y)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       real(WP), intent(in)       :: x
       real(WP)                   :: y
     end function y_1_

     function y_v_(this, x) result(y)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       real(WP), intent(in)       :: x(:)
       real(WP)                   :: y(SIZE(x))
     end function y_v_

     subroutine scaffold_(this, x)
       use core_kinds
       import model_t
       class(model_t), intent(in)         :: this
       real(WP), allocatable, intent(out) :: x(:)
     end subroutine scaffold_

  end interface

  ! Access specifiers

  private

  public :: model_seg_t

end module gyre_model_seg

