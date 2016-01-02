! Module   : gyre_model
! Purpose  : stellar model
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

module gyre_model

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  $define $PROC_DECL $sub
    $local $NAME $1
    procedure(f_1_), deferred :: ${NAME}_1_
    procedure(f_v_), deferred :: ${NAME}_v_
    generic, public           :: ${NAME} => ${NAME}_1_, ${NAME}_v_
  $endsub

  type, abstract :: model_t
     private
     integer, public  :: n_s
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
     $PROC_DECL(beta_rad)
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
     procedure(x_i), deferred, public :: x_i
     procedure(x_o), deferred, public :: x_o
     procedure(x), deferred, public   :: x
  end type model_t

  ! Interfaces

  abstract interface

     function f_1_ (this, s, x) result (f)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       integer, intent(in)        :: s
       real(WP), intent(in)       :: x
       real(WP)                   :: f
     end function f_1_

     function f_v_ (this, s, x) result (f)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       integer, intent(in)        :: s(:)
       real(WP), intent(in)       :: x(:)
       real(WP)                   :: f(SIZE(s))
     end function f_v_

     function x_i (this, s)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       integer, intent(in)        :: s
       real(WP)                   :: x_i
     end function x_i

     function x_o (this, s)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       integer, intent(in)        :: s
       real(WP)                   :: x_o
     end function x_o

     function x (this, s)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       integer, intent(in)        :: s
       real(WP), allocatable      :: x(:)
     end function x

  end interface

  ! Access specifiers

  private

  public :: model_t

end module gyre_model
