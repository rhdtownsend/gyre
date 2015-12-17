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

  type :: model_t
     private
     integer, public :: n_s
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
     procedure(scaffold_), deferred, public :: scaffold
  end type model_t

  ! Interfaces

  abstract interface

     function f_1_ (this, mc) result (f)
       use core_kinds
       use gyre_coords
       import model_t
       class(model_t), intent(in)  :: this
       class(coords_t), intent(in) :: co
       real(WP)                    :: f
     end function f_1_

     function f_v_ (this, mc) result (f)
       use core_kinds
       use gyre_coords
       import model_t
       class(model_t), intent(in) :: this
       type(coords_t), intent(in) :: co
       real(WP)                   :: f(SIZE(co))
     end function f_v_

     function scaffold_ (this) result (co)
       use gyre_coords
       import model_t
       class(model_t), intent(in)   :: this
       class(coords_t), allocatable :: co(:)
     end function scaffold_

  end interface

  ! Access specifiers

  private

  public :: model_t

end module gyre_model
