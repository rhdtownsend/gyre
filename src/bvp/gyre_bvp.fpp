! Module   : gyre_bvp
! Purpose  : solve boundary-value problems (interface)
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

module gyre_bvp

  ! Uses

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: bvp_t
   contains 
     private
     procedure(discrim_r_), deferred       :: discrim_r_
     procedure(discrim_c_), deferred       :: discrim_c_
     generic, public                       :: discrim => discrim_r_, discrim_c_
     procedure(mode_r_), deferred          :: mode_r_
     procedure(mode_c_), deferred          :: mode_c_
     generic, public                       :: mode => mode_r_, mode_c_
     procedure(model_), deferred, public   :: model
  end type bvp_t

  ! Interfaces

  abstract interface

     function discrim_r_ (this, omega) result (discrim)
       use core_kinds
       use gyre_ext_arith
       import bvp_t
       class(bvp_t), intent(inout) :: this
       real(WP), intent(in)        :: omega
       type(ext_real_t)            :: discrim
     end function discrim_r_

     function discrim_c_ (this, omega) result (discrim)
       use core_kinds
       use gyre_ext_arith
       import bvp_t
       class(bvp_t), intent(inout)   :: this
       complex(WP), intent(in)       :: omega
       type(ext_complex_t)           :: discrim
     end function discrim_c_

     function mode_r_ (this, omega) result (mode)
       use core_kinds
       use gyre_mode
       import bvp_t
       class(bvp_t), target, intent(inout) :: this
       real(WP), intent(in)                :: omega
       type(mode_t)                        :: mode
     end function mode_r_

     function mode_c_ (this, omega) result (mode)
       use core_kinds
       use gyre_mode
       import bvp_t
       class(bvp_t), target, intent(inout) :: this
       complex(WP), intent(in)             :: omega
       type(mode_t)                        :: mode
     end function mode_c_

     function model_ (this) result (ml)
       use gyre_model
       import bvp_t
       class(bvp_t), intent(in) :: this
       class(model_t), pointer  :: ml
     end function model_

  end interface

  ! Access specifiers

  private

  public :: bvp_t

end module gyre_bvp
