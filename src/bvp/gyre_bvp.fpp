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
     procedure(discrim_), deferred, public :: discrim
     procedure(mode_), deferred, public    :: mode
     procedure(model_), deferred, public   :: model
  end type bvp_t

  ! Interfaces

  abstract interface

     function discrim_ (this, omega, use_real) result (discrim)
       use core_kinds
       use gyre_ext_arith
       import bvp_t
       class(bvp_t), intent(inout)   :: this
       complex(WP), intent(in)       :: omega
       logical, optional, intent(in) :: use_real
       type(ext_complex_t)           :: discrim
     end function discrim_

     function mode_ (this, omega, discrim, use_real, omega_def) result (mode)
       use core_kinds
       use gyre_mode
       use gyre_ext_arith
       import bvp_t
       class(bvp_t), target, intent(inout)       :: this
       complex(WP), intent(in)                   :: omega(:)
       type(ext_complex_t), optional, intent(in) :: discrim(:)
       logical, optional, intent(in)             :: use_real
       complex(WP), optional, intent(in)         :: omega_def(:)
       type(mode_t)                              :: mode
     end function mode_

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
