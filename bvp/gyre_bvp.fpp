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
     procedure(init_i), deferred, public    :: init
     procedure(discrim_i), deferred, public :: discrim
     procedure(mode_i), deferred, public    :: mode
  end type bvp_t

  ! Interfaces

  abstract interface

     subroutine init_i (this, bc, op, np, shoot_gp, recon_gp, x_in, tc)
       use core_kinds
       use gyre_base_coeffs
       use gyre_therm_coeffs
       use gyre_oscpar
       use gyre_numpar
       use gyre_gridpar
       import bvp_t
       class(bvp_t), intent(out)                   :: this
       class(base_coeffs_t), intent(in)            :: bc
       type(oscpar_t), intent(in)                  :: op
       type(numpar_t), intent(in)                  :: np
       type(gridpar_t), intent(in)                 :: shoot_gp(:)
       type(gridpar_t), intent(in)                 :: recon_gp(:)
       real(WP), allocatable, intent(in)           :: x_in(:)
       class(therm_coeffs_t), intent(in), optional :: tc
     end subroutine init_i

     function discrim_i (this, omega, use_real) result (discrim)
       use core_kinds
       use gyre_ext_arith
       import bvp_t
       class(bvp_t), intent(inout)   :: this
       complex(WP), intent(in)       :: omega
       logical, intent(in), optional :: use_real
       type(ext_complex_t)           :: discrim
     end function discrim_i

     function mode_i (this, omega, discrim, use_real) result (mode)
       use core_kinds
       use gyre_mode
       use gyre_ext_arith
       import bvp_t
       class(bvp_t), intent(inout)               :: this
       complex(WP), intent(in)                   :: omega(:)
       type(ext_complex_t), intent(in), optional :: discrim(:)
       logical, intent(in), optional             :: use_real
       type(mode_t)                              :: mode
     end function mode_i

  end interface

  ! Access specifiers

  private

  public :: bvp_t

end module gyre_bvp
