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

  use core_kinds

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: bvp_t
   contains 
     private
     procedure(discrim_i), deferred, public :: discrim
     procedure(recon_i), deferred, public   :: recon
  end type bvp_t

  ! Interfaces

  abstract interface

     function discrim_i (this, omega, norm) result (discrim)
       use core_kinds
       use gyre_ext_arith
       import bvp_t
       class(bvp_t), intent(inout)   :: this
       complex(WP), intent(in)       :: omega
       logical, intent(in), optional :: norm
       type(ext_complex_t)           :: discrim
     end function discrim_i

     subroutine recon_i (this, omega, x, y)
       use core_kinds
       import bvp_t
       class(bvp_t), intent(inout)           :: this
       complex(WP), intent(in)               :: omega
       real(WP), allocatable, intent(out)    :: x(:)
       complex(WP), allocatable, intent(out) :: y(:,:)
     end subroutine recon_i
       
  end interface

  ! Access specifiers

  private

  public :: bvp_t

end module gyre_bvp
