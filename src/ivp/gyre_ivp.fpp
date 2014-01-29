! Module   : gyre_ivp
! Purpose  : solve initial-value problems (interface)
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

module gyre_ivp

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ivp_t
     private
     integer, public :: n_e
   contains
     private
     procedure(solve_), deferred, public    :: solve
     procedure(recon_), deferred, public    :: recon
     procedure(abscissa_), deferred, public :: abscissa
  end type ivp_t

  ! Interfaces

  abstract interface

     subroutine solve_ (this, omega, x_a, x_b, E_l, E_r, S, use_real)
       use core_kinds
       use gyre_jacobian
       use gyre_ext_arith
       import ivp_t
       class(ivp_t), intent(in)         :: this
       complex(WP), intent(in)          :: omega
       real(WP), intent(in)             :: x_a
       real(WP), intent(in)             :: x_b
       complex(WP), intent(out)         :: E_l(:,:)
       complex(WP), intent(out)         :: E_r(:,:)
       type(ext_complex_t), intent(out) :: S
       logical, optional, intent(in)    :: use_real
     end subroutine solve_

     subroutine recon_ (this, omega, x_a, x_b, y_a, y_b, x, y, use_real)
       use core_kinds
       use gyre_jacobian
       import ivp_t
       class(ivp_t), intent(in)      :: this
       complex(WP), intent(in)       :: omega
       real(WP), intent(in)          :: x_a
       real(WP), intent(in)          :: x_b
       complex(WP), intent(in)       :: y_a(:)
       complex(WP), intent(in)       :: y_b(:)
       real(WP), intent(in)          :: x(:)
       complex(WP), intent(out)      :: y(:,:)
       logical, optional, intent(in) :: use_real
     end subroutine recon_

     function abscissa_ (this, x_a, x_b) result (x)
       use core_kinds
       import ivp_t
       class(ivp_t), intent(in) :: this
       real(WP), intent(in)     :: x_a
       real(WP), intent(in)     :: x_b
       real(WP), allocatable    :: x(:)
     end function abscissa_

  end interface

  ! Access specifiers

  private

  public :: ivp_t

end module gyre_ivp
