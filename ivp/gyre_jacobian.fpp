! Module   : gyre_jacobian
! Purpose  : Jacobian evaluation (interface)
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

module gyre_jacobian

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: jacobian_t
     private
     integer, public :: n_e
   contains
     private
     procedure(eval_i), deferred, public :: eval
     procedure(eval_i), deferred, public :: eval_logx
  end type jacobian_t

  ! Interfaces

  abstract interface

     subroutine eval_i (this, omega, x, A)
       use core_kinds
       import jacobian_t
       class(jacobian_t), intent(in) :: this
       complex(WP), intent(in)       :: omega
       real(WP), intent(in)          :: x
       complex(WP), intent(out)      :: A(:,:)
     end subroutine eval_i

  end interface

  ! Access specifiers

  private

  public :: jacobian_t

end module gyre_jacobian
