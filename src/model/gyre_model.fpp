! Module   : gyre_model
! Purpose  : stellar model (interface)
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

  type, abstract :: model_t
     integer :: n_sg
   contains
     private
     procedure(sg_), deferred, public :: sg
  end type model_t

  ! Interfaces

  abstract interface

     function sg_ (this) result (sg)
       use gyre_model_seg
       import model_t
       class(model_t), intent(in)  :: this
       class(model_seg_t), pointer :: sg
     end function sg_

  end interface

 ! Access specifiers

  private

  public :: model_t

end module gyre_model
