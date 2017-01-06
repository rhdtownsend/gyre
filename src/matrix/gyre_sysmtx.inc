! Incfile  : gyre_sysmtx
! Purpose  : system matrix (template)
!
! Copyright 2013-2015 Rich Townsend
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

$define $MOD_DECL $sub

$local $T $1
$local $TYPE $2

module gyre_${T}_sysmtx

  ! Uses

  use core_kinds

  use gyre_ext

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: ${T}_sysmtx_t
     private
     integer, public :: n   ! Number of equation blocks
     integer, public :: n_e ! Number of equations per block
     integer, public :: n_i ! Number of inner boundary conditions
     integer, public :: n_o ! Number of outer boundary conditions
   contains
     private
     procedure(set_B_), deferred, public     :: set_B_i
     procedure(set_B_), deferred, public     :: set_B_o
     procedure(set_E_), deferred, public     :: set_E
     procedure(factorize_), deferred, public :: factorize
     procedure(det_), deferred, public       :: det
     procedure(soln_vec_hom_), deferred      :: soln_vec_hom_
     procedure(soln_vec_inhom_), deferred    :: soln_vec_inhom_
     generic, public                         :: soln_vec => soln_vec_hom_, soln_vec_inhom_
  end type ${T}_sysmtx_t

  ! Interfaces

  abstract interface

     subroutine set_B_ (this, B, scl)
       use core_kinds
       use gyre_ext
       import ${T}_sysmtx_t
       class(${T}_sysmtx_t), intent(inout) :: this
       $TYPE(WP), intent(in)               :: B(:,:)
       $TYPE(WP), intent(in)               :: scl(:)
     end subroutine set_B_
       
     subroutine set_E_ (this, k, E_l, E_r, scl)
       use core_kinds
       use gyre_ext
       import ${T}_sysmtx_t
       class(${T}_sysmtx_t), intent(inout) :: this
       integer, intent(in)                 :: k
       $TYPE(WP), intent(in)               :: E_l(:,:)
       $TYPE(WP), intent(in)               :: E_r(:,:)
       type(${T}_ext_t), intent(in)        :: scl
     end subroutine set_E_

     subroutine factorize_ (this)
       import ${T}_sysmtx_t
       class(${T}_sysmtx_t), intent(inout) :: this
     end subroutine factorize_

     function det_ (this) result (det)
       use gyre_ext
       import ${T}_sysmtx_t
       class(${T}_sysmtx_t), intent(in) :: this
       type(${T}_ext_t)                 :: det
     end function det_

     function soln_vec_hom_ (this) result (v)
       use core_kinds
       use gyre_ext
       import ${T}_sysmtx_t
       class(${T}_sysmtx_t), intent(in) :: this
       $TYPE(WP)                        :: v(this%n_e*(this%n+1))
     end function soln_vec_hom_

     function soln_vec_inhom_ (this, w_i, w_o) result (v)
       use core_kinds
       use gyre_ext
       import ${T}_sysmtx_t
       class(${T}_sysmtx_t), intent(in) :: this
       $TYPE(WP), intent(in)            :: w_i(:)
       $TYPE(WP), intent(in)            :: w_o(:)
       $TYPE(WP)                        :: v(this%n_e*(this%n+1))
     end function soln_vec_inhom_

  end interface

  ! Access specifiers

  private

  public :: ${T}_sysmtx_t

end module gyre_${T}_sysmtx

$endsub
