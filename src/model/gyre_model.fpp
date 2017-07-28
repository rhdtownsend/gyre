! Module   : gyre_model
! Purpose  : stellar model
!
! Copyright 2013-2017 Rich Townsend
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

  use gyre_grid

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: I_V_2 = 1
  integer, parameter :: I_As = 2
  integer, parameter :: I_U = 3
  integer, parameter :: I_C_1 = 4
  integer, parameter :: I_GAMMA_1 = 5
  integer, parameter :: I_DELTA = 6
  integer, parameter :: I_NABLA_AD = 7
  integer, parameter :: I_NABLA = 8
  integer, parameter :: I_BETA_RAD = 9
  integer, parameter :: I_C_LUM = 10
  integer, parameter :: I_C_RAD = 11
  integer, parameter :: I_C_THN = 13
  integer, parameter :: I_C_THK = 14
  integer, parameter :: I_C_EPS = 15
  integer, parameter :: I_EPS_RHO = 16
  integer, parameter :: I_EPS_T = 17
  integer, parameter :: I_KAP_RHO = 18
  integer, parameter :: I_KAP_T = 19
  integer, parameter :: I_F_LUAN_T = 20
  integer, parameter :: I_F_LUAN_C = 21
  integer, parameter :: I_OMEGA_ROT = 22

  integer, parameter :: I_LAST = I_OMEGA_ROT

  ! Derived-type definitions

  type, abstract :: model_t
   contains
     procedure(coeff), deferred      :: coeff
     procedure(coeff), deferred      :: dcoeff
     procedure(is_defined), deferred :: is_defined
     procedure(is_vacuum), deferred  :: is_vacuum
     procedure(Delta_p), deferred    :: Delta_p
     procedure(Delta_g), deferred    :: Delta_g
     procedure(grid), deferred       :: grid
  end type model_t

  ! Interfaces

  abstract interface

     function coeff (this, i, pt)
       use core_kinds
       use gyre_point
       import model_t
       class(model_t), intent(in) :: this
       integer, intent(in)        :: i
       type(point_t), intent(in)  :: pt
       real(WP)                   :: coeff
     end function coeff

     function dcoeff (this, i, pt)
       use core_kinds
       use gyre_point
       import model_t
       class(model_t), intent(in) :: this
       integer, intent(in)        :: i
       type(point_t), intent(in)  :: pt
       real(WP)                   :: dcoeff
     end function dcoeff

     function is_defined (this, i)
       import model_t
       class(model_t), intent(in) :: this
       integer, intent(in)        :: i
       logical                    :: is_defined
     end function is_defined

     function is_vacuum (this, pt)
       use gyre_point
       import model_t
       class(model_t), intent(in) :: this
       type(point_t), intent(in)  :: pt
       logical                    :: is_vacuum
     end function is_vacuum

     function Delta_p (this, x_i, x_o)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       real(WP), intent(in)       :: x_i
       real(WP), intent(in)       :: x_o
       real(WP)                   :: Delta_p
     end function Delta_p

     function Delta_g (this, x_i, x_o, lambda)
       use core_kinds
       import model_t
       class(model_t), intent(in) :: this
       real(WP), intent(in)       :: x_i
       real(WP), intent(in)       :: x_o
       real(WP), intent(in)       :: lambda
       real(WP)                   :: Delta_g
     end function Delta_g

     function grid (this)
       use gyre_grid
       import model_t
       class(model_t), intent(in) :: this
       type(grid_t)               :: grid
     end function grid

  end interface

  ! Access specifiers

  private

  public :: model_t
  public :: I_V_2
  public :: I_As
  public :: I_U
  public :: I_C_1
  public :: I_GAMMA_1
  public :: I_DELTA
  public :: I_NABLA_AD
  public :: I_NABLA
  public :: I_BETA_RAD
  public :: I_C_LUM
  public :: I_C_RAD
  public :: I_C_THN
  public :: I_C_THK
  public :: I_C_EPS
  public :: I_EPS_RHO
  public :: I_EPS_T
  public :: I_KAP_RHO
  public :: I_KAP_T
  public :: I_F_LUAN_T
  public :: I_F_LUAN_C
  public :: I_OMEGA_ROT
  public :: I_LAST

end module gyre_model
