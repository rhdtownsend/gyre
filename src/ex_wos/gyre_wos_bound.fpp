! Incfile  : gyre_wos_bound
! Purpose  : boundary conditions (adiabatic)
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

module gyre_wos_bound

  ! Uses

  use core_kinds

  use gyre_bound
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  ! Derived-type definitions

  type, extends (r_bound_t) :: wos_bound_t
     private
     real(WP) :: x_i
     real(WP) :: x_o
   contains 
     private
     procedure, public :: B_i => B_i_
     procedure, public :: B_o => B_o_
  end type wos_bound_t

  ! Interfaces

  interface wos_bound_t
     module procedure wos_bound_t_
  end interface wos_bound_t

  ! Access specifiers

  private

  public :: wos_bound_t

  ! Procedures

contains

  function wos_bound_t_ (x_i, x_o) result (bd)

    real(WP), intent(in) :: x_i
    real(WP), intent(in) :: x_o
    type(wos_bound_t)    :: bd

    ! Construct the wos_bound_t
    
    bd%x_i = x_i
    bd%x_o = x_o

    bd%n_i = 1
    bd%n_o = 1
    bd%n_e = bd%n_i + bd%n_o 

    ! Finish

    return
    
  end function wos_bound_t_

!****

  function B_i_ (this, omega) result (B_i)

    class(wos_bound_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP)                       :: B_i(this%n_i,this%n_e)

    ! Evaluate the inner boundary conditions

    B_i(1,1) = 1._WP
    B_i(1,2) = 0._WP

    ! Finish

    return

  end function B_i_

!****

  function B_o_ (this, omega) result (B_o)

    class(wos_bound_t), intent(in) :: this
    real(WP), intent(in)           :: omega
    real(WP)                       :: B_o(this%n_o,this%n_e)

    ! Evaluate the outer boundary conditions

    B_o(1,1) = 1._WP
    B_o(1,2) = 0._WP
    
    ! Finish

    return

  end function B_o_

end module gyre_wos_bound
