! Module   : gyre_wos_eqns
! Purpose  : differential equations evaluation (adiabatic)
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

module gyre_wos_eqns

  ! Uses

  use core_kinds

  use gyre_eqns

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  ! Derived-type definitions

  type, extends (r_eqns_t) :: wos_eqns_t
     private
     real(WP) :: c
   contains
     private
     procedure, public :: A => A_
     procedure, public :: xA => xA_
     procedure, public :: T => T_
  end type wos_eqns_t

  ! Interfaces

  interface wos_eqns_t
     module procedure wos_eqns_t_
  end interface wos_eqns_t

  ! Access specifiers

  private

  public :: wos_eqns_t

  ! Procedures

contains

  function wos_eqns_t_ (c) result (eq)

    real(WP), intent(in)   :: c
    type(wos_eqns_t)       :: eq

    ! Construct the wos_eqns_t

    eq%n_e = 2
    eq%c = c

    ! Finish

    return

  end function wos_eqns_t_

!****

  function A_ (this, x, omega) result (A)

    class(wos_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A(1,1) = 0._WP
    A(1,2) = 1._WP

    A(2,1) = 0.5*x**2 - omega**2 / this%c**2
    A(2,2) = 0._WP
		
    ! Finish

    return

  end function A_

!****

  function xA_ (this, x, omega) result (xA)

    class(wos_eqns_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: xA(this%n_e,this%n_e)
    
    ! Evaluate the log(x)-space RHS matrix (=x*A)

    xA = x * this%A(x, omega)
		
    ! Finish

    return

  end function xA_

!****

  function T_ (this, x, omega, to_canon) result (T)

    class(wos_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    logical, intent(in)           :: to_canon
    real(WP)                      :: T(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables to/from
    ! the canonical formulation

    T(1,1) = 1._WP
    T(1,2) = 0._WP

    T(2,1) = 0._WP
    T(2,2) = 1._WP
		
    ! Finish

    return

  end function T_

end module gyre_wos_eqns
