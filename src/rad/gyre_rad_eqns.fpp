! Module   : gyre_rad_eqns
! Purpose  : radial adiabatic differential equations
!
! Copyright 2013-2016 Rich Townsend
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

module gyre_rad_eqns

  ! Uses

  use core_kinds

  use gyre_eqns
  use gyre_model
  use gyre_mode_par
  use gyre_osc_par
  use gyre_rad_vars
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_eqns_t) :: rad_eqns_t
     private
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     type(rad_vars_t)            :: vr
     integer                     :: s
   contains
     private
     procedure, public :: A => A_
     procedure, public :: xA => xA_
  end type rad_eqns_t

  ! Interfaces

  interface rad_eqns_t
     module procedure rad_eqns_t_
  end interface rad_eqns_t

  ! Access specifiers

  private

  public :: rad_eqns_t

  ! Procedures

contains

  function rad_eqns_t_ (ml, s, md_p, os_p) result (eq)

    class(model_t), pointer, intent(in) :: ml
    integer, intent(in)                 :: s
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(rad_eqns_t)                    :: eq

    ! Construct the rad_eqns_t

    eq%ml => ml

    allocate(eq%rt, SOURCE=r_rot_t(ml, md_p, os_p))
    eq%vr = rad_vars_t(ml, md_p, os_p)

    eq%s = s

    eq%n_e = 2

    ! Finish

    return

  end function rad_eqns_t_

  !****

  function A_ (this, x, omega) result (A)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(x, omega)/x

    ! Finish

    return

  end function A_

!****

  function xA_ (this, x, omega) result (xA)

    class(rad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: xA(this%n_e,this%n_e)

    real(WP) :: V_g
    real(WP) :: U
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: omega_c
    
    ! Evaluate the log(x)-space RHS matrix

    associate (s => this%s)

      ! Calculate coefficients

      V_g = this%ml%V_2(s, x)*x**2/this%ml%Gamma_1(s, x)
      U = this%ml%U(s, x)
      As = this%ml%As(s, x)
      c_1 = this%ml%c_1(s, x)

      omega_c = this%rt%omega_c(s, x, omega)

      ! Set up the matrix

      xA(1,1) = V_g - 1._WP
      xA(1,2) = -V_g
      
      xA(2,1) = c_1*omega_c**2 + U - As
      xA(2,2) = As - U + 3._WP

      ! Apply the variables transformation

      xA = MATMUL(this%vr%G(s, x, omega), MATMUL(xA, this%vr%H(s, x, omega)) -&
           this%vr%dH(s, x, omega))

    end associate

    ! Finish

    return

  end function xA_

end module gyre_rad_eqns
