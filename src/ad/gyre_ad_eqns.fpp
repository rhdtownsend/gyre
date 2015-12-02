! Module   : gyre_ad_eqns
! Purpose  : differential equations for adibatic case
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

module gyre_ad_eqns

  ! Uses

  use core_kinds

  use gyre_ad_vars
  use gyre_eqns
  use gyre_model_seg
  use gyre_mode_par
  use gyre_osc_par
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (r_eqns_t) :: ad_eqns_t
     private
     class(model_seg_t), pointer :: ms => null()
     class(r_rot_t), allocatable :: rt
     type(ad_vars_t)             :: vr
     logical                     :: cowling_approx
   contains
     private
     procedure, public :: A => A_
     procedure, public :: xA => xA_
  end type ad_eqns_t

  ! Interfaces

  interface ad_eqns_t
     module procedure ad_eqns_t_
  end interface ad_eqns_t

  ! Access specifiers

  private

  public :: ad_eqns_t

  ! Procedures

contains

  function ad_eqns_t_ (ms, md_p, os_p) result (eq)

    class(model_seg_t), pointer, intent(in) :: ms
    type(mode_par_t), intent(in)            :: md_p
    type(osc_par_t), intent(in)             :: os_p
    type(ad_eqns_t)                         :: eq

    ! Construct the ad_eqns_t

    eq%ms => ms

    allocate(eq%rt, SOURCE=r_rot_t(ms, md_p, os_p))
    eq%vr = ad_vars_t(ms, md_p, os_p)

    eq%cowling_approx = os_p%cowling_approx

    eq%n_e = 4

    ! Finish

    return

  end function ad_eqns_t_

!****

  function A_ (this, x, omega) result (A)

    class(ad_eqns_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(x, omega)/x

    ! Finish

    return

  end function A_

!****

  function xA_ (this, x, omega) result (xA)

    class(ad_eqns_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: xA(this%n_e,this%n_e)

    real(WP) :: V_g
    real(WP) :: U
    real(WP) :: As
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: l_0
    real(WP) :: omega_c
    real(WP) :: alpha_gr
    
    ! Evaluate the log(x)-space RHS matrix

    ! Calculate coefficients

    V_g = this%ms%V_2(x)*x**2/this%ms%Gamma_1(x)
    U = this%ms%U(x)
    As = this%ms%As(x)
    c_1 = this%ms%c_1(x)

    lambda = this%rt%lambda(x, omega)
    l_0 = this%rt%l_0(0._WP, omega) ! THIS NEEDS FIXING; ORIGIN IS NOT IN MODEL SEGMENT!

    omega_c = this%rt%omega_c(x, omega)

    if (this%cowling_approx) then
       alpha_gr = 0._WP
    else
       alpha_gr = 1._WP
    endif

    ! Set up the matrix

    xA(1,1) = V_g - 1._WP - l_0
    xA(1,2) = lambda/(c_1*omega_c**2) - V_g
    xA(1,3) = alpha_gr*(V_g)
    xA(1,4) = alpha_gr*(0._WP)
      
    xA(2,1) = c_1*omega_c**2 - As
    xA(2,2) = As - U + 3._WP - l_0
    xA(2,3) = alpha_gr*(-As)
    xA(2,4) = alpha_gr*(0._WP)
      
    xA(3,1) = alpha_gr*(0._WP)
    xA(3,2) = alpha_gr*(0._WP)
    xA(3,3) = alpha_gr*(3._WP - U - l_0)
    xA(3,4) = alpha_gr*(1._WP)
      
    xA(4,1) = alpha_gr*(U*As)
    xA(4,2) = alpha_gr*(U*V_g)
    xA(4,3) = alpha_gr*(lambda - U*V_g)
    xA(4,4) = alpha_gr*(-U - l_0 + 2._WP)

    ! Apply the variables transformation

    xA = MATMUL(this%vr%S(x, omega), MATMUL(xA, this%vr%T(x, omega)) - this%vr%dT(x, omega))

    ! Finish

    return

  end function xA_
  
end module gyre_ad_eqns
