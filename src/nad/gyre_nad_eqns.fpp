! Module   : gyre_nad_eqns
! Purpose  : differential equations evaluation (nonadiabatic)
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

module gyre_nad_eqns

  ! Uses

  use core_kinds

  use gyre_eqns
  use gyre_linalg
  use gyre_model
  use gyre_nad_vars
  use gyre_osc_par
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_eqns_t) :: nad_eqns_t
     private
     class(model_t), pointer     :: ml => null()
     class(c_rot_t), allocatable :: rt
     type(nad_vars_t)            :: vr
     logical                     :: cowling_approx
     logical                     :: narf_approx
   contains
     private
     procedure, public :: A => A_
     procedure, public :: xA => xA_
  end type nad_eqns_t

  ! Interfaces

  interface nad_eqns_t
     module procedure nad_eqns_t_
  end interface nad_eqns_t

  ! Access specifiers

  private

  public :: nad_eqns_t

  ! Procedures

contains

  function nad_eqns_t_ (ml, rt, op) result (eq)

    class(model_t), pointer, intent(in)     :: ml
    class(c_rot_t), allocatable, intent(in) :: rt
    type(osc_par_t), intent(in)             :: op
    type(nad_eqns_t)                        :: eq

    ! Construct the nad_eqns_t

    eq%ml => ml
    allocate(eq%rt, SOURCE=rt)
    eq%vr = nad_vars_t(ml, rt, op)

    eq%cowling_approx = op%cowling_approx
    eq%narf_approx = op%narf_approx

    eq%n_e = 6

    ! Finish

    return

  end function nad_eqns_t_

!****

  function A_ (this, x, omega) result (A)

    class(nad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(x, omega)/x

    ! Finish

    return

  end function A_

!****

  function xA_ (this, x, omega) result (xA)

    class(nad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: xA(this%n_e,this%n_e)

    real(WP)    :: V
    real(WP)    :: V_g
    real(WP)    :: U
    real(WP)    :: As
    real(WP)    :: c_1
    real(WP)    :: nabla
    real(WP)    :: nabla_ad
    real(WP)    :: delta
    real(WP)    :: c_rad 
    real(WP)    :: dc_rad 
    real(WP)    :: c_thm
    real(WP)    :: c_dif
    real(WP)    :: c_eps_ad
    real(WP)    :: c_eps_S
    real(WP)    :: kappa_ad
    real(WP)    :: kappa_S
    complex(WP) :: lambda
    complex(WP) :: l_0
    complex(WP) :: omega_c
    real(WP)    :: alpha_gr
    real(WP)    :: alpha_hf
         
    ! Evaluate the log(x)-space RHS matrix

    ! Calculate coefficients

    V = this%ml%V_2(x)*x**2
    V_g = V/this%ml%Gamma_1(x)
    U = this%ml%U(x)
    As = this%ml%As(x)
    c_1 = this%ml%c_1(x)

    nabla = this%ml%nabla(x)
    nabla_ad = this%ml%nabla_ad(x)
    delta = this%ml%delta(x)
    c_rad = this%ml%c_rad(x)
    dc_rad = this%ml%dc_rad(x)
    c_thm = this%ml%c_thm(x)
    c_dif = this%ml%c_dif(x)
    c_eps_ad = this%ml%c_eps_ad(x)
    c_eps_S = this%ml%c_eps_S(x)
    kappa_ad = this%ml%kappa_ad(x)
    kappa_S = this%ml%kappa_S(x)

    lambda = this%rt%lambda(x, omega)
    l_0 = this%rt%l_0(omega)

    omega_c = this%rt%omega_c(x, omega)

    if (this%cowling_approx) then
       alpha_gr = 0._WP
    else
       alpha_gr = 1._WP
    endif

    if (this%narf_approx) then
       alpha_hf = 0._WP
    else
       alpha_hf = 1._WP
    endif

    ! Set up the matrix

    xA(1,1) = V_g - 1._WP - l_0
    xA(1,2) = lambda/(c_1*omega_c**2) - V_g
    xA(1,3) = alpha_gr*(V_g)
    xA(1,4) = alpha_gr*(0._WP)
    xA(1,5) = delta
    xA(1,6) = 0._WP

    xA(2,1) = c_1*omega_c**2 - As
    xA(2,2) = As - U + 3._WP - l_0
    xA(2,3) = alpha_gr*(-As)
    xA(2,4) = alpha_gr*(0._WP)
    xA(2,5) = delta
    xA(2,6) = 0._WP

    xA(3,1) = alpha_gr*(0._WP)
    xA(3,2) = alpha_gr*(0._WP)
    xA(3,3) = alpha_gr*(3._WP - U - l_0)
    xA(3,4) = alpha_gr*(1._WP)
    xA(3,5) = alpha_gr*(0._WP)
    xA(3,6) = alpha_gr*(0._WP)

    xA(4,1) = alpha_gr*(U*As)
    xA(4,2) = alpha_gr*(U*V_g)
    xA(4,3) = alpha_gr*(lambda - U*V_g)
    xA(4,4) = alpha_gr*(-U - l_0 + 2._WP)
    xA(4,5) = alpha_gr*(-U*delta)
    xA(4,6) = alpha_gr*(0._WP)

    xA(5,1) = V*(nabla_ad*(U - c_1*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif)
    xA(5,2) = V*(lambda/(c_1*omega_c**2)*(nabla_ad - nabla) - c_dif)
    xA(5,3) = alpha_gr*(V*c_dif)
    xA(5,4) = alpha_gr*(V*nabla_ad)
    xA(5,5) = V*nabla*(4._WP - kappa_S) - (l_0 - 2._WP)
    xA(5,6) = -V*nabla/c_rad

    xA(6,1) = alpha_hf*lambda*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
    xA(6,2) = V*c_eps_ad - lambda*c_rad*(alpha_hf*nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega_c**2))
    xA(6,3) = alpha_gr*(alpha_hf*lambda*nabla_ad/nabla*c_rad - V*c_eps_ad)
    xA(6,4) = alpha_gr*(0._WP)
    if (x > 0._WP) then
       xA(6,5) = c_eps_S - alpha_hf*lambda*c_rad/(nabla*V) + (0._WP,1._WP)*omega_c*c_thm
    else
       xA(6,5) = -alpha_hf*HUGE(0._WP)
    endif
    xA(6,6) = -1._WP - l_0

    ! Apply the variables transformation

    xA = MATMUL(this%vr%S(x, omega), MATMUL(xA, this%vr%T(x, omega)) - this%vr%dT(x, omega))

    ! Finish

    return

  end function xA_

end module gyre_nad_eqns
