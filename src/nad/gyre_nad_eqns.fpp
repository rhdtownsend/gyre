! Module   : gyre_nad_eqns
! Purpose  : nonadiabatic differential equations
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

module gyre_nad_eqns

  ! Uses

  use core_kinds

  use gyre_eqns
  use gyre_linalg
  use gyre_mode_par
  use gyre_model
  use gyre_nad_vars
  use gyre_osc_par
  use gyre_rot
  use gyre_rot_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_eqns_t) :: nad_eqns_t
     private
     class(model_t), pointer     :: ml => null()
     class(c_rot_t), allocatable :: rt
     type(nad_vars_t)            :: vr
     integer                     :: s
     logical                     :: cowling_approx
     logical                     :: narf_approx
   contains
     private
     procedure, public :: A
     procedure, public :: xA
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

  function nad_eqns_t_ (ml, s, md_p, os_p) result (eq)

    class(model_t), pointer, intent(in) :: ml
    integer, intent(in)                 :: s
    type(mode_par_t), intent(in)        :: md_p
    type(osc_par_t), intent(in)         :: os_p
    type(nad_eqns_t)                    :: eq

    ! Construct the nad_eqns_t

    eq%ml => ml

    allocate(eq%rt, SOURCE=c_rot_t(ml, md_p, os_p))
    eq%vr = nad_vars_t(ml, md_p, os_p)

    eq%s = s

    eq%cowling_approx = os_p%cowling_approx
    eq%narf_approx = os_p%narf_approx

    eq%n_e = 6

    ! Finish

    return

  end function nad_eqns_t_

  !****

  function A (this, x, omega)

    class(nad_eqns_t), intent(in) :: this
    real(WP), intent(in)          :: x
    complex(WP), intent(in)       :: omega
    complex(WP)                   :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(x, omega)/x

    ! Finish

    return

  end function A

  !****

  function xA (this, x, omega)

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
    real(WP)    :: dnabla_ad
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
    complex(WP) :: l_i
    complex(WP) :: omega_c
    real(WP)    :: alpha_gr
    real(WP)    :: alpha_hf
         
    ! Evaluate the log(x)-space RHS matrix

    associate (s => this%s)

      ! Calculate coefficients

      V = this%ml%V_2(s, x)*x**2
      V_g = V/this%ml%Gamma_1(s, x)
      U = this%ml%U(s, x)
      As = this%ml%As(s, x)
      c_1 = this%ml%c_1(s, x)

      nabla = this%ml%nabla(s, x)
      nabla_ad = this%ml%nabla_ad(s, x)
      dnabla_ad = this%ml%dnabla_ad(s, x)
      delta = this%ml%delta(s, x)
      c_rad = this%ml%c_rad(s, x)
      dc_rad = this%ml%dc_rad(s, x)
      c_thm = this%ml%c_thm(s, x)
      c_dif = this%ml%c_dif(s, x)
      c_eps_ad = this%ml%c_eps_ad(s, x)
      c_eps_S = this%ml%c_eps_S(s, x)
      kappa_ad = this%ml%kappa_ad(s, x)
      kappa_S = this%ml%kappa_S(s, x)

      lambda = this%rt%lambda(s, x, omega)
      l_i = this%rt%l_i(omega)

      omega_c = this%rt%omega_c(s, x, omega)

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

      xA(1,1) = V_g - 1._WP - l_i
      xA(1,2) = lambda/(c_1*omega_c**2) - V_g
      xA(1,3) = alpha_gr*(V_g)
      xA(1,4) = alpha_gr*(0._WP)
      xA(1,5) = delta
      xA(1,6) = 0._WP

      xA(2,1) = c_1*omega_c**2 - As
      xA(2,2) = As - U + 3._WP - l_i
      xA(2,3) = alpha_gr*(-As)
      xA(2,4) = alpha_gr*(0._WP)
      xA(2,5) = delta
      xA(2,6) = 0._WP

      xA(3,1) = alpha_gr*(0._WP)
      xA(3,2) = alpha_gr*(0._WP)
      xA(3,3) = alpha_gr*(3._WP - U - l_i)
      xA(3,4) = alpha_gr*(1._WP)
      xA(3,5) = alpha_gr*(0._WP)
      xA(3,6) = alpha_gr*(0._WP)

      xA(4,1) = alpha_gr*(U*As)
      xA(4,2) = alpha_gr*(U*V_g)
      xA(4,3) = alpha_gr*(lambda - U*V_g)
      xA(4,4) = alpha_gr*(-U - l_i + 2._WP)
      xA(4,5) = alpha_gr*(-U*delta)
      xA(4,6) = alpha_gr*(0._WP)

      xA(5,1) = V*(nabla_ad*(U - c_1*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif + nabla_ad*dnabla_ad)
      xA(5,2) = V*(lambda/(c_1*omega_c**2)*(nabla_ad - nabla) - (c_dif + nabla_ad*dnabla_ad))
      xA(5,3) = alpha_gr*(V*(c_dif + nabla_ad*dnabla_ad))
      xA(5,4) = alpha_gr*(V*nabla_ad)
      xA(5,5) = V*nabla*(4._WP - kappa_S) - (l_i - 2._WP)
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
      xA(6,6) = -1._WP - l_i

      ! Apply the variables transformation

      xA = MATMUL(this%vr%G(s, x, omega), MATMUL(xA, this%vr%H(s, x, omega)) - &
                                          this%vr%dH(s, x, omega))

    end associate

    ! Finish

    return

  end function xA

end module gyre_nad_eqns
