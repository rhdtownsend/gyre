! Module   : gyre_nad_jacob
! Purpose  : jacobian evaluation (nonadiabatic)
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

module gyre_nad_jacob

  ! Uses

  use core_kinds

  use gyre_jacob
  use gyre_linalg
  use gyre_model
  use gyre_osc_par
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: DZIEM_VARS = 1
  integer, parameter :: JCD_VARS = 2
  integer, parameter :: LAGP_VARS = 3

  ! Derived-type definitions

  type, extends (c_jacob_t) :: nad_jacob_t
     private
     class(model_t), pointer     :: ml => null()
     class(c_rot_t), allocatable :: rt
     integer                     :: vars
   contains
     private
     procedure, public :: A => A_
     procedure, public :: xA => xA_
     procedure         :: xA_dziem_
     procedure         :: xA_jcd_
     procedure         :: xA_lagp_
     procedure, public :: T => T_
     procedure         :: T_jcd_
     procedure         :: T_lagp_
  end type nad_jacob_t

  ! Interfaces

  interface nad_jacob_t
     module procedure nad_jacob_t_
  end interface nad_jacob_t

  ! Access specifiers

  private

  public :: nad_jacob_t

  ! Procedures

contains

  function nad_jacob_t_ (ml, rt, op) result (jc)

    class(model_t), pointer, intent(in)     :: ml
    class(c_rot_t), allocatable, intent(in) :: rt
    type(osc_par_t), intent(in)             :: op
    type(nad_jacob_t)                       :: jc

    ! Construct the nad_jacob_t

    jc%ml => ml
    allocate(jc%rt, SOURCE=rt)

    select case (op%variables_set)
    case ('DZIEM')
       jc%vars = DZIEM_VARS
    case ('JCD')
       jc%vars = JCD_VARS
    case ('LAGP')
       jc%vars = LAGP_VARS
    case default
       $ABORT(Invalid variables_set)
    end select

    jc%n_e = 6

    ! Finish

    return

  end function nad_jacob_t_

!****

  function A_ (this, x, omega) result (A)

    class(nad_jacob_t), intent(in) :: this
    real(WP), intent(in)           :: x
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: A(this%n_e,this%n_e)
    
    ! Evaluate the RHS matrix

    A = this%xA(x, omega)/x

    ! Finish

    return

  end function A_

!****

  function xA_ (this, x, omega) result (xA)

    class(nad_jacob_t), intent(in) :: this
    real(WP), intent(in)           :: x
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: xA(this%n_e,this%n_e)
    
    ! Evaluate the log(x)-space RHS matrix (=x*A)

    select case (this%vars)
    case (DZIEM_VARS)
       xA = this%xA_dziem_(x, omega)
    case (JCD_VARS)
       xA = this%xA_jcd_(x, omega)
    case (LAGP_VARS)
       xA = this%xA_lagp_(x, omega)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function xA_

!****

  function xA_dziem_ (this, x, omega) result (xA)

    class(nad_jacob_t), intent(in) :: this
    real(WP), intent(in)           :: x
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: xA(this%n_e,this%n_e)

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
         
    ! Evaluate the log(x)-space RHS matrix ([Dzi1971] formulation)

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

    ! Set up the matrix

    xA(1,1) = V_g - 1._WP - l_0
    xA(1,2) = lambda/(c_1*omega_c**2) - V_g
    xA(1,3) = V_g
    xA(1,4) = 0._WP
    xA(1,5) = delta
    xA(1,6) = 0._WP

    xA(2,1) = c_1*omega_c**2 - As
    xA(2,2) = As - U + 3._WP - l_0
    xA(2,3) = -As
    xA(2,4) = 0._WP
    xA(2,5) = delta
    xA(2,6) = 0._WP

    xA(3,1) = 0._WP
    xA(3,2) = 0._WP
    xA(3,3) = 3._WP - U - l_0
    xA(3,4) = 1._WP
    xA(3,5) = 0._WP
    xA(3,6) = 0._WP

    xA(4,1) = U*As
    xA(4,2) = U*V_g
    xA(4,3) = lambda - U*V_g
    xA(4,4) = -U - l_0 + 2._WP
    xA(4,5) = -U*delta
    xA(4,6) = 0._WP

    xA(5,1) = V*(nabla_ad*(U - c_1*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif)
    xA(5,2) = V*(lambda/(c_1*omega_c**2)*(nabla_ad - nabla) - c_dif)
    xA(5,3) = V*c_dif
    xA(5,4) = V*nabla_ad
    xA(5,5) = V*nabla*(4._WP - kappa_S) - (l_0 - 2._WP)
    xA(5,6) = -V*nabla/c_rad

    xA(6,1) = lambda*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
    xA(6,2) = V*c_eps_ad - lambda*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega_c**2))
    xA(6,3) = lambda*nabla_ad/nabla*c_rad - V*c_eps_ad
    xA(6,4) = 0._WP
    if (x > 0._WP) then
       xA(6,5) = c_eps_S - lambda*c_rad/(nabla*V) - (0._WP,1._WP)*omega_c*c_thm
    else
       xA(6,5) = -HUGE(0._WP)
    endif
    xA(6,6) = -1._WP - l_0

    ! Finish

    return

  end function xA_dziem_

!****

  function xA_jcd_ (this, x, omega) result (xA)

    class(nad_jacob_t), intent(in) :: this
    real(WP), intent(in)           :: x
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: xA(this%n_e,this%n_e)
    
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
    integer     :: l
    complex(WP) :: lambda
    complex(WP) :: l_0
    complex(WP) :: omega_c
         
    ! Evaluate the log(x)-space RHS matrix ([Chr2008] formulation)

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

    l = this%rt%mp%l
    lambda = this%rt%lambda(x, omega)
    l_0 = this%rt%l_0(omega)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix

    if (l /= 0) then

       xA(1,1) = V_g - 1._WP - l_0
       xA(1,2) = 1._WP - V_g*c_1*omega_c**2/lambda
       xA(1,3) = -V_g
       xA(1,4) = 0._WP
       xA(1,5) = delta
       xA(1,6) = 0._WP
      
       xA(2,1) = lambda - As*lambda/(c_1*omega_c**2)
       xA(2,2) = As - l_0
       xA(2,3) = As*lambda/(c_1*omega_c**2)
       xA(2,4) = 0._WP
       xA(2,5) = delta*lambda/(c_1*omega_c**2)
       xA(2,6) = 0._WP
      
       xA(3,1) = 0._WP
       xA(3,2) = 0._WP
       xA(3,3) = 2._WP - l_0
       xA(3,4) = 1._WP
       xA(3,5) = 0._WP
       xA(3,6) = 0._WP
      
       xA(4,1) = -U*As
       xA(4,2) = -U*V_g*c_1*omega_c**2/lambda
       xA(4,3) = lambda + U*(As - 2._WP)
       xA(4,4) = 2._WP*(1._WP-U) - (l_0 - 1._WP)
       xA(4,5) = U*delta
       xA(4,6) = 0._WP

       xA(5,1) = V*(nabla_ad*(U - c_1*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif)
       xA(5,2) = V*(lambda/(c_1*omega_c**2)*(nabla_ad - nabla) - c_dif)*c_1*omega_c**2/lambda
       xA(5,3) = -V*c_dif + V*nabla_ad*(1._WP-U)
       xA(5,4) = -V*nabla_ad
       xA(5,5) = V*nabla*(4._WP - kappa_S) - (l_0 - 2._WP)
       xA(5,6) = -V*nabla/c_rad

       xA(6,1) = l_0*(l_0+1._WP)*(nabla_ad/nabla - 1._WP)*c_rad - V*c_eps_ad
       xA(6,2) = (V*c_eps_ad - lambda*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega_c**2)))*c_1*omega_c**2/lambda
       xA(6,3) = -(lambda*nabla_ad/nabla*c_rad - V*c_eps_ad)
       xA(6,4) = 0._WP
       if (x > 0._WP) then
          xA(6,5) = c_eps_S - lambda*c_rad/(nabla*V) - (0._WP,1._WP)*omega_c*c_thm
       else
          xA(6,5) = -HUGE(0._WP)
       endif
       xA(6,6) = -1._WP - l_0

    else

       xA(1,1) = V_g - 1._WP
       xA(1,2) = -V_g*c_1*omega_c**2
       xA(1,3) = -V_g
       xA(1,4) = 0._WP
       xA(1,5) = delta
       xA(1,6) = 0._WP
      
       xA(2,1) = 1._WP - As/(c_1*omega_c**2)
       xA(2,2) = As
       xA(2,3) = As/(c_1*omega_c**2)
       xA(2,4) = 0._WP
       xA(2,5) = delta/(c_1*omega_c**2)
       xA(2,6) = 0._WP
      
       xA(3,1) = 0._WP
       xA(3,2) = 0._WP
       xA(3,3) = 2._WP
       xA(3,4) = 1._WP
       xA(3,5) = 0._WP
       xA(3,6) = 0._WP
      
       xA(4,1) = -U*As
       xA(4,2) = -U*V_g*c_1*omega_c**2
       xA(4,3) = U*(As - 2._WP)
       xA(4,4) = 2._WP*(1._WP-U) + 1._WP
       xA(4,5) = U*delta
       xA(4,6) = 0._WP

       xA(5,1) = V*(nabla_ad*(U - c_1*omega_c**2) - 4._WP*(nabla_ad - nabla) + c_dif)
       xA(5,2) = V*(-c_dif)*c_1*omega_c**2
       xA(5,3) = -V*c_dif + V*nabla_ad*(1._WP-U)
       xA(5,4) = -V*nabla_ad
       xA(5,5) = V*nabla*(4._WP - kappa_S) - (l_0 - 2._WP)
       xA(5,6) = -V*nabla/c_rad

       xA(6,1) = -V*c_eps_ad
       xA(6,2) = V*c_eps_ad*c_1*omega_c**2
       xA(6,3) = -V*c_eps_ad
       xA(6,4) = 0._WP
       xA(6,5) = c_eps_S - (0._WP,1._WP)*omega_c*c_thm
       xA(6,6) = -1._WP

    endif

    ! Finish

    return

  end function xA_jcd_

!****

  function xA_lagp_ (this, x, omega) result (xA)

    class(nad_jacob_t), intent(in) :: this
    real(WP), intent(in)           :: x
    complex(WP), intent(in)        :: omega
    complex(WP)                    :: xA(this%n_e,this%n_e)

    real(WP)    :: V_2
    real(WP)    :: V
    real(WP)    :: Gamma_1
    real(WP)    :: U
    real(WP)    :: D
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
         
    ! Evaluate the log(x)-space RHS matrix (Lagrangian pressure
    ! perturbation formulation)

    ! Calculate coefficients

    V_2 = this%ml%V_2(x)
    V = V_2*x**2
    Gamma_1 = this%ml%Gamma_1(x)
    U = this%ml%U(x)
    D = this%ml%D(x)
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

    ! Set up the matrix

    xA(1,1) = lambda/(c_1*omega_c**2) - 1._WP - l_0
    xA(1,2) = lambda/(c_1*omega_c**2)/V_2 - x**2/Gamma_1
    xA(1,3) = lambda/(c_1*omega_c**2)
    xA(1,4) = 0._WP
    xA(1,5) = delta
    xA(1,6) = 0._WP

    xA(2,1) = -V_2*(lambda/(c_1*omega_c**2) - c_1*omega_c**2 + U - 4._WP)
    xA(2,2) = -lambda/(c_1*omega_c**2) + V - l_0
    xA(2,3) = -V_2*lambda/(c_1*omega_c**2)
    xA(2,4) = -V_2
    xA(2,5) = 0._WP
    xA(2,6) = 0._WP

    xA(3,1) = 0._WP
    xA(3,2) = 0._WP
    xA(3,3) = 3._WP - U - l_0
    xA(3,4) = 1._WP
    xA(3,5) = 0._WP
    xA(3,6) = 0._WP

    xA(4,1) = -U*D
    xA(4,2) = U*x**2/Gamma_1
    xA(4,3) = lambda
    xA(4,4) = -U - l_0 + 2._WP
    xA(4,5) = -U*delta
    xA(4,6) = 0._WP

    xA(5,1) = V*(lambda/(c_1*omega_c**2)*(nabla_ad-nabla) + nabla_ad*(U - c_1*omega_c**2) - 4._WP*(nabla_ad - nabla))
    xA(5,2) = x**2*(lambda/(c_1*omega_c**2)*(nabla_ad - nabla) - c_dif)
    xA(5,3) = V*(lambda/(c_1*omega_c**2)*(nabla_ad-nabla))
    xA(5,4) = V*nabla_ad
    xA(5,5) = V*nabla*(4._WP - kappa_S) - (l_0 - 2._WP)
    xA(5,6) = -V*nabla/c_rad

    xA(6,1) = lambda*c_rad*((3._WP + dc_rad)/(c_1*omega_c**2) - 1._WP)
    xA(6,2) = x**2*c_eps_ad - lambda*c_rad*(nabla_ad/nabla - (3._WP + dc_rad)/(c_1*omega_c**2))/V_2
    xA(6,3) = lambda*c_rad*(3._WP + dc_rad)/(c_1*omega_c**2)
    xA(6,4) = 0._WP
    if (x > 0._WP) then
       xA(6,5) = c_eps_S - lambda*c_rad/(nabla*V) - (0._WP,1._WP)*omega_c*c_thm
    else
       xA(6,5) = -HUGE(0._WP)
    endif
    xA(6,6) = -1._WP - l_0

    ! Finish

    return

  end function xA_lagp_

!****

  function T_ (this, x, omega, to_canon) result (T)

    class(nad_jacob_t), intent(in) :: this
    real(WP), intent(in)           :: x
    complex(WP), intent(in)        :: omega
    logical, intent(in)            :: to_canon
    complex(WP)                    :: T(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables to/from
    ! the canonical (DZEIM) formulation

    select case (this%vars)
    case (DZIEM_VARS)
       T = identity_matrix(this%n_e)
    case (JCD_VARS)
       T = this%T_jcd_(x, omega, to_canon)
    case (LAGP_VARS)
       T = this%T_lagp_(x, omega, to_canon)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function T_

!****

  function T_jcd_ (this, x, omega, to_canon) result (T)

    class(nad_jacob_t), intent(in) :: this
    real(WP), intent(in)           :: x
    complex(WP), intent(in)        :: omega
    logical, intent(in)            :: to_canon
    complex(WP)                    :: T(this%n_e,this%n_e)

    real(WP)    :: U
    real(WP)    :: c_1
    integer     :: l
    complex(WP) :: lambda
    complex(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! to/from the canonical (DZIEM) formulation

    ! Calculate coefficients

    U = this%ml%U(x)
    c_1 = this%ml%c_1(x)

    l = this%rt%mp%l
    lambda = this%rt%lambda(x, omega)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix

    if (to_canon) then

       if (l /= 0) then

          T(1,1) = 1._WP
          T(1,2) = 0._WP
          T(1,3) = 0._WP
          T(1,4) = 0._WP
          T(1,5) = 0._WP
          T(1,6) = 0._WP

          T(2,1) = 0._WP
          T(2,2) = c_1*omega_c**2/lambda
          T(2,3) = 0._WP
          T(2,4) = 0._WP
          T(2,5) = 0._WP
          T(2,6) = 0._WP

          T(3,1) = 0._WP
          T(3,2) = 0._WP
          T(3,3) = -1._WP
          T(3,4) = 0._WP
          T(3,5) = 0._WP
          T(3,6) = 0._WP

          T(4,1) = 0._WP
          T(4,2) = 0._WP
          T(4,3) = 1._WP - U
          T(4,4) = -1._WP
          T(4,5) = 0._WP
          T(4,6) = 0._WP
         
          T(5,1) = 0._WP
          T(5,2) = 0._WP
          T(5,3) = 0._WP
          T(5,4) = 0._WP
          T(5,5) = 1._WP
          T(5,6) = 0._WP

          T(6,1) = 0._WP
          T(6,2) = 0._WP
          T(6,3) = 0._WP
          T(6,4) = 0._WP
          T(6,5) = 0._WP
          T(6,6) = 1._WP

       else

          T(1,1) = 1._WP
          T(1,2) = 0._WP
          T(1,3) = 0._WP
          T(1,4) = 0._WP
          T(1,5) = 0._WP
          T(1,6) = 0._WP

          T(2,1) = 0._WP
          T(2,2) = c_1*omega_c**2
          T(2,3) = 0._WP
          T(2,4) = 0._WP
          T(2,5) = 0._WP
          T(2,6) = 0._WP

          T(3,1) = 0._WP
          T(3,2) = 0._WP
          T(3,3) = -1._WP
          T(3,4) = 0._WP
          T(3,5) = 0._WP
          T(3,6) = 0._WP

          T(4,1) = 0._WP
          T(4,2) = 0._WP
          T(4,3) = 1._WP - U
          T(4,4) = -1._WP
          T(4,5) = 0._WP
          T(4,6) = 0._WP
         
          T(5,1) = 0._WP
          T(5,2) = 0._WP
          T(5,3) = 0._WP
          T(5,4) = 0._WP
          T(5,5) = 1._WP
          T(5,6) = 0._WP

          T(6,1) = 0._WP
          T(6,2) = 0._WP
          T(6,3) = 0._WP
          T(6,4) = 0._WP
          T(6,5) = 0._WP
          T(6,6) = 1._WP

       end if

    else

       if (l /= 0) then

          T(1,1) = 1._WP
          T(1,2) = 0._WP
          T(1,3) = 0._WP
          T(1,4) = 0._WP
          T(1,5) = 0._WP
          T(1,6) = 0._WP

          T(2,1) = 0._WP
          T(2,2) = lambda/(c_1*omega_c**2)
          T(2,3) = 0._WP
          T(2,4) = 0._WP
          T(2,5) = 0._WP
          T(2,6) = 0._WP

          T(3,1) = 0._WP
          T(3,2) = 0._WP
          T(3,3) = -1._WP
          T(3,4) = 0._WP
          T(3,5) = 0._WP
          T(3,6) = 0._WP

          T(4,1) = 0._WP
          T(4,2) = 0._WP
          T(4,3) = -(1._WP - U)
          T(4,4) = -1._WP
          T(4,5) = 0._WP
          T(4,6) = 0._WP

          T(5,1) = 0._WP
          T(5,2) = 0._WP
          T(5,3) = 0._WP
          T(5,4) = 0._WP
          T(5,5) = 1._WP
          T(5,6) = 0._WP

          T(6,1) = 0._WP
          T(6,2) = 0._WP
          T(6,3) = 0._WP
          T(6,4) = 0._WP
          T(6,5) = 0._WP
          T(6,6) = 1._WP
          
       else

          T(1,1) = 1._WP
          T(1,2) = 0._WP
          T(1,3) = 0._WP
          T(1,4) = 0._WP
          T(1,5) = 0._WP
          T(1,6) = 0._WP

          T(2,1) = 0._WP
          T(2,2) = 1._WP/(c_1*omega_c**2)
          T(2,3) = 0._WP
          T(2,4) = 0._WP
          T(2,5) = 0._WP
          T(2,6) = 0._WP

          T(3,1) = 0._WP
          T(3,2) = 0._WP
          T(3,3) = -1._WP
          T(3,4) = 0._WP
          T(3,5) = 0._WP
          T(3,6) = 0._WP

          T(4,1) = 0._WP
          T(4,2) = 0._WP
          T(4,3) = -(1._WP - U)
          T(4,4) = -1._WP
          T(4,5) = 0._WP
          T(4,6) = 0._WP

          T(5,1) = 0._WP
          T(5,2) = 0._WP
          T(5,3) = 0._WP
          T(5,4) = 0._WP
          T(5,5) = 1._WP
          T(5,6) = 0._WP
          
          T(6,1) = 0._WP
          T(6,2) = 0._WP
          T(6,3) = 0._WP
          T(6,4) = 0._WP
          T(6,5) = 0._WP
          T(6,6) = 1._WP

       end if

    end if

    ! Finish

    return

  end function T_jcd_

!****

  function T_lagp_ (this, x, omega, to_canon) result (T)

    class(nad_jacob_t), intent(in) :: this
    real(WP), intent(in)           :: x
    complex(WP), intent(in)        :: omega
    logical, intent(in)            :: to_canon
    complex(WP)                    :: T(this%n_e,this%n_e)

    real(WP) :: V_2

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to/from the canonical (DZIEM) formulation

    ! Calculate coefficients

    V_2 = this%ml%V_2(x)

    ! Set up the matrix

    if (to_canon) then

       T(1,1) = 1._WP
       T(1,2) = 0._WP
       T(1,3) = 0._WP
       T(1,4) = 0._WP
       T(1,5) = 0._WP
       T(1,6) = 0._WP

       T(2,1) = 1._WP
       T(2,2) = 1._WP/V_2
       T(2,3) = 1._WP
       T(2,4) = 0._WP
       T(2,5) = 0._WP
       T(2,6) = 0._WP

       T(3,1) = 0._WP
       T(3,2) = 0._WP
       T(3,3) = 1._WP
       T(3,4) = 0._WP
       T(3,5) = 0._WP
       T(3,6) = 0._WP

       T(4,1) = 0._WP
       T(4,2) = 0._WP
       T(4,3) = 0._WP
       T(4,4) = 1._WP
       T(4,5) = 0._WP
       T(4,6) = 0._WP

       T(5,1) = 0._WP
       T(5,2) = 0._WP
       T(5,3) = 0._WP
       T(5,4) = 0._WP
       T(5,5) = 1._WP
       T(5,6) = 0._WP

       T(6,1) = 0._WP
       T(6,2) = 0._WP
       T(6,3) = 0._WP
       T(6,4) = 0._WP
       T(6,5) = 0._WP
       T(6,6) = 1._WP

    else

       T(1,1) = 1._WP
       T(1,2) = 0._WP
       T(1,3) = 0._WP
       T(1,4) = 0._WP
       T(1,5) = 0._WP
       T(1,6) = 0._WP

       T(2,1) = -V_2
       T(2,2) = V_2
       T(2,3) = -V_2
       T(2,4) = 0._WP
       T(2,5) = 0._WP
       T(2,6) = 0._WP

       T(3,1) = 0._WP
       T(3,2) = 0._WP
       T(3,3) = 1._WP
       T(3,4) = 0._WP
       T(3,5) = 0._WP
       T(3,6) = 0._WP

       T(4,1) = 0._WP
       T(4,2) = 0._WP
       T(4,3) = 0._WP
       T(4,4) = 1._WP
       T(4,5) = 0._WP
       T(4,6) = 0._WP

       T(5,1) = 0._WP
       T(5,2) = 0._WP
       T(5,3) = 0._WP
       T(5,4) = 0._WP
       T(5,5) = 1._WP
       T(5,6) = 0._WP

       T(6,1) = 0._WP
       T(6,2) = 0._WP
       T(6,3) = 0._WP
       T(6,4) = 0._WP
       T(6,5) = 0._WP
       T(6,6) = 1._WP

    end if

    ! Finish

    return

  end function T_lagp_

end module gyre_nad_jacob
