! Module   : gyre_ad_eqns
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

module gyre_ad_eqns

  ! Uses

  use core_kinds

  use gyre_eqns
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
  integer, parameter :: MIX_VARS = 3
  integer, parameter :: LAGP_VARS = 4

  ! Derived-type definitions

  type, extends (r_eqns_t) :: ad_eqns_t
     private
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     integer                     :: vars
   contains
     private
     procedure, public :: A => A_
     procedure, public :: xA => xA_
     procedure         :: xA_dziem_
     procedure         :: xA_jcd_
     procedure         :: xA_mix_
     procedure         :: xA_lagp_
     procedure, public :: T => T_
     procedure         :: T_jcd_
     procedure         :: T_mix_
     procedure         :: T_lagp_
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

  function ad_eqns_t_ (ml, rt, op) result (eq)

    class(model_t), pointer, intent(in) :: ml
    class(r_rot_t), intent(in)          :: rt
    type(osc_par_t), intent(in)         :: op
    type(ad_eqns_t)                     :: eq

    ! Construct the ad_eqns_t

    eq%ml => ml
    allocate(eq%rt, SOURCE=rt)

    select case (op%variables_set)
    case ('DZIEM')
       eq%vars = DZIEM_VARS
    case ('JCD')
       eq%vars = JCD_VARS
    case ('MIX')
       eq%vars = MIX_VARS
    case ('LAGP')
       eq%vars = LAGP_VARS
    case default
       $ABORT(Invalid variables_set)
    end select

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
    
    ! Evaluate the log(x)-space RHS matrix (=x*A)

    select case (this%vars)
    case (DZIEM_VARS)
       xA = this%xA_dziem_(x, omega)
    case (JCD_VARS)
       xA = this%xA_jcd_(x, omega)
    case (MIX_VARS)
       xA = this%xA_mix_(x, omega)
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
    
    ! Evaluate the log(x)-space RHS matrix ([Dzi1971] formulation)

    ! Calculate coefficients

    V_g = this%ml%V_2(x)*x**2/this%ml%Gamma_1(x)
    U = this%ml%U(x)
    As = this%ml%As(x)
    c_1 = this%ml%c_1(x)

    lambda = this%rt%lambda(x, omega)
    l_0 = this%rt%l_0(omega)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix

    xA(1,1) = V_g - 1._WP - l_0
    xA(1,2) = lambda/(c_1*omega_c**2) - V_g
    xA(1,3) = V_g
    xA(1,4) = 0._WP
      
    xA(2,1) = c_1*omega_c**2 - As
    xA(2,2) = As - U + 3._WP - l_0
    xA(2,3) = -As
    xA(2,4) = 0._WP
      
    xA(3,1) = 0._WP
    xA(3,2) = 0._WP
    xA(3,3) = 3._WP - U - l_0
    xA(3,4) = 1._WP
      
    xA(4,1) = U*As
    xA(4,2) = U*V_g
    xA(4,3) = lambda - U*V_g
    xA(4,4) = -U - l_0 + 2._WP

    ! Finish

    return

  end function xA_dziem_

!****

  function xA_jcd_ (this, x, omega) result (xA)

    class(ad_eqns_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: xA(this%n_e,this%n_e)

    real(WP) :: V_g
    real(WP) :: U
    real(WP) :: As
    real(WP) :: c_1
    integer  :: l
    real(WP) :: lambda
    real(WP) :: l_0
    real(WP) :: omega_c
    
    ! Evaluate the log(x)-space RHS matrix ([Chr2008]
    ! formulation)

    ! Calculate coefficients

    V_g = this%ml%V_2(x)*x**2/this%ml%Gamma_1(x)
    U = this%ml%U(x)
    As = this%ml%As(x)
    c_1 = this%ml%c_1(x)

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
      
       xA(2,1) = lambda - As*lambda/(c_1*omega_c**2)
       xA(2,2) = As - l_0
       xA(2,3) = As*lambda/(c_1*omega_c**2)
       xA(2,4) = 0._WP
      
       xA(3,1) = 0._WP
       xA(3,2) = 0._WP
       xA(3,3) = 2._WP - l_0
       xA(3,4) = 1._WP
      
       xA(4,1) = -U*As
       xA(4,2) = -U*V_g*c_1*omega_c**2/lambda
       xA(4,3) = lambda + U*(As - 2._WP)
       xA(4,4) = 2._WP*(1._WP-U) - (l_0 - 1._WP)

    else

       xA(1,1) = V_g - 1._WP
       xA(1,2) = -V_g*c_1*omega_c**2
       xA(1,3) = -V_g
       xA(1,4) = 0._WP
      
       xA(2,1) = 1._WP - As/(c_1*omega_c**2)
       xA(2,2) = As
       xA(2,3) = As/(c_1*omega_c**2)
       xA(2,4) = 0._WP
      
       xA(3,1) = 0._WP
       xA(3,2) = 0._WP
       xA(3,3) = 2._WP
       xA(3,4) = 1._WP
      
       xA(4,1) = -U*As
       xA(4,2) = -U*V_g*c_1*omega_c**2
       xA(4,3) = U*(As - 2._WP)
       xA(4,4) = 2._WP*(1._WP-U) + 1._WP

    endif

    ! Finish

    return

  end function xA_jcd_

!****

  function xA_mix_ (this, x, omega) result (xA)

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
    
    ! Evaluate the log(x)-space RHS matrix (mixed formulation)

    ! Calculate coefficients

    V_g = this%ml%V_2(x)*x**2/this%ml%Gamma_1(x)
    U = this%ml%U(x)
    As = this%ml%As(x)
    c_1 = this%ml%c_1(x)

    lambda = this%rt%lambda(x, omega)
    l_0 = this%rt%l_0(omega)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix

    xA(1,1) = V_g - 1._WP - l_0
    xA(1,2) = lambda/(c_1*omega_c**2) - V_g
    xA(1,3) = -V_g
    xA(1,4) = 0._WP
      
    xA(2,1) = c_1*omega_c**2 - As
    xA(2,2) = As - U + 3._WP - l_0
    xA(2,3) = As
    xA(2,4) = 0._WP
      
    xA(3,1) = 0._WP
    xA(3,2) = 0._WP
    xA(3,3) = 2._WP - l_0
    xA(3,4) = 1._WP
      
    xA(4,1) = -U*As
    xA(4,2) = -U*V_g
    xA(4,3) = lambda + U*(As - 2._WP)
    xA(4,4) = 2._WP*(1._WP-U) - (l_0 - 1._WP)

    ! Finish

    return

  end function xA_mix_

!****

  function xA_lagp_ (this, x, omega) result (xA)

    class(ad_eqns_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    real(WP)                     :: xA(this%n_e,this%n_e)

    real(WP) :: V_2
    real(WP) :: V
    real(WP) :: Gamma_1
    real(WP) :: U
    real(WP) :: D
    real(WP) :: c_1
    real(WP) :: lambda
    real(WP) :: l_0
    real(WP) :: omega_c
    
    ! Evaluate the log(x)-space RHS matrix (Lagrangian pressure
    ! perturbation formulation)

    ! Calculate coefficients

    V_2 = this%ml%V_2(x)
    V = V_2*x**2
    Gamma_1 = this%ml%Gamma_1(x)
    U = this%ml%U(x)
    D = this%ml%D(x)
    c_1 = this%ml%c_1(x)

    lambda = this%rt%lambda(x, omega)
    l_0 = this%rt%l_0(omega)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix

    xA(1,1) = lambda/(c_1*omega_c**2) - 1._WP - l_0
    xA(1,2) = lambda/(c_1*omega_c**2)/V_2 - x**2/Gamma_1
    xA(1,3) = lambda/(c_1*omega_c**2)
    xA(1,4) = 0._WP
      
    xA(2,1) = -V_2*(lambda/(c_1*omega_c**2) - c_1*omega_c**2 + U - 4._WP)
    xA(2,2) = -lambda/(c_1*omega_c**2) + V - l_0
    xA(2,3) = -V_2*lambda/(c_1*omega_c**2)
    xA(2,4) = -V_2
      
    xA(3,1) = 0._WP
    xA(3,2) = 0._WP
    xA(3,3) = 3._WP - U - l_0
    xA(3,4) = 1._WP
      
    xA(4,1) = -U*D
    xA(4,2) = U*x**2/Gamma_1
    xA(4,3) = lambda
    xA(4,4) = -U - l_0 + 2._WP

    ! Finish

    return

  end function xA_lagp_
  
!****

  function T_ (this, x, omega, to_canon) result (T)

    class(ad_eqns_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    logical, intent(in)          :: to_canon
    real(WP)                     :: T(this%n_e,this%n_e)

    ! Evaluate the transformation matrix to convert variables to/from
    ! the canonical (DZEIM) formulation

    select case (this%vars)
    case (DZIEM_VARS)
       T = identity_matrix(this%n_e)
    case (JCD_VARS)
       T = this%T_jcd_(x, omega, to_canon)
    case (MIX_VARS)
       T = this%T_mix_(x, omega, to_canon)
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

    class(ad_eqns_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    logical, intent(in)          :: to_canon
    real(WP)                     :: T(this%n_e,this%n_e)

    real(WP) :: U
    real(WP) :: c_1
    integer  :: l
    real(WP) :: lambda
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! to/from the canonical (DZEIM) formulation

    ! Calculate coefficients

    U = this%ml%U(x)
    c_1 = this%ml%c_1(x)

    l = this%rt%mp%l
    lambda = this%rt%lambda(x, omega)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix
      
    if (to_canon) then

       if (l /= 0._WP) then

          T(1,1) = 1._WP
          T(1,2) = 0._WP
          T(1,3) = 0._WP
          T(1,4) = 0._WP

          T(2,1) = 0._WP
          T(2,2) = c_1*omega_c**2/lambda
          T(2,3) = 0._WP
          T(2,4) = 0._WP

          T(3,1) = 0._WP
          T(3,2) = 0._WP
          T(3,3) = -1._WP
          T(3,4) = 0._WP

          T(4,1) = 0._WP
          T(4,2) = 0._WP
          T(4,3) = 1._WP - U
          T(4,4) = -1._WP

       else

          T(1,1) = 1._WP
          T(1,2) = 0._WP
          T(1,3) = 0._WP
          T(1,4) = 0._WP

          T(2,1) = 0._WP
          T(2,2) = c_1*omega_c**2
          T(2,3) = 0._WP
          T(2,4) = 0._WP

          T(3,1) = 0._WP
          T(3,2) = 0._WP
          T(3,3) = -1._WP
          T(3,4) = 0._WP

          T(4,1) = 0._WP
          T(4,2) = 0._WP
          T(4,3) = 1._WP - U
          T(4,4) = -1._WP

       endif

    else

       if (l /= 0) then

          T(1,1) = 1._WP
          T(1,2) = 0._WP
          T(1,3) = 0._WP
          T(1,4) = 0._WP
          
          T(2,1) = 0._WP
          T(2,2) = lambda/(c_1*omega_c**2)
          T(2,3) = 0._WP
          T(2,4) = 0._WP
          
          T(3,1) = 0._WP
          T(3,2) = 0._WP
          T(3,3) = -1._WP
          T(3,4) = 0._WP

          T(4,1) = 0._WP
          T(4,2) = 0._WP
          T(4,3) = -(1._WP - U)
          T(4,4) = -1._WP

       else

          T(1,1) = 1._WP
          T(1,2) = 0._WP
          T(1,3) = 0._WP
          T(1,4) = 0._WP

          T(2,1) = 0._WP
          T(2,2) = 1._WP/(c_1*omega_c**2)
          T(2,3) = 0._WP
          T(2,4) = 0._WP

          T(3,1) = 0._WP
          T(3,2) = 0._WP
          T(3,3) = -1._WP
          T(3,4) = 0._WP

          T(4,1) = 0._WP
          T(4,2) = 0._WP
          T(4,3) = -(1._WP - U)
          T(4,4) = -1._WP

       endif

    end if

    ! Finish

    return

  end function T_jcd_

!****

  function T_mix_ (this, x, omega, to_canon) result (T)

    class(ad_eqns_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    logical, intent(in)          :: to_canon
    real(WP)                     :: T(this%n_e,this%n_e)

    real(WP) :: U

    ! Evaluate the transformation matrix to convert MIX variables
    ! to/from the canonical (DZEIM) formulation

    ! Calculate coefficients

    U = this%ml%U(x)

    ! Set up the matrix

    if (to_canon) then

       T(1,1) = 1._WP
       T(1,2) = 0._WP
       T(1,3) = 0._WP
       T(1,4) = 0._WP

       T(2,1) = 0._WP
       T(2,2) = 1._WP
       T(2,3) = 0._WP
       T(2,4) = 0._WP

       T(3,1) = 0._WP
       T(3,2) = 0._WP
       T(3,3) = -1._WP
       T(3,4) = 0._WP

       T(4,1) = 0._WP
       T(4,2) = 0._WP
       T(4,3) = 1._WP - U
       T(4,4) = -1._WP

    else

       T(1,1) = 1._WP
       T(1,2) = 0._WP
       T(1,3) = 0._WP
       T(1,4) = 0._WP

       T(2,1) = 0._WP
       T(2,2) = 1._WP
       T(2,3) = 0._WP
       T(2,4) = 0._WP

       T(3,1) = 0._WP
       T(3,2) = 0._WP
       T(3,3) = -1._WP
       T(3,4) = 0._WP

       T(4,1) = 0._WP
       T(4,2) = 0._WP
       T(4,3) = -(1._WP - U)
       T(4,4) = -1._WP

    end if

    ! Finish

    return

  end function T_mix_

!****

  function T_lagp_ (this, x, omega, to_canon) result (T)

    class(ad_eqns_t), intent(in) :: this
    real(WP), intent(in)         :: x
    real(WP), intent(in)         :: omega
    logical, intent(in)          :: to_canon
    real(WP)                     :: T(this%n_e,this%n_e)

    real(WP) :: V
    real(WP) :: V_2

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to/from the canonical (DZEIM) formulation

    ! Calculate coefficients

    V_2 = this%ml%V_2(x)

    ! Set up the matrix

    if (to_canon) then

       T(1,1) = 1._WP
       T(1,2) = 0._WP
       T(1,3) = 0._WP
       T(1,4) = 0._WP

       T(2,1) = 1._WP
       T(2,2) = 1._WP/V_2
       T(2,3) = 1._WP
       T(2,4) = 0._WP

       T(3,1) = 0._WP
       T(3,2) = 0._WP
       T(3,3) = 1._WP
       T(3,4) = 0._WP

       T(4,1) = 0._WP
       T(4,2) = 0._WP
       T(4,3) = 0._WP
       T(4,4) = 1._WP

    else

       T(1,1) = 1._WP
       T(1,2) = 0._WP
       T(1,3) = 0._WP
       T(1,4) = 0._WP

       T(2,1) = -V_2
       T(2,2) = V_2
       T(2,3) = -V_2
       T(2,4) = 0._WP

       T(3,1) = 0._WP
       T(3,2) = 0._WP
       T(3,3) = 1._WP
       T(3,4) = 0._WP

       T(4,1) = 0._WP
       T(4,2) = 0._WP
       T(4,3) = 0._WP
       T(4,4) = 1._WP

    end if

    ! Finish

    return

  end function T_lagp_

end module gyre_ad_eqns
