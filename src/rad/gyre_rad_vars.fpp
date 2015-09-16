! Module   : gyre_rad_vars
! Purpose  : variables transformations (adiabatic radial)
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

module gyre_rad_vars

  ! Uses

  use core_kinds

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

  type :: rad_vars_t
     private
     class(model_t), pointer     :: ml => null()
     class(r_rot_t), allocatable :: rt
     integer                     :: vars
   contains
     private
     procedure, public :: S => S_
     procedure         :: S_jcd_
     procedure         :: S_lagp_
     procedure, public :: T => T_
     procedure         :: T_jcd_
     procedure         :: T_lagp_
     procedure, public :: dT => dT_
     procedure         :: dT_jcd_
     procedure         :: dT_lagp_
  end type rad_vars_t

  ! Interfaces

  interface rad_vars_t
     module procedure rad_vars_t_
  end interface rad_vars_t

  ! Access specifiers

  private

  public :: rad_vars_t

  ! Procedures

contains

  function rad_vars_t_ (ml, rt, op) result (vr)

    class(model_t), pointer, intent(in) :: ml
    class(r_rot_t), intent(in)          :: rt
    type(osc_par_t), intent(in)         :: op
    type(rad_vars_t)                    :: vr

    ! Construct the ad_vars_t

    vr%ml => ml
    allocate(vr%rt, SOURCE=rt)

    select case (op%variables_set)
    case ('DZIEM')
       vr%vars = DZIEM_VARS
    case ('JCD')
       vr%vars = JCD_VARS
    case ('MIX')
       vr%vars = MIX_VARS
    case ('LAGP')
       vr%vars = LAGP_VARS
    case default
       $ABORT(Invalid variables_set)
    end select

    ! Finish

    return

  end function rad_vars_t_

!****

  function S_ (this, x, omega) result (S)

    class(rad_vars_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: S(2,2)

    ! Evaluate the transformation matrix to convert variables from
    ! the canonical form

    select case (this%vars)
    case (DZIEM_VARS)
       S = identity_matrix(2)
    case (JCD_VARS)
       S = this%S_jcd_(x, omega)
    case (MIX_VARS)
       S = identity_matrix(2)
    case (LAGP_VARS)
       S = this%S_lagp_(x, omega)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function S_

!****

  function S_jcd_ (this, x, omega) result (S)

    class(rad_vars_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: S(2,2)

    real(WP) :: U
    real(WP) :: c_1
    integer  :: l
    real(WP) :: lambda
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! from the canonical form

    ! Calculate coefficients

    U = this%ml%U(x)
    c_1 = this%ml%c_1(x)

    l = this%rt%mp%l
    lambda = this%rt%lambda(x, omega)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix

    S(1,1) = 1._WP
    S(1,2) = 0._WP

    S(2,1) = 0._WP
    S(2,2) = 1._WP/(c_1*omega_c**2)

    ! Finish

    return

  end function S_jcd_

!****

  function S_lagp_ (this, x, omega) result (S)

    class(rad_vars_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: S(2,2)

    real(WP) :: V_2

    ! Evaluate the transformation matrix to convert LAGP variables
    ! from the canonical form

    ! Calculate coefficients

    V_2 = this%ml%V_2(x)

    ! Set up the matrix

    S(1,1) = 1._WP
    S(1,2) = 0._WP

    S(2,1) = -V_2
    S(2,2) = V_2

    ! Finish

    return

  end function S_lagp_

!****

  function T_ (this, x, omega) result (T)

    class(rad_vars_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: T(2,2)

    ! Evaluate the transformation matrix to convert variables to
    ! canonical form

    select case (this%vars)
    case (DZIEM_VARS)
       T = identity_matrix(2)
    case (JCD_VARS)
       T = this%T_jcd_(x, omega)
    case (MIX_VARS)
       T = identity_matrix(2)
    case (LAGP_VARS)
       T = this%T_lagp_(x, omega)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function T_

!****

  function T_jcd_ (this, x, omega) result (T)

    class(rad_vars_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: T(2,2)

    real(WP) :: c_1
    real(WP) :: omega_c

    ! Evaluate the transformation matrix to convert JCD variables
    ! to the canonical form

    ! Calculate coefficients

    c_1 = this%ml%c_1(x)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix
      
    T(1,1) = 1._WP
    T(1,2) = 0._WP

    T(2,1) = 0._WP
    T(2,2) = c_1*omega_c**2

    ! Finish

    return

  end function T_jcd_

!****

  function T_lagp_ (this, x, omega) result (T)

    class(rad_vars_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: T(2,2)

    real(WP) :: V_2

    ! Evaluate the transformation matrix to convert LAGP variables
    ! to the canonical form

    ! Calculate coefficients

    V_2 = this%ml%V_2(x)

    ! Set up the matrix

    T(1,1) = 1._WP
    T(1,2) = 0._WP

    T(2,1) = 1._WP
    T(2,2) = 1._WP/V_2

    ! Finish

    return

  end function T_lagp_

!****

  function dT_ (this, x, omega) result (dT)

    class(rad_vars_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: dT(2,2)

    ! Evaluate the derivative x dT/dx of the transformation matrix T

    select case (this%vars)
    case (DZIEM_VARS)
       dT = 0._WP
    case (JCD_VARS)
       dT = this%dT_jcd_(x, omega)
    case (MIX_VARS)
       dT = 0._WP
    case (LAGP_VARS)
       dT = this%dT_lagp_(x, omega)
    case default
       $ABORT(Invalid vars)
    end select

    ! Finish

    return

  end function dT_

!****

  function dT_jcd_ (this, x, omega) result (dT)

    class(rad_vars_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: dT(2,2)

    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U
    real(WP) :: c_1
    real(WP) :: omega_c

    ! Evaluate the derivative x dT/dx of the JCD-variables
    ! transformation matrix T

    ! Calculate coefficients

    V_g = this%ml%V_2(x)*x**2/this%ml%Gamma_1(x)
    As = this%ml%As(x) 
    U = this%ml%U(x)
    c_1 = this%ml%c_1(x)

    omega_c = this%rt%omega_c(x, omega)

    ! Set up the matrix (nb: the derivative of omega_c is neglected;
    ! this is incorrect when rotation is non-zero)
      
    dT(1,1) = 0._WP
    dT(1,2) = 0._WP

    dT(2,1) = 0._WP
    dT(2,2) = c_1*(3._WP - U)*omega_c**2

    ! Finish

    return

  end function dT_jcd_

!****

  function dT_lagp_ (this, x, omega) result (dT)

    class(rad_vars_t), intent(in) :: this
    real(WP), intent(in)          :: x
    real(WP), intent(in)          :: omega
    real(WP)                      :: dT(2,2)

    real(WP) :: V_2
    real(WP) :: V
    real(WP) :: V_g
    real(WP) :: As
    real(WP) :: U

    ! Evaluate the derivative x dT/dx of the LAGP-variables
    ! transformation matrix T

    ! Calculate coefficients

    V_2 = this%ml%V_2(x)
    V = V_2*x**2
    V_g = V/this%ml%Gamma_1(x)
    As = this%ml%As(x)
    U = this%ml%U(x)

    ! Set up the matrix

    dT(1,1) = 0._WP
    dT(1,2) = 0._WP

    dT(2,1) = 0._WP
    dT(2,2) = -(-V_g - As + U + V - 3)/V_2

    ! Finish

    return

  end function dT_lagp_
  
end module gyre_rad_vars
